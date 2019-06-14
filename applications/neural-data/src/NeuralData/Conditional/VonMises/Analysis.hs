{-# LANGUAGE
    TypeApplications,
    DeriveGeneric,
    OverloadedStrings
    #-}

module NeuralData.Conditional.VonMises.Analysis
    (
    -- * Analyses
      analyzePopulationCurves
    , preferredStimulusHistogram
    , precisionsHistogram
    , gainsHistograms
    -- ** Cross-Validation
    , CrossValidationStats (CrossValidationStats)
    , crossValidateMixtureLikelihood
    , crossValidateEmpiricalCovariances
    -- ** IO
    , runPopulationParameterAnalyses
    -- * Miscellaneous
    , kFoldDataset
    , kFoldConditionalDataset
    , dataCovarianceFunction
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Foreign.Storable
import qualified Goal.Core.Vector.Storable as S

import NeuralData.Conditional.VonMises
import NeuralData.Conditional.VonMises.Training

import NeuralData
import Data.String

import qualified Data.List as L
import qualified Data.Map as M
import Data.Tuple

--- Types ---


--- Inference ---


dataCovarianceFunction :: KnownNat k => [(Response k,Double)] -> Double -> Source # MultivariateNormal k
{-# INLINE dataCovarianceFunction #-}
dataCovarianceFunction zxs x =
    let mpnrms = M.fromList $ swap <$> dataCovariances zxs
        kys = M.keys mpnrms
        ky = last $ last kys : takeWhile (<= x) kys
     in mpnrms M.! ky


--- Analyses ---


runPopulationParameterAnalyses
    :: forall k m . (KnownNat k, KnownNat m)
    => Experiment
    -> String
    -> [Double]
    -> Int
    -> String
    -> String
    -> Maybe (Natural # LogNormal)
    -> Maybe [Natural # LogNormal]
    -> [(Response k, Double)]
    -> Natural #> ConditionalMixture (Neurons k) m VonMises
    -> IO ()
runPopulationParameterAnalyses expmnt dst xsmps nbns rltv ttl mprcsdst mgndsts zxs mlkl = do

    --let mlkl = sortVonMisesMixturePopulationEncoder mlkl0
    let msbexp = Just $ Analysis ("population-parameters/" ++ ttl) dst
        (tcs,gps,ccrvs,ctcrvs) = analyzePopulationCurves xsmps mlkl

    goalExportNamed True expmnt msbexp tcs
    goalExportNamed False expmnt msbexp gps
    goalExportNamed False expmnt msbexp ccrvs
    goalExportNamed False expmnt msbexp ctcrvs

    let tcgpi = rltv ++ "tuning-curves.gpi"
        gpgpi = rltv ++ "gain-profiles.gpi"
        ccgpi = rltv ++ "conjugacy-curves.gpi"
        ctgpi = rltv ++ "category-dependence.gpi"

    mapM_ (runGnuplot expmnt msbexp defaultGnuplotOptions) [tcgpi,gpgpi,ccgpi,ctgpi]

    let (pfshst,pfsft,_) = preferredStimulusHistogram nbns mlkl Nothing

    goalExportNamed False expmnt msbexp pfshst
    goalExportNamed False expmnt msbexp pfsft

    let (prcshst,prcsft,_) = precisionsHistogram nbns mlkl mprcsdst

    goalExportNamed False expmnt msbexp prcshst
    goalExportNamed False expmnt msbexp prcsft

    let (gnhsts,gnfts,_) = gainsHistograms nbns mlkl $ map breakPoint <$> mgndsts

    mapM_ (goalExportNamed False expmnt msbexp) gnhsts
    mapM_ (goalExportNamed False expmnt msbexp) gnfts

    let phgpi = rltv ++ "population-histogram.gpi"
    runGnuplot expmnt msbexp defaultGnuplotOptions phgpi

    let crsbexp = Just $ Analysis ("noise-correlations-vs-data/" ++ ttl) dst
        mvn = dataCovarianceFunction zxs

    let (mtxln:mtxlns) = do
            x <- xsmps
            let mdlcrs = mixturePopulationNoiseCorrelations (mlkl >.>* x)
                mvncrs = multivariateNormalCorrelations $ mvn x
                mlticrs = S.combineTriangles (S.replicate 1) mdlcrs mvncrs
            return $ S.toList <$> S.toList (S.toRows mlticrs)

    goalExport True expmnt crsbexp mtxln
    mapM_ (goalExport False expmnt crsbexp) mtxlns

    let aniopts = defaultGnuplotOptions { whetherPNG = False, whetherAnimate = True }
        crgpi = rltv ++ "noise-correlations.gpi"

    runGnuplot expmnt crsbexp aniopts crgpi

    let expnm = experimentName expmnt

    mtrulkl0 <- getMixtureLikelihood expnm dst

    unless (null mtrulkl0) $ do

        let (_,m',cs) = fromJust mtrulkl0

        case someNatVal m'
            of SomeNat (Proxy :: Proxy m') -> do

                let trulkl :: Natural #> ConditionalMixture (Neurons k) m' VonMises
                    trulkl = strengthenMixtureLikelihood cs

                    --trulkl = sortVonMisesMixturePopulationEncoder trulkl0

                    trucrsbexp = Just $ Analysis ("noise-correlations-vs-true/" ++ ttl) dst

                let (mtxln':mtxlns') = do
                        x <- xsmps
                        let mdlcrs = mixturePopulationNoiseCorrelations (mlkl >.>* x)
                            trucrs = mixturePopulationNoiseCorrelations (trulkl >.>* x)
                            mlticrs = S.combineTriangles (S.replicate 1) mdlcrs trucrs
                        return $ S.toList <$> S.toList (S.toRows mlticrs)

                goalExport True expmnt trucrsbexp mtxln'
                mapM_ (goalExport False expmnt trucrsbexp) mtxlns'

                runGnuplot expmnt trucrsbexp aniopts crgpi


--- Analysis ---


analyzePopulationCurves
    :: forall m k . (KnownNat m, KnownNat k)
    => Sample VonMises
    -> (Natural #> ConditionalMixture (Neurons k) m VonMises)
    -> ([TuningCurves k],[GainProfiles (m+1)],[ConjugacyCurves (m+1)],[CategoryDependence (m+1)])
{-# INLINE analyzePopulationCurves #-}
analyzePopulationCurves smps mlkl =
    let (_,ngnss,tcrws) = splitVonMisesMixturePopulationEncoder mlkl
        nrmlkl = joinVonMisesPopulationEncoder (Left 1) tcrws
        tcs = zipWith TuningCurves smps $ coordinates . toSource <$> nrmlkl >$>* smps
        mus = head . listCoordinates . toSource <$> S.toList tcrws
        gps = map (uncurry GainProfiles) . L.sortOn fst . zip mus
            . S.toList . S.toRows . S.fromColumns $ S.map (coordinates . toSource) ngnss
        mgnxs = mlkl >$>* smps
        cnjs = potential <$> mgnxs
        (cts,stcs) = unzip $ do
            (rts,ct) <- splitMixture <$> mgnxs
            let sct = coordinates $ toSource ct
            return (S.cons (1 - S.sum sct) sct, S.map potential rts)
        ccrvs = L.zipWith3 ConjugacyCurves smps stcs cnjs
        ctcrvs = L.zipWith CategoryDependence smps cts
     in (tcs,gps,ccrvs,ctcrvs)

preferredStimulusHistogram
    :: (KnownNat k, KnownNat m)
    => Int
    -> (Natural #> ConditionalMixture (Neurons k) m VonMises)
    -> Maybe (Natural # VonMises)
    -> ([PreferredStimuli], [ParameterDistributionFit], Natural # VonMises)
{-# INLINE preferredStimulusHistogram #-}
preferredStimulusHistogram nbns mlkl mtrudns =
    let (_,_,nxs) = splitVonMisesMixturePopulationEncoder mlkl
        mus =  head . listCoordinates . toSource <$> S.toList nxs
        (bns,[cnts],[wghts]) = histograms nbns Nothing [mus]
        vm0 = Point $ S.doubleton 0.01 0.01
        vm :: Natural # VonMises
        vm = vanillaGradientSequence (logLikelihoodDifferential mus) 0.1 defaultAdamPursuit vm0 !! 500
        xs = range mnx mxx 1000
        dnss = density vm <$> xs
        mdnss = [ density <$> mtrudns <*> Just x | x <- xs ]
     in (zipWith3 PreferredStimuli bns cnts wghts, zipWith3 ParameterDistributionFit xs dnss mdnss, vm)

logNormalHistogram
    :: Int
    -> Maybe (Natural # LogNormal)
    -> [Double]
    -> ( ([Double],[Int],[Double])
       , ([Double],[Double],[Maybe Double])
       , Natural # LogNormal )
{-# INLINE logNormalHistogram #-}
logNormalHistogram nbns mtrudns prms =
    let (bns,[cnts],[wghts]) = histograms nbns Nothing [prms]
        dx = head (tail bns) - head bns
        lgnrm :: Natural # LogNormal
        lgnrm = mle $ filter (/= 0) prms
        xs = range 0 (last bns + dx/2) 1000
        dnss = density lgnrm <$> xs
        mdnss = [ density <$> mtrudns <*> Just x | x <- xs ]
     in ((bns,cnts,wghts),(xs,dnss,mdnss),lgnrm)

precisionsHistogram
    :: (KnownNat k, KnownNat m)
    => Int
    -> (Natural #> ConditionalMixture (Neurons k) m VonMises)
    -> Maybe (Natural # LogNormal)
    -> ([Precisions], [ParameterDistributionFit], Natural # LogNormal)
{-# INLINE precisionsHistogram #-}
precisionsHistogram nbns mlkl mtrudns =
    let (_,_,nxs) = splitVonMisesMixturePopulationEncoder mlkl
        rhos = head . tail . listCoordinates . toSource <$> S.toList nxs
        ((bns,cnts,wghts),(xs,dnss,mdnss),lgnrm) = logNormalHistogram nbns mtrudns rhos
     in (zipWith3 Precisions bns cnts wghts, zipWith3 ParameterDistributionFit xs dnss mdnss, lgnrm)

gainsHistograms
    :: (KnownNat k, KnownNat m)
    => Int
    -> Natural #> ConditionalMixture (Neurons k) m VonMises
    -> Maybe [Natural # LogNormal]
    -> ([[Gains]], [[ParameterDistributionFit]], [Natural # LogNormal])
{-# INLINE gainsHistograms #-}
gainsHistograms nbns mlkl mtrudnss0 =
    let (_,ngnss,_) = splitVonMisesMixturePopulationEncoder mlkl
        gnss = listCoordinates . toSource <$> S.toList ngnss
        mtrudnss = maybe (repeat Nothing) (map Just) mtrudnss0
     in unzip3 $ do
            (mtrudns,gns) <- zip mtrudnss gnss
            let ((bns,cnts,wghts),(xs,dnss,mdnss),lgnrm) = logNormalHistogram nbns mtrudns gns
            return ( zipWith3 Gains bns cnts wghts
                   , zipWith3 ParameterDistributionFit xs dnss mdnss
                   , lgnrm )

crossValidateMixtureLikelihood
    :: forall k n r . (KnownNat k, KnownNat n)
    => Initialization k n
    -> Optimization
    -> Int -- ^ Number of parallel fits
    -> Int -- ^ Number of epochs
    -> [([(Response k, Double)],[(Response k, Double)])] -- ^ K-folded data
    -> Random r (Maybe Double)
{-# INLINE crossValidateMixtureLikelihood #-}
crossValidateMixtureLikelihood intl optm nshtgn nepchs vtzxss = do
    (cstss,_,_) <- unzip3 <$> mapM (uncurry (shotgunFitMixtureLikelihood intl optm nshtgn nepchs)) vtzxss
    let csts' = catMaybes $ do
            csts <- cstss
            return $ descentMinimum <$> maybeLast csts
    return $ if null csts'
                then Nothing
                else Just $ average csts'
        where maybeLast [] = Nothing
              maybeLast xs = Just $ last xs

crossValidateEmpiricalCovariances
    :: KnownNat k
    => [([(Response k, Double)],[(Response k, Double)])] -- ^ K-folded data
    -> Double
{-# INLINE crossValidateEmpiricalCovariances #-}
crossValidateEmpiricalCovariances vtzxss = average $ do
    (vzxs,tzxs) <- vtzxss
    let vnrms = transition . fst <$> dataCovariances vzxs
        tnrms = transition . fst <$> dataCovariances tzxs
    return . average $ zipWith relativeEntropy tnrms vnrms



--- CSV ---


data PreferredStimuli = PreferredStimuli
    { psBinCentre :: Double
    , preferredStimulusCount :: Int
    , preferredStimulusDensity :: Double }
    deriving (Show,Generic)

instance ToNamedRecord PreferredStimuli where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered PreferredStimuli where
    headerOrder = goalCSVOrder

data Precisions = Precisions
    { pBinCentre :: Double
    , precisionCount :: Int
    , precisionDensity :: Double }
    deriving (Show,Generic)

instance ToNamedRecord Precisions where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered Precisions where
    headerOrder = goalCSVOrder

data ParameterDistributionFit = ParameterDistributionFit
    { parameterValue :: Double
    , densityFit :: Double
    , trueDensity :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord ParameterDistributionFit where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered ParameterDistributionFit where
    headerOrder = goalCSVOrder

data Gains = Gains
    { gBinCentre :: Double
    , gainCount :: Int
    , gainDensity :: Double }
    deriving (Show,Generic)

instance ToNamedRecord Gains where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered Gains where
    headerOrder = goalCSVOrder

data DataDependenceStats = DataDependenceStats
    { trainingSetSize :: Int
    , dataDependenceAverage :: Double
    , dataDependenceMinimum :: Double
    , dataDependenceMaximum :: Double }
    deriving (Show,Generic)

instance ToNamedRecord DataDependenceStats where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered DataDependenceStats where
    headerOrder = goalCSVOrder

data LikelihoodBounds = LikelihoodBounds
    { independentUpperBound :: Double
    , trueLowerBound :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord LikelihoodBounds where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered LikelihoodBounds where
    headerOrder = goalCSVOrder

data CrossValidationStats = CrossValidationStats
    { numberOfMixtureComponents :: Int
    , cvNegativeLogLikelihood :: Double
    , informationGain :: Double }
    deriving (Show,Generic)

--data CrossValidationStats = CrossValidationStats
--    { numberOfMixtureComponents :: Int
--    , modelNegativeLogLikelihood :: Double
--    , modelCovarianceDivergence :: Double
--    , empiricalCovarianceDivergence :: Double }
--    deriving (Show,Generic)

instance ToNamedRecord CrossValidationStats where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered CrossValidationStats where
    headerOrder = goalCSVOrder

data TuningCurves k = TuningCurves
    { tcStimulus :: Double
    , tuningCurves :: S.Vector k Double }
    deriving (Show,Generic)

instance KnownNat k => ToNamedRecord (TuningCurves k) where
    toNamedRecord (TuningCurves stm tcs) =
        let stmrc = "Stimulus" .= stm
            tcrcs = countRecords "Tuning Curve" tcs
         in namedRecord $ stmrc : tcrcs

instance KnownNat k => DefaultOrdered (TuningCurves k) where
    headerOrder _ =
        orderedHeader $ "Stimulus" : countHeaders "Tuning Curve" (Proxy @ k)

data GainProfiles m = GainProfiles
    { gpStimulus :: Double
    , gainProfiles :: S.Vector m Double }
    deriving (Show,Generic)

instance KnownNat m => ToNamedRecord (GainProfiles m) where
    toNamedRecord (GainProfiles stm gns) =
        let stmrc = "Stimulus" .= stm
            gnrcs = countRecords "Gain" gns
         in namedRecord $ stmrc : gnrcs

instance KnownNat m => DefaultOrdered (GainProfiles m) where
    headerOrder _ =
        orderedHeader $ "Stimulus" : countHeaders "Gain" (Proxy @ m)

data ConjugacyCurves m = ConjugacyCurves
    { rcStimulus :: Double
    , sumOfTuningCurves :: S.Vector m Double
    , conjugacyCurve :: Double }
    deriving (Show,Generic)

instance KnownNat m => ToNamedRecord (ConjugacyCurves m) where
    toNamedRecord (ConjugacyCurves stm stcs cnj) =
        let stmrc = "Stimulus" .= stm
            cnjrc = "Conjugacy Curve" .= cnj
            stcrc = countRecords "Sum of Tuning Curves" stcs
         in namedRecord $ stmrc : stcrc ++ [cnjrc]

instance KnownNat m => DefaultOrdered (ConjugacyCurves m) where
    headerOrder _ =
         orderedHeader $ "Stimulus" : countHeaders "Sum of Tuning Curves" (Proxy @ m) ++ ["Conjugacy Curve"]

data CategoryDependence m = CategoryDependence
    { cdStimulus :: Double
    , categoryStimulusDependence :: S.Vector m Double }
    deriving (Show,Generic)

instance KnownNat m => ToNamedRecord (CategoryDependence m) where
    toNamedRecord (CategoryDependence stm cts) =
        let stmrc = "Stimulus" .= stm
            ctrcs = countRecords "Category" cts
         in namedRecord $ stmrc : ctrcs

instance KnownNat m => DefaultOrdered (CategoryDependence m) where
    headerOrder _ = orderedHeader $ "Stimulus" : countHeaders "Category" (Proxy @ m)


--- Internal ---


countHeaders :: KnownNat n => String -> Proxy n -> [ByteString]
countHeaders ttl prxn = [ fromString (ttl ++ ' ' : show i) | i <- take (natValInt prxn) [(0 :: Int)..] ]

countRecords
    :: forall x n . (ToField x, KnownNat n, Storable x)
    => String
    -> S.Vector n x
    -> [(ByteString,ByteString)]
countRecords ttl = zipWith (.=) (countHeaders ttl $ Proxy @ n) . S.toList

-- | Returns (validation,training) pairs
kFoldDataset :: Int -> [x] -> [([x],[x])]
{-# INLINE kFoldDataset #-}
kFoldDataset k xs =
    let nvls = ceiling . (/(fromIntegral k :: Double)) . fromIntegral $ length xs
     in L.unfoldr unfoldFun ([], breakEvery nvls xs)
    where unfoldFun (_,[]) = Nothing
          unfoldFun (hds,tl:tls) = Just ((tl,concat $ hds ++ tls),(tl:hds,tls))

kFoldConditionalDataset
    :: KnownNat k
    => Int -> [(Response k, Double)] -> [([(Response k, Double)],[(Response k, Double)])]
{-# INLINE kFoldConditionalDataset #-}
kFoldConditionalDataset k zxs =
    let mp = M.fromListWith (++) . zip (snd <$> zxs) $ (:[]) <$> zxs
     in [ foldr1 folder vtzxss | vtzxss <- L.transpose . M.elems $ kFoldDataset k <$> mp ]
    where folder (vzxs',tzxs') (vzxs,tzxs) = (vzxs' ++ vzxs,tzxs' ++ tzxs)

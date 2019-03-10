{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    Arrows,
    TypeApplications,
    DeriveGeneric,
    OverloadedStrings,
    FlexibleContexts,
    TypeOperators
    #-}

module NeuralData.Mixture
    ( -- * Mixture
      getFittedMixtureLikelihood
    , strengthenMixtureLikelihood
    , randomMixtureLikelihood
    -- * Fitting
    , fitMixtureLikelihood
    , shotgunFitMixtureLikelihood
    , crossValidateMixtureLikelihood
    , mixtureLikelihoodDataDependence
    -- * Analyses
    , kFoldDataset
    , analyzePopulationCurves
    , preferredStimulusHistogram
    , precisionsHistogram
    , gainsHistograms
    , runPopulationParameterAnalyses
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Foreign.Storable
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import NeuralData
import NeuralData.VonMises hiding (ParameterCounts,ParameterDistributionFit)
import Data.String

import qualified Data.List as L

--- Types ---


--- Inference ---


getFittedMixtureLikelihood
    :: String
    -> String
    -> IO (NatNumber,NatNumber,[Double])
getFittedMixtureLikelihood expnm dst = do
    (k,n,xs) <- read . fromJust <$> goalReadDataset (Experiment prjnm expnm) (dst ++ "-parameters")
    return (k,n,xs)

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Mean #> Natural # MixtureGLM n (Neurons k) VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs

randomMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Source # LogNormal
    -> Source # VonMises
    -> Source # LogNormal
    -> Random r (Mean #> Natural # MixtureGLM n (Neurons k) VonMises)
{-# INLINE randomMixtureLikelihood #-}
randomMixtureLikelihood rmxs sgns sprf sprcs = do
    mxs <- samplePoint rmxs
    gns <- S.replicateM $ randomGains sgns
    tcs <- randomTuningCurves sprf sprcs
    let nctgl = toNatural . Point @ Source $ S.tail mxs
    return $ joinVonMisesMixturePopulationEncoder nctgl (S.map toNatural gns) tcs


--- Fitting ---


fitMixtureLikelihood
    :: forall r k n . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Natural # LogNormal -- ^ Initial Precision Distribution
    -> [(Response k,Double)]
    -> Random r [Mean #> Natural # MixtureGLM n (Neurons k) VonMises]
{-# INLINE fitMixtureLikelihood #-}
fitMixtureLikelihood eps nbtch nepchs rmxs rkp zxs = do
    kps <- S.replicateM $ samplePoint rkp
    mxs <- samplePoint rmxs
    let (zs,xs) = unzip zxs
        mus = S.generate $ \fnt ->
            let zis = fromIntegral . (`B.index` fnt) <$> zs
             in weightedCircularAverage $ zip zis xs
        sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps mus
    let gns0 = transition . sufficientStatisticT $ fst <$> zxs
        mtx = snd . splitAffine $ joinVonMisesPopulationEncoder (Right gns0) sps
    gnss' <- S.replicateM $ Point <$> S.replicateM (uniformR (-10,0))
    let gnss = S.map (gns0 <+>) gnss'
        nctgl = toNatural . Point @ Source $ S.tail mxs
        mxmdl = joinMixtureModel gnss nctgl
        mxlkl0 = joinBottomSubLinear mxmdl mtx
        grdcrc = loopCircuit' mxlkl0 $ proc (zxs',mxlkl) -> do
            let (zs',xs') = unzip zxs'
                dmxlkl = vanillaGradient
                    $ mixtureStochasticConditionalCrossEntropyDifferential xs' zs' mxlkl
            gradientCircuit eps defaultAdamPursuit -< joinTangentPair mxlkl dmxlkl
    streamCircuit grdcrc . take nepchs . breakEvery nbtch $ cycle zxs

data CrossEntropyDescentStats = CrossEntropyDescentStats
    { descentAverage :: Double
    , descentMinimum :: Double
    , descentMaximum :: Double }
    deriving (Show,Generic)

instance ToNamedRecord CrossEntropyDescentStats where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered CrossEntropyDescentStats where
    headerOrder = goalCSVOrder

costsToCrossEntropyDescentStats :: [Double] -> CrossEntropyDescentStats
{-# INLINE costsToCrossEntropyDescentStats #-}
costsToCrossEntropyDescentStats csts =
     CrossEntropyDescentStats (average csts) (minimum csts) (maximum csts)

data CrossValidationStats = CrossValidationStats
    { numberOfMixtureComponents :: Int
    , crossValidationAverage :: Double
    , crossValidationMinimum :: Double
    , crossValidationMaximum :: Double }
    deriving (Show,Generic)

instance ToNamedRecord CrossValidationStats where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered CrossValidationStats where
    headerOrder = goalCSVOrder

kFoldDataset :: Int -> [x] -> [([x],[x])]
kFoldDataset k xs =
    let nvls = ceiling . (/(fromIntegral k :: Double)) . fromIntegral $ length xs
     in L.unfoldr unfoldFun ([], breakEvery nvls xs)
    where unfoldFun (_,[]) = Nothing
          unfoldFun (hds,tl:tls) = Just ((tl,concat $ hds ++ tls),(tl:hds,tls))

shotgunFitMixtureLikelihood
    :: (KnownNat k, KnownNat m)
    => Int -- ^ Number of parallel fits
    -> Double -- ^ Learning rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> Natural # LogNormal -- ^ Initial Precisions
    -> [(Response k, Double)] -- ^ Training data
    -> Random r ([CrossEntropyDescentStats],Int, [Mean #> Natural # MixtureGLM m (Neurons k) VonMises])
{-# INLINE shotgunFitMixtureLikelihood #-}
shotgunFitMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm zxs = do
    let (zs,xs) = unzip zxs
    mlklss <- replicateM nshtgn $ fitMixtureLikelihood eps nbtch nepchs drch lgnrm zxs
    let cost = mixtureStochasticConditionalCrossEntropy xs zs
        mlklcstss = do
            mlkls <- mlklss
            return . zip mlkls $ cost <$> mlkls
        (nanmlklcstss,mlklcstss') = L.partition (any (\x -> isNaN x || isInfinite x) . map snd) mlklcstss
        mxmlkls = map fst . L.minimumBy (comparing (snd . last)) $ mlklcstss'
        sgdnrms = costsToCrossEntropyDescentStats <$> L.transpose (map (map snd) mlklcstss')
    return (sgdnrms,length nanmlklcstss, mxmlkls)

validateMixtureLikelihood
    :: (KnownNat k, KnownNat m)
    => Int -- ^ Number of parallel fits
    -> Double -- ^ Learning rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> Natural # LogNormal -- ^ Initial Pre,Validation data
    -> ([(Response k, Double)],[(Response k, Double)]) -- ^ Validation/Training data
    -> Random r Double
{-# INLINE validateMixtureLikelihood #-}
validateMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm (vzxs,tzxs) = do
    (_,_,mxmlkl) <- shotgunFitMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm tzxs
    let (vzs,vxs) = unzip vzxs
    return . mixtureStochasticConditionalCrossEntropy vxs vzs $ last mxmlkl

crossValidateMixtureLikelihood
    :: forall k m r . (KnownNat k, KnownNat m)
    => Int -- ^ Number of cross validation folds
    -> Int -- ^ Number of parallel fits
    -> Double -- ^ Learning rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> Natural # LogNormal -- ^ Initial Precisions
    -> [(Response k, Double)] -- ^ Data
    -> Random r CrossValidationStats
{-# INLINE crossValidateMixtureLikelihood #-}
crossValidateMixtureLikelihood kfld nshtgn eps nbtch nepchs drch lgnrm zxs = do
    let tvxzss = kFoldDataset kfld zxs
    csts <- mapM (validateMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm) tvxzss
    return $ CrossValidationStats (natValInt $ Proxy @ (m+1)) (average csts) (minimum csts) (maximum csts)

upperBoundsToCrossEntropyDescentStats :: Int -> [Double] -> DataDependenceStats
{-# INLINE upperBoundsToCrossEntropyDescentStats #-}
upperBoundsToCrossEntropyDescentStats sz csts =
     DataDependenceStats sz (average csts) (minimum csts) (maximum csts)


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

mixtureLikelihoodDataDependence
    :: forall k m r . (KnownNat k, KnownNat m)
    => Int -- ^ Number of data sets to generate
    -> Int -- ^ Number of shotgun fits
    -> Double -- ^ Learning rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> Natural # LogNormal -- ^ Initial Precisions
    -> Mean #> Natural # MixtureGLM m (Neurons k) VonMises -- ^ True Mixture Model
    -> Sample VonMises
    -> Int -- ^ Data set size
    -> Random r ([DataDependenceStats], Mean #> Natural # MixtureGLM m (Neurons k) VonMises)
{-# INLINE mixtureLikelihoodDataDependence #-}
mixtureLikelihoodDataDependence nsts nshtgn eps nbtch nepchs drch lgnrm plkl xsmps sz = do
    zxss <- replicateM nsts $ sampleMixtureLikelihood plkl sz
    (_,_,qlklss) <- unzip3 <$> mapM (shotgunFitMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm) zxss
    let qlkls = fst . L.minimumBy (comparing snd) $ do
            qlkls' <- qlklss
            return (qlkls', conditionalMixtureRelativeEntropyUpperBound xsmps plkl $ last qlkls')
    let ddsts = map (upperBoundsToCrossEntropyDescentStats sz) . L.transpose
            $ map (map $ conditionalMixtureRelativeEntropyUpperBound xsmps plkl) qlklss
    return (ddsts, last qlkls)
    -- = kFoldDataset kfld zxs
    --csts <- mapM (validateMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm) tvxzss
    --return $ CrossValidationStats (natValInt $ Proxy @ (m+1)) (average csts) (minimum csts) (maximum csts)

sampleMixtureLikelihood
    :: (KnownNat k, KnownNat m)
    => Mean #> Natural # MixtureGLM m (Neurons k) VonMises
    -> Int
    -> Random r [(Response k, Double)]
sampleMixtureLikelihood mlkl sz = do
    xs <- sample sz uni
    zs <- mapM (fmap hHead . samplePoint) $ mlkl >$>* xs
    return $ zip zs xs
        where uni :: Source # VonMises
              uni = Point $ S.doubleton 0 0

runPopulationParameterAnalyses
    :: (KnownNat k, KnownNat m)
    => Experiment
    -> String
    -> [Double]
    -> Int
    -> String
    -> String
    -> Maybe (Natural # LogNormal)
    -> Maybe [Natural # LogNormal]
    -> ((Mean #> Natural) # MixtureGLM m (Neurons k) VonMises)
    -> IO ()
runPopulationParameterAnalyses expmnt dst xsmps nbns rltv ttl mprcsdst mgndsts mlkl = do

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

    let crsbexp = Just $ Analysis ("noise-correlations/" ++ ttl) dst

    let (mtxln:mtxlns) = do
            let smlkl = sortVonMisesMixturePopulationEncoder mlkl
            mtx <- mixturePopulationNoiseCorrelations smlkl <$> xsmps
            return $ S.toList <$> S.toList (S.toRows mtx)

    goalExport True expmnt crsbexp mtxln
    mapM_ (goalExport False expmnt crsbexp) mtxlns

    let aniopts = defaultGnuplotOptions { whetherPNG = False, whetherAnimate = True }
        crgpi = rltv ++ "noise-correlations.gpi"
    runGnuplot expmnt crsbexp aniopts crgpi



--- Analysis ---


countHeaders :: KnownNat n => String -> Proxy n -> [ByteString]
countHeaders ttl prxn = [ fromString (ttl ++ ' ' : show i) | i <- take (natValInt prxn) [(0 :: Int)..] ]

countRecords
    :: forall x n . (ToField x, KnownNat n, Storable x)
    => String
    -> S.Vector n x
    -> [(ByteString,ByteString)]
countRecords ttl = zipWith (.=) (countHeaders ttl $ Proxy @ n) . S.toList


--- Population Curves ---


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

analyzePopulationCurves
    :: forall m k . (KnownNat m, KnownNat k)
    => Sample VonMises
    -> (Mean #> Natural # MixtureGLM m (Neurons k) VonMises)
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
            (rts,ct) <- splitMixtureModel <$> mgnxs
            let sct = coordinates $ toSource ct
            return (S.cons (1 - S.sum sct) sct, S.map potential rts)
        ccrvs = L.zipWith3 ConjugacyCurves smps stcs cnjs
        ctcrvs = L.zipWith CategoryDependence smps cts
     in (tcs,gps,ccrvs,ctcrvs)


--- Population Histrograms ---


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


preferredStimulusHistogram
    :: (KnownNat k, KnownNat m)
    => Int
    -> (Mean #> Natural # MixtureGLM m (Neurons k) VonMises)
    -> Maybe (Natural # VonMises)
    -> ([PreferredStimuli], [ParameterDistributionFit], Natural # VonMises)
{-# INLINE preferredStimulusHistogram #-}
preferredStimulusHistogram nbns mlkl mtrudns =
    let (_,_,nxs) = splitVonMisesMixturePopulationEncoder mlkl
        mus =  head . listCoordinates . toSource <$> S.toList nxs
        (bns,[cnts],[wghts]) = histograms nbns Nothing [mus]
        backprop = pairTangentFunction $ stochasticCrossEntropyDifferential mus
        vm0 = Point $ S.doubleton 0.01 0.01
        vm :: Natural # VonMises
        vm = vanillaGradientSequence backprop (-0.1) defaultAdamPursuit vm0 !! 500
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
    -> (Mean #> Natural # MixtureGLM m (Neurons k) VonMises)
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
    -> Mean #> Natural # MixtureGLM m (Neurons k) VonMises
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


--multiFitMixtureLikelihood
--    :: (KnownNat k, KnownNat m)
--    => Int -- ^ Number of parallel fits
--    -> Double -- ^ Learning rate
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of epochs
--    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
--    -> Natural # LogNormal -- ^ Initial Precisions
--    -> [(Response k, Double)] -- ^ Training data
--    -> [(Response k, Double)] -- ^ Validation data
--    -> Random r ([CrossEntropyDescentStats], Int, Mean #> Natural # MixtureGLM m (Neurons k) VonMises)
--{-# INLINE multiFitMixtureLikelihood #-}
--multiFitMixtureLikelihood npop eps nbtch nepchs drch lgnrm tzxs vzxs = do
--    let (vzs,vxs) = unzip vzxs
--    mlklss <- replicateM npop $ fitMixtureLikelihood eps nbtch nepchs drch lgnrm tzxs
--    let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
--        mlklcstss = do
--            mlkls <- mlklss
--            return . zip mlkls $ cost <$> mlkls
--        (nanmlklcstss,mlklcstss') = L.partition (any (\x -> isNaN x || isInfinite x) . map snd) mlklcstss
--        sgdnrms = costsToCrossEntropyDescentStats <$> L.transpose (map (map snd) mlklcstss')
--        mxmlkl = fst . L.minimumBy (comparing snd) $ last <$> mlklcstss'
--    return (sgdnrms, length nanmlklcstss, mxmlkl)

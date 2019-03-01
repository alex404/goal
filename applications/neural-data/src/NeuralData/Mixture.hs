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
      fitMixtureLikelihood
    , getFittedMixtureLikelihood
    , strengthenMixtureLikelihood
    , randomMixtureLikelihood
    -- * Analyses
    , analyzePopulationCurves
    , preferredStimulusHistogram
    , precisionsHistogram
    , gainsHistograms
    , multiFitMixtureLikelihood
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
    (k,n,xs) <- read . fromJust <$> goalReadDataset (Experiment prjnm expnm) dst
    return (k,n,xs)

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs

randomMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Source # LogNormal
    -> Source # VonMises
    -> Source # LogNormal
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
randomMixtureLikelihood rmxs sgns sprf sprcs = do
    mxs <- samplePoint rmxs
    gns <- S.replicateM $ randomGains sgns
    tcs <- randomTuningCurves sprf sprcs
    let nctgl = toNatural . Point @ Source $ S.tail mxs
    return $ joinVonMisesMixturePopulationEncoder nctgl (S.map toNatural gns) tcs

fitMixtureLikelihood
    :: forall r k n . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Natural # LogNormal -- ^ Initial Precision Distribution
    -> [(Response k,Double)]
    -> Random r [Mean #> Natural # MixtureGLM (Neurons k) n VonMises]
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
            mxlkl' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair mxlkl dmxlkl
            returnA -< mxlkl'
    streamCircuit grdcrc . take nepchs . breakEvery nbtch $ cycle zxs


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
    , conjugacyCurve :: Double
    , sumOfTuningCurves :: S.Vector m Double }
    deriving (Show,Generic)

instance KnownNat m => ToNamedRecord (ConjugacyCurves m) where
    toNamedRecord (ConjugacyCurves stm cnj stcs) =
        let stmrc = "Stimulus" .= stm
            cnjrc = "Conjugacy Curve" .= cnj
            stcrc = countRecords "Sum of Tuning Curves" stcs
         in namedRecord $ stmrc : cnjrc : stcrc

instance KnownNat m => DefaultOrdered (ConjugacyCurves m) where
    headerOrder _ =
         orderedHeader $ "Stimulus" : "Conjugacy Curve" : countHeaders "Sum of Tuning Curves" (Proxy @ m)

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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
    -> ([TuningCurves k],[GainProfiles (m+1)],[ConjugacyCurves (m+1)],[CategoryDependence (m+1)])
analyzePopulationCurves smps mlkl =
    let (_,ngnss,tcrws) = splitVonMisesMixturePopulationEncoder mlkl
        nrmlkl = joinVonMisesPopulationEncoder (Left 1) tcrws
        tcs = zipWith TuningCurves smps $ coordinates . toSource <$> nrmlkl >$>* smps
        mus = head . listCoordinates . toSource <$> S.toList tcrws
        gps = zipWith GainProfiles mus . S.toList . S.toRows . S.fromColumns
            $ S.map (coordinates . toSource) ngnss
        mgnxs = mlkl >$>* smps
        cnjs = potential <$> mgnxs
        (cts,stcs) = unzip $ do
            (rts,ct) <- splitMixtureModel <$> mgnxs
            let sct = coordinates $ toSource ct
            return (S.cons (1 - S.sum sct) sct, S.map potential rts)
        ccrvs = L.zipWith3 ConjugacyCurves smps cnjs stcs
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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
    -> Maybe (Natural # VonMises)
    -> ([PreferredStimuli], [ParameterDistributionFit], Natural # VonMises)
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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
    -> Maybe (Natural # LogNormal)
    -> ([Precisions], [ParameterDistributionFit], Natural # LogNormal)
precisionsHistogram nbns mlkl mtrudns =
    let (_,_,nxs) = splitVonMisesMixturePopulationEncoder mlkl
        rhos = head . tail . listCoordinates . toSource <$> S.toList nxs
        ((bns,cnts,wghts),(xs,dnss,mdnss),lgnrm) = logNormalHistogram nbns mtrudns rhos
     in (zipWith3 Precisions bns cnts wghts, zipWith3 ParameterDistributionFit xs dnss mdnss, lgnrm)

gainsHistograms
    :: (KnownNat k, KnownNat m)
    => Int
    -> Mean #> Natural # MixtureGLM (Neurons k) m VonMises
    -> Maybe [Natural # LogNormal]
    -> ([[Gains]], [[ParameterDistributionFit]], [Natural # LogNormal])
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


--- Mixture Cross Validation ---


data SGDDistribution = SGDDistribution
    { descentMean :: Double
    , descentStandardDeviation :: Double }
    deriving (Show,Generic)

instance ToNamedRecord SGDDistribution where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered SGDDistribution where
    headerOrder = goalCSVOrder

normalToSGDDistribution :: Source # Normal -> SGDDistribution
normalToSGDDistribution nrm =
    let (mu,var) = S.toPair $ coordinates nrm
     in SGDDistribution mu $ sqrt var

multiFitMixtureLikelihood
    :: (KnownNat k, KnownNat m)
    => Int -- ^ Number of parallel fits
    -> Double -- ^ Learning rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> Natural # LogNormal -- ^ Initial Precisions
    -> [(Response k, Double)] -- ^ Training data
    -> [(Response k, Double)] -- ^ Validation data
    -> Prob (ST r) ([SGDDistribution],Int, Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
multiFitMixtureLikelihood npop eps nbtch nepchs drch lgnrm tzxs vzxs = do
    let (vzs,vxs) = unzip vzxs
    mlklss <- replicateM npop $ fitMixtureLikelihood eps nbtch nepchs drch lgnrm tzxs
    let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
        mlklcstss = do
            mlkls <- mlklss
            return . zip mlkls $ cost <$> mlkls
        (nanmlklcstss,mlklcstss') = L.partition (any isNaN . map snd) mlklcstss
        sgdnrms = normalToSGDDistribution . mle <$> L.transpose (map (map snd) mlklcstss')
        mxmlkl = fst . L.minimumBy (comparing snd) $ last <$> mlklcstss'
    return (sgdnrms, length nanmlklcstss, mxmlkl)

--        tracer mlkls =
--            let nanbl = any isNaN . listCoordinates $ last mlkls
--                weighter = listCoordinates . toSource . snd
--                    . splitMixtureModel . fst . splitBottomSubLinear
--                iwghts = weighter $ head mlkls
--                lwghts = weighter $ last mlkls
--                trcstr = concat [ "\nAny NaNs?\n"
--                                , show nanbl
--                                , "\nInitial Mixture Weights:\n"
--                                , show iwghts
--                                , "\nFinal Mixture Weights:\n"
--                                , show lwghts ]
--             in trace trcstr mlkls


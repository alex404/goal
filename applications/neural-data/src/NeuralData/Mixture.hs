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
    (
    -- * Initialization
      Initialization
          ( RandomInitialization
          , DataInitialization
          , PartitionedDataInitialization )
    , initializeMixtureLikelihood
    -- * Fitting
    , Optimization
        ( StochasticGradientDescent
        , ExpectationMaximization
        , Alternating
        , Hybrid
        , Hybrid2 )
    , fitMixtureLikelihood
    , shotgunFitMixtureLikelihood
    , LikelihoodBounds (LikelihoodBounds)
    -- ** Stimulus Independent
    , shotgunFitMixture
    , shotgunEMFitMixture
    -- * Analyses
    , analyzePopulationCurves
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
    , multipartitionList
    , getMixtureLikelihood
    , strengthenMixtureLikelihood
    , dataCovarianceFunction
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Foreign.Storable
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import NeuralData
import NeuralData.VonMises hiding (ParameterCounts,ParameterDistributionFit)
import Data.String

import qualified Data.List as L
import qualified Data.Map as M
import Data.Tuple

--- Types ---


--- Inference ---


getMixtureLikelihood
    :: String -- ^ Experiment name
    -> String -- ^ Dataset
    -> IO (Maybe (NatNumber,NatNumber,[Double]))
getMixtureLikelihood expnm dst =
    fmap read <$> goalReadDataset (Experiment prjnm expnm) (dst ++ "-parameters")

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs

dataCovarianceFunction :: KnownNat k => [(Response k,Double)] -> Double -> Source # MultivariateNormal k
{-# INLINE dataCovarianceFunction #-}
dataCovarianceFunction zxs x =
    let mpnrms = M.fromList $ swap <$> dataCovariances zxs
        kys = M.keys mpnrms
        ky = last $ last kys : takeWhile (<= x) kys
     in mpnrms M.! ky


--- Initialization ---


data Initialization k n =
    RandomInitialization
        (Natural # Dirichlet (n+1))
        (Source # LogNormal)
        (Source # VonMises)
        (Source # LogNormal)
    | DataInitialization
        (Natural # Dirichlet (n+1)) Double Int Int (Double,Double) [(Response k,Double)]
    | PartitionedDataInitialization
        (Natural # Dirichlet (n+1)) Double Int Int [(Response k,Double)]

initializeCategorical
    :: KnownNat n
    => Natural # Dirichlet (n+1)
    -> Random r (Natural # Categorical Int n)
initializeCategorical drch =
    toNatural . Point @ Source . S.tail <$> samplePoint drch

initializeMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Initialization k n
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
initializeMixtureLikelihood (RandomInitialization drch sgns sprf sprcs) = do
    gns <- S.replicateM $ randomGains sgns
    tcs <- randomTuningCurves sprf sprcs
    ncts <- initializeCategorical drch
    return $ joinVonMisesMixturePopulationEncoder ncts (S.map toNatural gns) tcs
initializeMixtureLikelihood (DataInitialization drch eps nepchs nbtch mnmx zxs) = do
    gnss <- gainsFromData mnmx zxs
    let tcs = tuningCurvesFromData eps nepchs nbtch zxs
    ncts <- initializeCategorical drch
    return $ joinVonMisesMixturePopulationEncoder ncts gnss tcs
initializeMixtureLikelihood (PartitionedDataInitialization drch eps nepchs nbtch zxs) = do
    ncts <- initializeCategorical drch
    let gnss = gainsFromPartitionedData ncts zxs
    let tcs = tuningCurvesFromData eps nepchs nbtch zxs
    return $ joinVonMisesMixturePopulationEncoder ncts gnss tcs

preferredStimuliFromData
    :: forall k . KnownNat k
    => [(Response k,Double)]
    -> S.Vector k Double
{-# INLINE preferredStimuliFromData #-}
preferredStimuliFromData zxs =
    let (zs,xs) = unzip zxs
     in S.generate $ \fnt ->
         let zis = fromIntegral . (`B.index` fnt) <$> zs
          in weightedCircularAverage $ zip zis xs

tuningCurvesFromData
    :: forall k . KnownNat k
    => Double
    -> Int
    -> Int
    -> [(Response k,Double)]
    -> S.Vector k (Natural # VonMises)
tuningCurvesFromData eps nepchs nbtch zxs =
    let mus = preferredStimuliFromData zxs
        tcs0 = S.map (\mu -> toNatural . Point @ Source $ S.fromTuple (mu,1)) mus
        lkl0 = joinVonMisesPopulationEncoder (Left 1) tcs0
        grdcrc = loopCircuit' lkl0 $ proc (zxs',lkl) -> do
            let (zs',xs') = unzip zxs'
                dlkl = vanillaGradient $ stochasticConditionalCrossEntropyDifferential xs' zs' lkl
            lkl' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair lkl dlkl
            returnA -< lkl'
        lkl1 = runIdentity . iterateCircuit grdcrc . take nepchs . breakEvery nbtch $ cycle zxs
     in snd $ splitVonMisesPopulationEncoder lkl1

gainsFromData
    :: forall r k n . (KnownNat k, KnownNat n)
    => (Double,Double) -- ^ Uniform gain noise
    -> [(Response k,Double)]
    -> Random r (S.Vector (n+1) (Natural # Neurons k))
{-# INLINE gainsFromData #-}
gainsFromData mnmx zxs = do
    let gns0 = transition $ sufficientStatisticT $ fst <$> zxs
    gnss' <- S.replicateM $ Point <$> S.replicateM (uniformR mnmx)
    return $ S.map (gns0 <+>) gnss'

gainsFromPartitionedData
    :: forall k n . (KnownNat k, KnownNat n)
    => Natural # Categorical Int n -- ^ Categories
    -> [(Response k,Double)]
    -> S.Vector (n+1) (Natural # Neurons k)
{-# INLINE gainsFromPartitionedData #-}
gainsFromPartitionedData ctgl zxs =
    let zs = fst <$> L.sortOn snd zxs
        zss = multipartitionList ctgl zs
     in G.convert $ toNatural . sufficientStatisticT <$> zss

-- | Given a categorical distribution, partitions a list into sublists with
-- sizes given by the weights of the distribution.
multipartitionList
    :: (KnownNat n, Transition c Source (Categorical Int n))
    => c # Categorical Int n
    -> [x]
    -> B.Vector (n+1) [x]
multipartitionList cts0 xs0 = B.unfoldrN unfolder (xs0,0)
    where
        ctgl = toSource cts0
        ttl = fromIntegral $ length xs0
        unfolder (xs,k) =
              let spltn = round $ ttl * density ctgl k
                  (hdxs,tlxs) = splitAt spltn xs
               in (hdxs,(tlxs,k+1))


--- Fitting ---


data Optimization =
    StochasticGradientDescent Double Double Int
    | ExpectationMaximization Double Int Int
    | Alternating Int Double Double Int
    | Hybrid Double Double Int
    | Hybrid2 Double Int Int

fitMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Optimization -- ^ Training Data
    -> [(Response k,Double)] -- ^ Training Data
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises -- ^ Initial Likelihood
    -> Chain (Random r) (Mean #> Natural # MixtureGLM (Neurons k) n VonMises) -- ^ CE Descent
fitMixtureLikelihood (StochasticGradientDescent eps dcy nbtch) zxs mlkl0 =
    chain (sgdMixtureLikelihoodEpoch eps dcy nbtch zxs) mlkl0
fitMixtureLikelihood (ExpectationMaximization eps nbtch nstps) zxs mlkl0 =
    let (zs,xs) = unzip zxs
        gp = defaultAdamPursuit
     in chain (mixtureStochasticConditionalEMAscent eps gp nbtch nstps xs zs) mlkl0
fitMixtureLikelihood (Alternating nrpts eps dcy nbtch) zxs mlkl0 =
    chain (alternatingMixtureLikelihoodEpoch nrpts eps dcy nbtch zxs) ((True,0),mlkl0) >>^ snd
fitMixtureLikelihood (Hybrid eps dcy nbtch) zxs mlkl0 =
    chain (hybridMixtureLikelihoodEpoch eps dcy nbtch zxs) (False,mlkl0) >>^ snd
fitMixtureLikelihood (Hybrid2 eps nbtch nstps) zxs mlkl0 =
    chain (hybrid2MixtureLikelihoodEpoch eps nbtch nstps zxs) (False,mlkl0) >>^ snd


shotgunFitMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Initialization k n
    -> Optimization
    -> Int -- ^ Number of parallel fits
    -> Int -- ^ Number of epochs
    -> [(Response k, Double)] -- ^ Training data
    -> [(Response k, Double)] -- ^ Validation Data
    -> Random r ([CrossEntropyDescentStats],Int,Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
{-# INLINE shotgunFitMixtureLikelihood #-}
shotgunFitMixtureLikelihood intl optm nshtgn nepchs tzxs vzxs = do
    mlkl0s <- replicateM nshtgn $ initializeMixtureLikelihood intl
    mlklss <- sequence $ streamChain nepchs . fitMixtureLikelihood optm tzxs <$> mlkl0s
    let (vzs,vxs) = unzip vzxs
    let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
    return $ filterShotgun mlklss cost

hybridMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Double -- ^ Weight decay rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> (Bool, Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
    -> Random r (Bool, Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
{-# INLINE hybridMixtureLikelihoodEpoch #-}
hybridMixtureLikelihoodEpoch _ _ _ zxs (True,mlkl) = do
    let mlkl' = mixturePopulationPartialExpectationMaximization zxs mlkl
    return (False,mlkl')
hybridMixtureLikelihoodEpoch eps dcy nbtch zxs (False,mlkl) = do
    mlkl' <- sgdMixtureLikelihoodEpoch eps dcy nbtch zxs mlkl
    return (True,mlkl')

hybrid2MixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of steps
    -> [(Response k,Double)] -- ^ Training Data
    -> (Bool, Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
    -> Random r (Bool, Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
{-# INLINE hybrid2MixtureLikelihoodEpoch #-}
hybrid2MixtureLikelihoodEpoch _ _ _ zxs (True,mlkl) = do
    let mlkl' = mixturePopulationPartialExpectationMaximization zxs mlkl
    return (False,mlkl')
hybrid2MixtureLikelihoodEpoch eps nbtch nstps zxs (False,mlkl) = do
    let (zs,xs) = unzip zxs
        gp = defaultAdamPursuit
    mlkl' <- mixtureStochasticConditionalEMAscent eps gp nbtch nstps xs zs mlkl
    return (True,mlkl')

alternatingMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Int -- & Number of repeats of each algorithm
    -> Double -- ^ Learning Rate
    -> Double -- ^ Weight decay rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> ((Bool, Int), Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
    -> Random r ((Bool, Int), Mean #> Natural # MixtureGLM (Neurons k) n VonMises)
{-# INLINE alternatingMixtureLikelihoodEpoch #-}
alternatingMixtureLikelihoodEpoch nrpts _ _ _ zxs ((True,stp),mlkl) = do
    let mlkl' = mixturePopulationPartialExpectationMaximization zxs mlkl
        (bl',stp') = if stp == nrpts
                        then (False,0)
                        else (True,stp+1)
    return ((bl',stp'),mlkl')
alternatingMixtureLikelihoodEpoch nrpts eps dcy nbtch zxs ((False,stp),mlkl) = do
    mlkl' <- sgdSubMixtureLikelihoodEpoch eps dcy nbtch zxs mlkl
    let (bl',stp') = if stp == nrpts
                        then (True,0)
                        else (False,stp+1)
    return ((bl',stp'),mlkl')

-- | Note that the way it is at the moment might be problematic because it
-- constantly restarts the Adam variables.
sgdMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Double -- ^ Weight decay rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises -- ^ Initial Likelihood
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) n VonMises) -- ^ CE Descent
{-# INLINE sgdMixtureLikelihoodEpoch #-}
sgdMixtureLikelihoodEpoch eps dcy nbtch zxs0 mlkl0 = do
    zxs <- shuffleList zxs0
    let grdcrc = loopCircuit' mlkl0 $ proc (zxs',mlkl) -> do
            let (zs',xs') = unzip zxs'
                dmlkl = vanillaGradient
                    $ mixtureStochasticConditionalCrossEntropyDifferential xs' zs' mlkl
                dcymlkl = weightDecay dcy mlkl
            gradientCircuit eps defaultAdamPursuit -< joinTangentPair dcymlkl dmlkl
    return . runIdentity . iterateCircuit grdcrc . breakEvery nbtch $ zxs

sgdSubMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Double -- ^ Weight decay rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises -- ^ Initial Likelihood
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) n VonMises) -- ^ CE Descent
{-# INLINE sgdSubMixtureLikelihoodEpoch #-}
sgdSubMixtureLikelihoodEpoch eps dcy nbtch zxs0 mlkl0 = do
    zxs <- shuffleList zxs0
    let (hrm0,tcs0) = splitBottomSubLinear mlkl0
        (gnss,cat0) = splitBottomHarmonium hrm0
        grdcrc = loopCircuit' (tcs0,fromOneHarmonium cat0) $ proc (zxs',(tcs,cat)) -> do
            let mlkl = joinBottomSubLinear (joinBottomHarmonium gnss $ toOneHarmonium cat) tcs
                (zs',xs') = unzip zxs'
                dmlkls = dualIsomorphism
                    <$> zipWith stochasticMixtureDifferential ((:[]) <$> zs') (mlkl >$>* xs')
                (dgnss,dcats) = unzip $ splitBottomHarmonium <$> dmlkls
                dzs = fst . splitAffine <$> dgnss
                dtcs = primalIsomorphism . fst $ propagate dzs (sufficientStatistic <$> xs') tcs
                dcat = primalIsomorphism . fromOneHarmonium $ averagePoint dcats
                dcytcs = weightDecay dcy tcs
                dcycat = weightDecay dcy cat
            tcs' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair dcytcs $ vanillaGradient dtcs
            cat' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair dcycat $ vanillaGradient dcat
            returnA -< (tcs',cat')
    let (tcs1,cat1) = runIdentity . iterateCircuit grdcrc . breakEvery nbtch $ zxs
    return $ joinBottomSubLinear (joinBottomHarmonium gnss $ toOneHarmonium cat1) tcs1


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
    -> Mean #> Natural # MixtureGLM (Neurons k) m VonMises
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

                let trulkl :: Mean #> Natural # MixtureGLM (Neurons k) m' VonMises
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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
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
    -> (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
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
    -> Mean #> Natural # MixtureGLM (Neurons k) m VonMises
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
--    let xs = L.sort . L.nub $ snd <$> snd (head vtzxss)
--    let nnrmcv = average $ do
--            (mlkl,vzxs) <- zip mlkls $ fst <$> vtzxss
--            let vnrms = transition . fst <$> dataCovariances vzxs
--                tnrms = transition . mixturePopulationCovariance <$> mlkl >$>* xs
--            return . average $ zipWith relativeEntropy tnrms vnrms
    let csts' = catMaybes $ do
            csts <- cstss
            return $ descentMinimum <$> maybeLast csts
        cst = if null csts'
                 then Nothing
                 else Just $ average csts'
    return cst
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






--validateMixtureLikelihood
--    :: (KnownNat k, KnownNat m)
--    => Int -- ^ Number of parallel fits
--    -> Double -- ^ Learning rate
--    -> Double -- ^ Weight decay
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of epochs
--    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
--    -> Natural # LogNormal -- ^ Initial Pre,Validation data
--    -> ([(Response k, Double)],[(Response k, Double)]) -- ^ Validation/Training data
--    -> Random r Double
--{-# INLINE validateMixtureLikelihood #-}
--validateMixtureLikelihood shotgunFitMixtureLikelihood intl optm nshtgn nepchs zxs zxs = do
--    mlkl <- shotgunFitMixtureLikelihood intl optm nshtgn nepchs zxs zxs
--    mlkl0s <- replicateM nshtgn $ initializeMixtureLikelihood intl
--    mlklss <- sequence $ streamChain nepchs . fitMixtureLikelihood optm tzxs <$> mlkl0s
--    let (vzs,vxs) = unzip vzxs
--    let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
--    return $ filterShotgun mlklss cost
--
--
--    (_,_,mxmlkl) <- shotgunFitMixtureLikelihood intl optm nshtgn nepchs tzxs vzxs
--    let (vzs,vxs) = unzip vzxs
--    return . mixtureStochasticConditionalCrossEntropy vxs vzs $ last mxmlkl



--- Stimulus Independent ---


initializeMixture
    :: forall r k n . (KnownNat k, KnownNat n)
    => Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> [Response k]
    -> Random r (Natural # Mixture (Neurons k) n)
{-# INLINE initializeMixture #-}
initializeMixture rmxs zs = do
    mxs <- samplePoint rmxs
    let gns0 = transition $ sufficientStatisticT zs
    gnss' <- S.replicateM $ Point <$> S.replicateM (uniformR (-0.01,0.01))
    let gnss = S.map (gns0 <+>) gnss'
        nctgl = toNatural . Point @ Source $ S.tail mxs
    return $ joinMixture gnss nctgl

fitMixture
    :: forall r k n . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Double -- ^ Weight decay rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> [Response k]
    -> Random r [Natural # Mixture (Neurons k) n]
{-# INLINE fitMixture #-}
fitMixture eps dcy nbtch nepchs rmxs zs = do
    mmdl0 <- initializeMixture rmxs zs
    let grdcrc = loopCircuit' mmdl0 $ proc (zs',mmdl) -> do
            let dmmdl = vanillaGradient $ stochasticMixtureDifferential zs' mmdl
                dcymmdl = weightDecay dcy mmdl
            gradientCircuit eps defaultAdamPursuit -< joinTangentPair dcymmdl dmmdl
    streamCircuit grdcrc . take nepchs . breakEvery nbtch $ cycle zs

emFitMixture
    :: forall r k n . (KnownNat k, KnownNat n)
    => Int -- ^ Number of epochs
    -> Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> [Response k]
    -> Random r [Natural # Mixture (Neurons k) n]
{-# INLINE emFitMixture #-}
emFitMixture nepchs rmxs zs = do
    mmdl0 <- initializeMixture rmxs zs
    return . take nepchs $ iterate (mixtureExpectationMaximization zs) mmdl0

shotgunFitMixture
    :: (KnownNat k, KnownNat m)
    => Int -- ^ Number of parallel fits
    -> Double -- ^ Learning rate
    -> Double -- ^ Weight Decay
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> [Response k] -- ^ Training data
    -> Random r ([CrossEntropyDescentStats],Int,Natural # Mixture (Neurons k) m)
{-# INLINE shotgunFitMixture #-}
shotgunFitMixture nshtgn eps dcy nbtch nepchs drch zs = do
    mmdlss <- replicateM nshtgn $ fitMixture eps dcy nbtch nepchs drch zs
    let cost hrm = negate . average $ logMixtureDensity hrm <$> zs
    return $ filterShotgun mmdlss cost

shotgunEMFitMixture
    :: (KnownNat k, KnownNat m)
    => Int -- ^ Number of parallel fits
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
    -> [Response k] -- ^ Training data
    -> Random r ([CrossEntropyDescentStats],Int,Natural # Mixture (Neurons k) m)
{-# INLINE shotgunEMFitMixture #-}
shotgunEMFitMixture nshtgn nepchs drch zs = do
    mmdlss <- replicateM nshtgn $ emFitMixture nepchs drch zs
    let cost hrm = negate . average $ logMixtureDensity hrm <$> zs
    return $ filterShotgun mmdlss cost


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


filterShotgun :: Show a => [[a]] -> (a -> Double) -> ([CrossEntropyDescentStats], Int, a)
{-# INLINE filterShotgun #-}
filterShotgun mdlss cost =
    let mdlcstss = do
            mdls <- mdlss
            return . zip mdls $ cost <$> mdls
        (nanmdlcstss,mdlcstss') = L.partition (any (\x -> isNaN x || isInfinite x) . map snd) mdlcstss
        mxmdls = map fst . L.minimumBy (comparing (snd . last)) $ mdlcstss'
        sgdnrms = costsToCrossEntropyDescentStats <$> L.transpose (map (map snd) mdlcstss')
     in (sgdnrms, length nanmdlcstss, last mxmdls)

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

--- Graveyard ---


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
--    -> Random r ([CrossEntropyDescentStats], Int, Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
--{-# INLINE multiFitMixtureLikelihood #-}
--multiFitMixtureLikelihood npop eps nbtch nepchs drch lgnrm tzxs vzxs = do
--    let (vzs,vxs) = unzip vzxs
--    mlklss <- replicateM npop $ sgdFitMixtureLikelihood eps nbtch nepchs drch lgnrm tzxs
--    let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
--        mlklcstss = do
--            mlkls <- mlklss
--            return . zip mlkls $ cost <$> mlkls
--        (nanmlklcstss,mlklcstss') = L.partition (any (\x -> isNaN x || isInfinite x) . map snd) mlklcstss
--        sgdnrms = costsToCrossEntropyDescentStats <$> L.transpose (map (map snd) mlklcstss')
--        mxmlkl = fst . L.minimumBy (comparing snd) $ last <$> mlklcstss'
--    return (sgdnrms, length nanmlklcstss, mxmlkl)

--mixtureLikelihoodDataDependence
--    :: forall k m r . (KnownNat k, KnownNat m)
--    => Int -- ^ Number of data sets to generate
--    -> Int -- ^ Number of shotgun fits
--    -> Double -- ^ Learning rate
--    -> Double -- ^ Weight decay
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of epochs
--    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
--    -> Natural # LogNormal -- ^ Initial Precisions
--    -> Mean #> Natural # MixtureGLM (Neurons k) m VonMises -- ^ True Mixture Model
--    -> Sample VonMises
--    -> Int -- ^ Data set size
--    -> Random r ([DataDependenceStats], Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
--{-# INLINE mixtureLikelihoodDataDependence #-}
--mixtureLikelihoodDataDependence nsts nshtgn eps dcy nbtch nepchs drch lgnrm plkl xsmps sz = do
--    zxss <- replicateM nsts $ sampleMixtureLikelihood plkl sz
--    (_,_,qlklss) <- unzip3
--        <$> mapM (shotgunFitMixtureLikelihood nshtgn eps dcy nbtch nepchs drch lgnrm) zxss
--    let qlkls = fst . L.minimumBy (comparing snd) $ do
--            qlkls' <- qlklss
--            return (qlkls', conditionalMixtureRelativeEntropyUpperBound xsmps plkl $ last qlkls')
--    let ddsts = map (upperBoundsToCrossEntropyDescentStats sz) . L.transpose
--            $ map (map $ conditionalMixtureRelativeEntropyUpperBound xsmps plkl) qlklss
--    return (ddsts, last qlkls)
--    -- = kFoldDataset kfld zxs
--    --csts <- mapM (validateMixtureLikelihood nshtgn eps nbtch nepchs drch lgnrm) tvxzss
--    --return $ CrossValidationStats (natValInt $ Proxy @ (m+1)) (average csts) (minimum csts) (maximum csts)
--
--sampleMixtureLikelihood
--    :: (KnownNat k, KnownNat m)
--    => Mean #> Natural # MixtureGLM (Neurons k) m VonMises
--    -> Int
--    -> Random r [(Response k, Double)]
--sampleMixtureLikelihood mlkl sz = do
--    xs <- sample sz uni
--    zs <- mapM (fmap hHead . samplePoint) $ mlkl >$>* xs
--    return $ zip zs xs
--        where uni :: Source # VonMises
--              uni = Point $ S.doubleton 0 0

-- | Given a list of percentages, breaks a list into sublists of approximate
-- chunks with the given percenteges. If the sum of the percentages is smaller
-- than one, elements will be dropped, and if the sum is greater than one,
-- additional percentages will be ignored.
--multipartitionList
--    :: [Double]
--    -> [x]
--    -> [[x]]
--multipartitionList wghts xs0 =
--    let nsmps = length xs0
--        stps0 = round . (fromIntegral nsmps *) <$> wghts
--     in L.unfoldr unfolder (xs0,stps0)
--    where unfolder (x:xs,stp:stps) =
--              let (hdxs,tlxs) = splitAt stp (x:xs)
--               in Just (hdxs,(tlxs,stps))
--          unfolder _ = Nothing

--emFitMixtureLikelihood
--    :: forall k n . (KnownNat k, KnownNat n)
--    => [(Response k,Double)]
--    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises
--    -> [Mean #> Natural # MixtureGLM (Neurons k) n VonMises]
--{-# INLINE emFitMixtureLikelihood #-}
--emFitMixtureLikelihood zxs =
--    iterate (mixturePopulationExpectationMaximization zxs)
--
--shotgunFitMixtureLikelihood
--    :: (KnownNat k, KnownNat m)
--    => Int -- ^ Number of steps between performance evaluation
--    -> Int -- ^ Number of parallel fits
--    -> Double -- ^ Learning rate
--    -> Double -- ^ Weight Decay
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of epochs
--    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
--    -> Natural # LogNormal -- ^ Initial Precisions
--    -> [(Response k, Double)] -- ^ Training data
--    -> Random r ([CrossEntropyDescentStats],Int,[Mean #> Natural # MixtureGLM (Neurons k) m VonMises])
--{-# INLINE shotgunFitMixtureLikelihood #-}
--shotgunFitMixtureLikelihood skp nshtgn eps dcy nbtch nepchs drch lgnrm zxs = do
--    let (zs,xs) = unzip zxs
--    mlklss <- replicateM nshtgn $ sgdFitMixtureLikelihood eps dcy nbtch nepchs drch lgnrm zxs
--    let cost = mixtureStochasticConditionalCrossEntropy (take 1000 xs) (take 1000 zs)
--    return $ filterShotgun skp mlklss cost
--
--validateMixtureLikelihood
--    :: (KnownNat k, KnownNat m)
--    => Int -- ^ Number of parallel fits
--    -> Double -- ^ Learning rate
--    -> Double -- ^ Weight decay
--    -> Int -- ^ Batch size
--    -> Int -- ^ Number of epochs
--    -> Natural # Dirichlet (m + 1) -- ^ Initial mixture parameters
--    -> Natural # LogNormal -- ^ Initial Pre,Validation data
--    -> ([(Response k, Double)],[(Response k, Double)]) -- ^ Validation/Training data
--    -> Random r Double
--{-# INLINE validateMixtureLikelihood #-}
--validateMixtureLikelihood nshtgn eps dcy nbtch nepchs drch lgnrm (vzxs,tzxs) = do
--    (_,_,mxmlkl) <- shotgunFitMixtureLikelihood 100 nshtgn eps dcy nbtch nepchs drch lgnrm tzxs
--    let (vzs,vxs) = unzip vzxs
--    return . mixtureStochasticConditionalCrossEntropy vxs vzs $ last mxmlkl
--


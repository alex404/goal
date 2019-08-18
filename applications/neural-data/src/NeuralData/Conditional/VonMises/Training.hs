{-# LANGUAGE Arrows,TypeApplications,DeriveGeneric #-}

module NeuralData.Conditional.VonMises.Training
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
        , Hybrid
        , Hybrid2 )
    , fitIPLikelihood
    , fitMixtureLikelihood
    , shotgunFitMixtureLikelihood
    -- ** Cross Entropy calculations
    , CrossEntropyDescentStats (descentMinimum, descentMaximum, descentAverage)
    , costsToCrossEntropyDescentStats
    ) where


--- Imports ---


-- Goal --

import NeuralData.Conditional
import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Data.List as L


--- Initialization ---


-- Synthetic --

randomGains
    :: KnownNat k => Source # LogNormal -> Random r (Source # Neurons k)
randomGains = initialize

randomPreferredStimuli
    :: KnownNat k => Source # VonMises -> Random r (S.Vector k Double)
randomPreferredStimuli = S.replicateM . samplePoint

randomPrecisions
    :: KnownNat k => Source # LogNormal -> Random r (S.Vector k Double)
randomPrecisions = S.replicateM . samplePoint

randomTuningCurves
    :: KnownNat k
    => Source # VonMises
    -> Source # LogNormal
    -> Random r (S.Vector k (Natural # VonMises))
randomTuningCurves sprf sprcs = do
    mus <- randomPreferredStimuli sprf
    kps <- randomPrecisions sprcs
    let mukps = S.zipWith S.doubleton mus kps
    return $ S.map (toNatural . Point @ Source) mukps

-- Algorithms --

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
    -> Random r (Natural # Categorical n)
initializeCategorical drch =
    toNatural . Point @ Source . S.tail <$> samplePoint drch

initializeMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Initialization k n
    -> Random r (Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
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
         let zis = fromIntegral . (`S.index` fnt) <$> zs
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
            let dlkl = vanillaGradient $ conditionalLogLikelihoodDifferential zxs' lkl
            gradientCircuit eps defaultAdamPursuit -< (lkl,dlkl)
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
    => Natural # Categorical n -- ^ Categories
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
    :: (KnownNat n, Transition c Source (Categorical n))
    => c # Categorical n
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

--- Training ---


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


filterShotgun :: Show a => [[a]] -> (a -> Double) -> ([CrossEntropyDescentStats], Int, a)
{-# INLINE filterShotgun #-}
filterShotgun mdlss cost =
    let mdlcstss = do
            mdls <- mdlss
            return . zip mdls $ cost <$> mdls
        (nanmdlcstss,mdlcstss') = L.partition (any (\x -> isNaN x || isInfinite x) . map snd) mdlcstss
        mxmdls = map fst . L.maximumBy (comparing (snd . last)) $ mdlcstss'
        sgdnrms = costsToCrossEntropyDescentStats <$> L.transpose (map (map snd) mdlcstss')
     in (sgdnrms, length nanmdlcstss, last mxmdls)


data Optimization =
    StochasticGradientDescent Double Int
    | ExpectationMaximization Double Int Int
    | Hybrid Double Int
    | Hybrid2 Double Int Int

fitIPLikelihood
    :: KnownNat k
    => Double -- ^ Learning rate
    -> Int -- ^ Number of Iterations
    -> [(Response k,Double)] -- ^ Training Data
    -> Natural #> Affine Tensor (Neurons k) VonMises -- ^ CE Descent
fitIPLikelihood eps nepchs zxs =
    let gns = transition $ sufficientStatisticT $ fst <$> zxs
        mus = preferredStimuliFromData zxs
        tcs = S.map (\mu -> toNatural . Point @ Source $ S.fromTuple (mu,1)) mus
        lkl0 = joinVonMisesPopulationEncoder (Right gns) tcs
     in vanillaGradientSequence (conditionalLogLikelihoodDifferential zxs)
            eps defaultAdamPursuit lkl0 !! nepchs


fitMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Optimization -- ^ Training Data
    -> [(Response k,Double)] -- ^ Training Data
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises -- ^ Initial Likelihood
    -> Chain (Random r) (Natural #> ConditionalMixture (Neurons k) n Tensor VonMises) -- ^ CE Descent
fitMixtureLikelihood (StochasticGradientDescent eps nbtch) zxs mlkl0 =
    chain (sgdMixtureLikelihoodEpoch eps nbtch zxs) mlkl0
fitMixtureLikelihood (ExpectationMaximization eps nbtch nstps) zxs mlkl0 =
     chain (conditionalExpectationMaximizationAscent eps defaultAdamPursuit nbtch nstps zxs) mlkl0
fitMixtureLikelihood (Hybrid eps nbtch) zxs mlkl0 =
    chain (hybridMixtureLikelihoodEpoch eps nbtch zxs) (False,mlkl0) >>^ snd
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
    -> Random r ([CrossEntropyDescentStats],Int,Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
{-# INLINE shotgunFitMixtureLikelihood #-}
shotgunFitMixtureLikelihood intl optm nshtgn nepchs tzxs vzxs = do
    mlkl0s <- replicateM nshtgn $ initializeMixtureLikelihood intl
    mlklss <- sequence $ streamChain nepchs . fitMixtureLikelihood optm tzxs <$> mlkl0s
    let cost = conditionalLogLikelihood vzxs
    return $ filterShotgun mlklss cost

hybridMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> (Bool, Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
    -> Random r (Bool, Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
{-# INLINE hybridMixtureLikelihoodEpoch #-}
hybridMixtureLikelihoodEpoch _ _ zxs (True,mlkl) = do
    let mlkl' = mixturePopulationPartialExpectationMaximization zxs mlkl
    return (False,mlkl')
hybridMixtureLikelihoodEpoch eps nbtch zxs (False,mlkl) = do
    mlkl' <- sgdMixtureLikelihoodEpoch eps nbtch zxs mlkl
    return (True,mlkl')

hybrid2MixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of steps
    -> [(Response k,Double)] -- ^ Training Data
    -> (Bool, Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
    -> Random r (Bool, Natural #> ConditionalMixture (Neurons k) n Tensor VonMises)
{-# INLINE hybrid2MixtureLikelihoodEpoch #-}
hybrid2MixtureLikelihoodEpoch _ _ _ zxs (True,mlkl) = do
    let mlkl' = mixturePopulationPartialExpectationMaximization zxs mlkl
    return (False,mlkl')
hybrid2MixtureLikelihoodEpoch eps nbtch nstps zxs (False,mlkl) = do
    mlkl' <- conditionalExpectationMaximizationAscent eps defaultAdamPursuit nbtch nstps zxs mlkl
    return (True,mlkl')

-- | Note that the way it is at the moment might be problematic because it
-- constantly restarts the Adam variables.
sgdMixtureLikelihoodEpoch
    :: forall k n r . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> [(Response k,Double)] -- ^ Training Data
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises -- ^ Initial Likelihood
    -> Random r (Natural #> ConditionalMixture (Neurons k) n Tensor VonMises) -- ^ CE Descent
{-# INLINE sgdMixtureLikelihoodEpoch #-}
sgdMixtureLikelihoodEpoch eps nbtch zxs0 mlkl0 = do
    zxs <- shuffleList zxs0
    let grdcrc = loopCircuit' mlkl0 $ proc (zxs',mlkl) -> do
            let dmlkl = vanillaGradient $ conditionalLogLikelihoodDifferential zxs' mlkl
            gradientCircuit eps defaultAdamPursuit -< (mlkl,dmlkl)
    return . runIdentity . iterateCircuit grdcrc . breakEvery nbtch $ zxs

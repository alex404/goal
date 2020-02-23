{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    RankNTypes,
    TypeOperators,
    FlexibleContexts,
    ScopedTypeVariables
#-}
-- | The main module of goal-probability. Import this module to use all the
-- types, functions, and classes provided by goal-probability.
module Goal.Probability
    ( -- * Package Exports
      module Goal.Probability.Statistical
    , module Goal.Probability.ExponentialFamily
    , module Goal.Probability.Distributions
    , module Goal.Probability.LatentVariable
    , module Goal.Probability.ExponentialFamily.Harmonium
    , module Goal.Probability.ExponentialFamily.Conditional
    , module Goal.Probability.ExponentialFamily.Harmonium.Learning
    , module Goal.Probability.ExponentialFamily.Harmonium.Inference
      -- * Stochastic Operations
    , shuffleList
    , resampleVector
    , subsampleVector
    , noisyFunction
    -- ** Circuits
    , minibatcher
    -- * Statistics
    , estimateMeanVariance
    , estimateFanoFactor
    , estimateCoefficientOfVariation
    , estimateCorrelations
    , histograms
    -- ** Model Selection
    , akaikesInformationCriterion
    , bayesianInformationCriterion
    , conditionalAkaikesInformationCriterion
    , conditionalBayesianInformationCriterion
    -- * External Exports
    , module System.Random.MWC
    , module System.Random.MWC.Probability
    ) where


--- Imports ---


-- Re-exports --

import System.Random.MWC (Seed,save,restore)
import qualified System.Random.MWC as MWC

import System.Random.MWC.Probability hiding (initialize,sample)
import System.Random.MWC.Distributions (uniformShuffle)

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.LatentVariable
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Conditional
import Goal.Probability.ExponentialFamily.Harmonium.Learning
import Goal.Probability.ExponentialFamily.Harmonium.Inference


-- Package --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Mutable as M
import qualified Goal.Core.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable.Base as MV
import qualified Data.Vector as V

import qualified Statistics.Sample as STAT hiding (range)
import qualified Statistics.Sample.Histogram as STAT
import qualified Data.Vector.Storable as VS


--- Statistics ---


-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanVariance
    :: Traversable f
    => f Double
    -> (Double,Double)
{-# INLINE estimateMeanVariance #-}
estimateMeanVariance xs = STAT.meanVarianceUnb . VS.fromList $ toList xs

-- | Estimate the Fano Factor of a sample.
estimateFanoFactor
    :: Traversable f
    => f Double
    -> Double
{-# INLINE estimateFanoFactor #-}
estimateFanoFactor xs =
    let (mu,vr) = estimateMeanVariance xs
     in vr / mu

-- | Estimate the coefficient of variation from a sample.
estimateCoefficientOfVariation :: Traversable f => f Double -> Double
{-# INLINE estimateCoefficientOfVariation #-}
estimateCoefficientOfVariation zs =
    let (mu,vr) = estimateMeanVariance zs
     in sqrt vr / mu

-- | Computes the empirical covariance matrix given a sample if iid random vectors.
estimateCorrelations
    :: forall k x v . (G.VectorClass v x, G.VectorClass v Double, KnownNat k, Real x)
    => [G.Vector v k x]
    -> S.Matrix k k Double
{-# INLINE estimateCorrelations #-}
estimateCorrelations zs =
    let mnrm :: Source # MultivariateNormal k
        mnrm = mle $ G.convert . G.map realToFrac <$> zs
     in multivariateNormalCorrelations mnrm

---- | Stimulus Dependent Noise Correlations, ordered by preferred stimulus.
--mixturePopulationNoiseCorrelations
--    :: forall k n x . ( KnownNat k, KnownNat n, ExponentialFamily x )
--    => Mean #> Natural # MixtureGLM n (Neurons k) x -- ^ Mixture Encoder
--    -> SamplePoint x
--    -> S.Matrix k k Double -- ^ Mean Parameter Correlations
--{-# INLINE mixturePopulationNoiseCorrelations #-}
--mixturePopulationNoiseCorrelations mlkl x =
--    let mxmdl = mlkl >.>* x
--        (ngnss, nwghts) = splitMixtureModel mxmdl
--        wghts0 = coordinates $ toSource nwghts
--        wghts = 1 - S.sum wghts0 : S.toList wghts0
--        gnss = toMean <$> S.toList ngnss
--        mgns = weightedAveragePoint $ zip wghts gnss
--        mgns2 :: Natural #> Mean # Tensor (Neurons k) (Neurons k)
--        mgns2 = mgns >.< mgns
--        mmgns = weightedAveragePoint $ zip wghts [ gns >.< gns | gns <- gnss ]
--        cvgns = mmgns <-> mgns2
--        cvnrns = cvgns <+> (fromMatrix . S.diagonalMatrix $ coordinates mgns)
--        sdnrns = S.map sqrt $ S.takeDiagonal (toMatrix cvgns) `S.add` coordinates mgns
--        sdmtx = S.outerProduct sdnrns sdnrns
--     in G.Matrix $ S.zipWith (/) (G.toVector $ toMatrix cvnrns) (G.toVector sdmtx)

-- | Computes histograms (and densities) with the given number of bins for the
-- given list of samples. Bounds can be given or computed automatically. The
-- returned values are the list of bin centres and the binned samples. If bounds
-- are given but are not greater than all given sample points, then an error
-- will be thrown.
histograms
    :: Int -- ^ Number of Bins
    -> Maybe (Double, Double) -- ^ Maybe bin bounds
    -> [[Double]] -- ^ Datasets
    -> ([Double],[[Int]],[[Double]]) -- ^ Bin centres, counts, and densities for each dataset
{-# INLINE histograms #-}
histograms nbns mmnmx smpss =
    let (mn,mx) = case mmnmx of
                    Just (mn0,mx0) -> (mn0,mx0)
                    Nothing -> STAT.range nbns . VS.fromList $ concat smpss
        stp = (mx - mn) / fromIntegral nbns
        bns = take nbns [ mn + stp/2 + stp * fromIntegral n | n <- [0 :: Int,1..] ]
        hsts = VS.toList . STAT.histogram_ nbns mn mx . VS.fromList <$> smpss
        ttls = sum <$> hsts
        dnss = do
            (hst,ttl) <- zip hsts ttls
            return $ if ttl == 0
                        then []
                        else (/(fromIntegral ttl * stp)) . fromIntegral <$> hst
     in (bns,hsts,dnss)


--- Stochastic Functions ---


-- | Shuffle the elements of a list.
shuffleList :: [a] -> Random r [a]
{-# INLINE shuffleList #-}
shuffleList xs = fmap V.toList . Prob $ uniformShuffle (V.fromList xs)

minibatcher :: Int -> [x] -> Chain (Random r) [x]
{-# INLINE minibatcher #-}
minibatcher nbtch xs0 = accumulateFunction [] $ \() xs ->
    if (length xs < nbtch)
       then do
           xs1 <- shuffleList xs0
           let (hds',tls') = splitAt nbtch (xs ++ xs1)
           return (hds',tls')
       else do
           let (hds',tls') = splitAt nbtch xs
           return (hds',tls')

-- My idea here only works if all the conditions are the same size, which isn't
-- generally a good assmption. So maybe this is worthless? Unless of course
-- sometimes it's not...
--conditionalMinibatcher :: Int -> [(z,x)] -> Chain (Random r) (M.Map x [y])
--{-# INLINE conditionalMinibatcher #-}
--conditionalMinibatcher nbtch zxs0 =
--    let xzmp = conditionalDataMap zxs0
--     in accumulateFunction [] $ \() xs ->
--         if (length xs < nbtch)
--           then do
--               xs1 <- shuffleList xs0
--               let (hds',tls') = splitAt nbtch (xs ++ xs1)
--               return (hds',tls')
--           else do
--               let (hds',tls') = splitAt nbtch xs
--               return (hds',tls')





-- | Returns a uniform sample of elements from the given vector with replacement.
resampleVector :: (KnownNat n, KnownNat k) => B.Vector n x -> Random s (B.Vector k x)
{-# INLINE resampleVector #-}
resampleVector xs = do
    ks <- B.replicateM $ uniformR (0, B.length xs-1)
    return $ B.backpermute xs ks

-- | Returns a sample from the given function with added noise.
noisyFunction
    :: (Generative c x, Num (SamplePoint x))
    => Point c x -- ^ Noise model
    -> (y -> SamplePoint x) -- ^ Function
    -> y -- ^ Input
    -> Random s (SamplePoint x) -- ^ Stochastic Output
{-# INLINE noisyFunction #-}
noisyFunction m f x = do
    ns <- samplePoint m
    return $ f x + ns

-- | Take a random, unordered subset of a list.
subsampleVector
    :: forall k m v x r . (KnownNat k, KnownNat m, G.VectorClass v x)
    => G.Vector v (k + m) x
    -> Random r (G.Vector v k x)
{-# INLINE subsampleVector #-}
subsampleVector v = Prob $ \gn -> do
    let k = natValInt (Proxy :: Proxy k)
    mv <- G.thaw v
    randomSubSample0 k mv gn
    v' <- G.unsafeFreeze mv
    let foo :: (G.Vector v k x, G.Vector v m x)
        foo = G.splitAt v'
    return $ fst foo

randomSubSample0
    :: (KnownNat n, PrimMonad m, MV.MVector v a)
    => Int -> G.MVector v n (PrimState m) a -> Gen (PrimState m) -> m ()
{-# INLINE randomSubSample0 #-}
randomSubSample0 k v gn = looper 0
    where n = M.length v
          looper i
            | i == k = return ()
            | otherwise = do
                j <- MWC.uniformR (i,n-1) gn
                M.unsafeSwap v i j
                looper (i+1)


-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion
    :: forall c x . (Manifold x, AbsolutelyContinuous c x)
    => c # x
    -> Sample x
    -> Double
{-# INLINE akaikesInformationCriterion #-}
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
     in 2 * fromIntegral d - 2 * sum (log <$> densities p xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion
    :: forall c x . (AbsolutelyContinuous c x, Manifold x)
    => c # x
    -> Sample x
    -> Double
{-# INLINE bayesianInformationCriterion #-}
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum (log <$> densities p xs)

-- | Calculate the conditional AIC for a given model and sample.
conditionalAkaikesInformationCriterion
    :: forall d f x y
    . (AbsolutelyContinuous d y, ExponentialFamily x, Map Mean d f y x)
    => Function Mean d # f y x
    -> Sample (y,x)
    -> Double
{-# INLINE conditionalAkaikesInformationCriterion #-}
conditionalAkaikesInformationCriterion f yxs =
    let (ys,xs) = unzip yxs
        d = natVal (Proxy :: Proxy (Dimension y))
        yhts = f >$>* xs
     in 2 * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]

-- | Calculate the conditional BIC for a given model and sample.
conditionalBayesianInformationCriterion
    :: forall d f x y
    . (AbsolutelyContinuous d y, ExponentialFamily x, Map Mean d f y x)
    => Function Mean d # f y x
    -> Sample (y,x)
    -> Double
{-# INLINE conditionalBayesianInformationCriterion #-}
conditionalBayesianInformationCriterion f yxs =
    let (ys,xs) = unzip yxs
        d = natVal (Proxy :: Proxy (Dimension y))
        yhts = f >$>* xs
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]



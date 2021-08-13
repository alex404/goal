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
    , module Goal.Probability.Conditional
    , module Goal.Probability.Distributions
    , module Goal.Probability.Distributions.Gaussian
    , module Goal.Probability.Distributions.CoMPoisson
      -- * Stochastic Operations
    , shuffleList
    , resampleVector
    , subsampleVector
    , noisyFunction
    -- ** Circuits
    , minibatcher
    -- * Statistics
    , estimateMeanVariance
    , estimateMeanSD
    , estimateFanoFactor
    , estimateCoefficientOfVariation
    , estimateCorrelation
    , estimateCorrelations
    , histograms
    -- ** Model Selection
    , akaikesInformationCriterion
    , bayesianInformationCriterion
    --, conditionalAkaikesInformationCriterion
    --, conditionalBayesianInformationCriterion
    ) where


--- Imports ---


-- Re-exports --

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Conditional
import Goal.Probability.Distributions
import Goal.Probability.Distributions.Gaussian
import Goal.Probability.Distributions.CoMPoisson

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

import qualified System.Random.MWC as R
import qualified System.Random.MWC.Distributions as R


--- Statistics ---


-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanVariance
    :: Traversable f
    => f Double
    -> (Double,Double)
estimateMeanVariance xs = STAT.meanVarianceUnb . VS.fromList $ toList xs

-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanSD
    :: Traversable f
    => f Double
    -> (Double,Double)
estimateMeanSD xs =
    let (mu,vr) = estimateMeanVariance xs
     in (mu,sqrt vr)

-- | Estimate the Fano Factor of a sample.
estimateFanoFactor
    :: Traversable f
    => f Double
    -> Double
estimateFanoFactor xs =
    let (mu,vr) = estimateMeanVariance xs
     in vr / mu

-- | Estimate the coefficient of variation from a sample.
estimateCoefficientOfVariation :: Traversable f => f Double -> Double
estimateCoefficientOfVariation zs =
    let (mu,vr) = estimateMeanVariance zs
     in sqrt vr / mu

-- | Computes the empirical covariance matrix given a sample if iid random vectors.
estimateCorrelations
    :: forall k x v . (G.VectorClass v x, G.VectorClass v Double, KnownNat k, Real x)
    => [G.Vector v k x]
    -> S.Matrix k k Double
estimateCorrelations zs =
    let mnrm :: Source # MultivariateNormal k
        mnrm = mle $ G.convert . G.map realToFrac <$> zs
     in multivariateNormalCorrelations mnrm

-- | Computes the empirical covariance matrix given a sample from a bivariate random variable.
estimateCorrelation
    :: [(Double,Double)]
    -> Double
estimateCorrelation zs = STAT.correlation $ V.fromList zs

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
shuffleList :: [a] -> Random [a]
shuffleList xs = V.toList <$> Random (R.uniformShuffle (V.fromList xs))

-- | A 'Circuit' that helps fitting data based on minibatches. Essentially, it
-- creates an infinite list out of shuffled versions of the input list, and
-- breaks down and returns the result in chunks of the specified size.
minibatcher :: Int -> [x] -> Chain Random [x]
minibatcher nbtch xs0 = accumulateFunction [] $ \() xs ->
    if length (take nbtch xs) < nbtch
       then do
           xs1 <- shuffleList xs0
           let (hds',tls') = splitAt nbtch (xs ++ xs1)
           return (hds',tls')
       else do
           let (hds',tls') = splitAt nbtch xs
           return (hds',tls')

-- | Returns a uniform sample of elements from the given vector with replacement.
resampleVector :: (KnownNat n, KnownNat k) => B.Vector n x -> Random (B.Vector k x)
resampleVector xs = do
    ks <- B.replicateM $ Random (R.uniformR (0, B.length xs-1))
    return $ B.backpermute xs ks

-- | Returns a sample from the given function with added noise.
noisyFunction
    :: (Generative c x, Num (SamplePoint x))
    => Point c x -- ^ Noise model
    -> (y -> SamplePoint x) -- ^ Function
    -> y -- ^ Input
    -> Random (SamplePoint x) -- ^ Stochastic Output
noisyFunction m f x = do
    ns <- samplePoint m
    return $ f x + ns

-- | Take a random, unordered subset of a list.
subsampleVector
    :: forall k m v x . (KnownNat k, KnownNat m, G.VectorClass v x)
    => G.Vector v (k + m) x
    -> Random (G.Vector v k x)
subsampleVector v = Random $ \gn -> do
    let k = natValInt (Proxy :: Proxy k)
    mv <- G.thaw v
    randomSubSample0 k mv gn
    v' <- G.unsafeFreeze mv
    let foo :: (G.Vector v k x, G.Vector v m x)
        foo = G.splitAt v'
    return $ fst foo

randomSubSample0
    :: (KnownNat n, PrimMonad m, MV.MVector v a)
    => Int -> G.MVector v n (PrimState m) a -> R.Gen (PrimState m) -> m ()
randomSubSample0 k v gn = looper 0
    where n = M.length v
          looper i
            | i == k = return ()
            | otherwise = do
                j <- R.uniformR (i,n-1) gn
                M.unsafeSwap v i j
                looper (i+1)


-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion
    :: forall c x s . (Manifold x, LogLikelihood c x s)
    => c # x
    -> [s]
    -> Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
     in 2 * fromIntegral d - 2 * logLikelihood xs p * fromIntegral (length xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion
    :: forall c x s . (LogLikelihood c x s, Manifold x)
    => c # x
    -> [s]
    -> Double
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * logLikelihood xs p * fromIntegral (length xs)

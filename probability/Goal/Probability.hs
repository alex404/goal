{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | The main module of goal-probability. Import this module to use all the
types, functions, and classes provided by goal-probability.
-}
module Goal.Probability (
    -- * Package Exports
    module Goal.Probability.Statistical,
    module Goal.Probability.ExponentialFamily,
    module Goal.Probability.Conditional,
    module Goal.Probability.Distributions,
    module Goal.Probability.Distributions.Gaussian,
    module Goal.Probability.Distributions.CoMPoisson,

    -- * Stochastic Operations
    shuffleList,
    resampleVector,
    subsampleVector,
    noisyFunction,

    -- ** Circuits
    minibatcher,

    -- * Statistics
    estimateMeanVariance,
    estimateMeanSD,
    estimateFanoFactor,
    estimateCoefficientOfVariation,
    estimateCorrelation,
    estimateCorrelations,

    -- ** Model Selection
    akaikesInformationCriterion,
    bayesianInformationCriterion,
    -- , conditionalAkaikesInformationCriterion
    -- , conditionalBayesianInformationCriterion
) where

--- Imports ---

--- Re-exports

import Goal.Probability.Conditional
import Goal.Probability.Distributions
import Goal.Probability.Distributions.CoMPoisson
import Goal.Probability.Distributions.Gaussian
import Goal.Probability.ExponentialFamily
import Goal.Probability.Statistical

--- Goal

import Goal.Core
import Goal.Geometry

import Goal.Core.Vector.Boxed qualified as B
import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Generic.Mutable qualified as M
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Vector qualified as V
import Data.Vector.Generic.Mutable (PrimMonad, PrimState)
import Data.Vector.Generic.Mutable.Base qualified as MV
import Data.Vector.Storable qualified as VS

import Statistics.Sample qualified as STAT hiding (range)

import System.Random.MWC qualified as R
import System.Random.MWC.Distributions qualified as R

import Data.Foldable (toList)
import Data.Proxy

--- Statistics ---

-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanVariance ::
    (Traversable f) =>
    f Double ->
    (Double, Double)
estimateMeanVariance xs = STAT.meanVarianceUnb . VS.fromList $ toList xs

-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanSD ::
    (Traversable f) =>
    f Double ->
    (Double, Double)
estimateMeanSD xs =
    let (mu, vr) = estimateMeanVariance xs
     in (mu, sqrt vr)

-- | Estimate the Fano Factor of a sample.
estimateFanoFactor ::
    (Traversable f) =>
    f Double ->
    Double
estimateFanoFactor xs =
    let (mu, vr) = estimateMeanVariance xs
     in vr / mu

-- | Estimate the coefficient of variation from a sample.
estimateCoefficientOfVariation :: (Traversable f) => f Double -> Double
estimateCoefficientOfVariation zs =
    let (mu, vr) = estimateMeanVariance zs
     in sqrt vr / mu

-- | Computes the empirical covariance matrix given a sample if iid random vectors.
estimateCorrelations ::
    forall k x v.
    ( G.VectorClass v x
    , G.VectorClass v Double
    , KnownCovariance L.PositiveDefinite k
    , Real x
    ) =>
    [G.Vector v k x] ->
    Source # Tensor (StandardNormal k) (StandardNormal k)
estimateCorrelations zs =
    let mnrm :: Source # FullNormal k
        mnrm = mle $ G.convert . G.map realToFrac <$> zs
     in multivariateNormalCorrelations mnrm

-- | Computes the empirical covariance matrix given a sample from a bivariate random variable.
estimateCorrelation ::
    [(Double, Double)] ->
    Double
estimateCorrelation zs = STAT.correlation $ V.fromList zs

--- Stochastic Functions ---

-- | Shuffle the elements of a list.
shuffleList :: [a] -> Random [a]
shuffleList xs = V.toList <$> Random (R.uniformShuffle (V.fromList xs))

{- | A 'Circuit' that helps fitting data based on minibatches. Essentially, it
creates an infinite list out of shuffled versions of the input list, and
breaks down and returns the result in chunks of the specified size.
-}
minibatcher :: Int -> [x] -> Chain Random [x]
minibatcher nbtch xs0 = accumulateFunction [] $ \() xs ->
    if length (take nbtch xs) < nbtch
        then do
            xs1 <- shuffleList xs0
            let (hds', tls') = splitAt nbtch (xs ++ xs1)
            return (hds', tls')
        else do
            let (hds', tls') = splitAt nbtch xs
            return (hds', tls')

-- | Returns a uniform sample of elements from the given vector with replacement.
resampleVector :: (KnownNat n, KnownNat k) => B.Vector n x -> Random (B.Vector k x)
resampleVector xs = do
    ks <- B.replicateM $ Random (R.uniformR (0, B.length xs - 1))
    return $ B.backpermute xs ks

-- | Returns a sample from the given function with added noise.
noisyFunction ::
    (Generative c x, Num (SamplePoint x)) =>
    -- | Noise model
    Point c x ->
    -- | Function
    (y -> SamplePoint x) ->
    -- | Input
    y ->
    -- | Stochastic Output
    Random (SamplePoint x)
noisyFunction m f x = do
    ns <- samplePoint m
    return $ f x + ns

-- | Take a random, unordered subset of a list.
subsampleVector ::
    forall k m v x.
    (KnownNat k, KnownNat m, G.VectorClass v x) =>
    G.Vector v (k + m) x ->
    Random (G.Vector v k x)
subsampleVector v = Random $ \gn -> do
    let k = natValInt (Proxy :: Proxy k)
    mv <- G.thaw v
    randomSubSample0 k mv gn
    v' <- G.unsafeFreeze mv
    let foo :: (G.Vector v k x, G.Vector v m x)
        foo = G.splitAt v'
    return $ fst foo

randomSubSample0 ::
    (KnownNat n, PrimMonad m, MV.MVector v a) =>
    Int ->
    G.MVector v n (PrimState m) a ->
    R.Gen (PrimState m) ->
    m ()
randomSubSample0 k v gn = looper 0
  where
    n = M.length v
    looper i
        | i == k = return ()
        | otherwise = do
            j <- R.uniformR (i, n - 1) gn
            M.unsafeSwap v i j
            looper (i + 1)

-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion ::
    forall c x s.
    (Manifold x, LogLikelihood c x s) =>
    c # x ->
    [s] ->
    Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
     in 2 * fromIntegral d - 2 * logLikelihood xs p * fromIntegral (length xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion ::
    forall c x s.
    (LogLikelihood c x s, Manifold x) =>
    c # x ->
    [s] ->
    Double
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * logLikelihood xs p * fromIntegral (length xs)

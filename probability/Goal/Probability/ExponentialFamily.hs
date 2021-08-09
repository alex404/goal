{-# LANGUAGE UndecidableInstances,TypeApplications #-}
-- | Definitions for working with exponential families.
module Goal.Probability.ExponentialFamily
    ( -- * Exponential Families
    ExponentialFamily (sufficientStatistic, averageSufficientStatistic, logBaseMeasure)
    , LegendreExponentialFamily
    , DuallyFlatExponentialFamily
    , exponentialFamilyLogDensities
    , unnormalizedLogDensities
    -- ** Coordinate Systems
    , Natural
    , Mean
    , Source
    -- ** Coordinate Transforms
    , toNatural
    , toMean
    , toSource
    -- ** Entropies
    , relativeEntropy
    , crossEntropy
    -- ** Differentials
    , relativeEntropyDifferential
    , stochasticRelativeEntropyDifferential
    , stochasticInformationProjectionDifferential
    -- *** Maximum Likelihood Instances
    , exponentialFamilyLogLikelihood
    , exponentialFamilyLogLikelihoodDifferential
    ) where

--- Imports ---


-- Package --

import Goal.Probability.Statistical

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Foreign.Storable

--- Exponential Families ---


-- | A parameterization which represents the standard or typical parameterization of
-- the given manifold, e.g. the Poisson rate or Normal mean and standard deviation.
data Source

-- | A parameterization in terms of the natural parameters of an exponential family.
data Natural

-- | A parameterization in terms of the mean 'sufficientStatistic' of an exponential family.
data Mean

instance Primal Natural where
    type Dual Natural = Mean

instance Primal Mean where
    type Dual Mean = Natural

-- | Expresses an exponential family distribution in 'Natural' coordinates.
toNatural :: (Transition c Natural x) => c # x -> Natural # x
toNatural = transition

-- | Expresses an exponential family distribution in 'Mean' coordinates.
toMean :: (Transition c Mean x) => c # x -> Mean # x
toMean = transition

-- | Expresses an exponential family distribution in 'Source' coordinates.
toSource :: (Transition c Source x) => c # x -> Source # x
toSource = transition

-- | An 'ExponentialFamily' is a 'Statistical' 'Manifold' \( \mathcal M \)
-- determined by a fixed-length 'sufficientStatistic' \(s_i\) and a
-- 'logBaseMeasure' \(\mu\). Each distribution \(P \in \mathcal M\) may then be
-- identified with 'Natural' parameters \(\theta_i\) such that
-- \(p(x) \propto e^{\sum_{i=1}^n \theta_i s_i(x)}\mu(x)\).  'ExponentialFamily'
-- distributions theoretically have a 'Riemannian' geometry, with 'metric'
-- 'Tensor' given by the Fisher information metric. However, not all
-- distributions (e.g. the von Mises distribution) afford closed-form
-- expressions for all the relevant structures.
class Statistical x => ExponentialFamily x where
    sufficientStatistic :: SamplePoint x -> Mean # x
    averageSufficientStatistic :: Sample x -> Mean # x
    averageSufficientStatistic = average . map sufficientStatistic
    logBaseMeasure :: Proxy x -> SamplePoint x -> Double

-- | When the log-partition function and its derivative of the given
-- 'ExponentialFamily' may be computed in closed-form, then we refer to it as a
-- 'LegendreExponentialFamily'.
--
-- Note that the log-partition function is the 'potential' of the 'Legendre'
-- class, and its derivative maps 'Natural' coordinates to 'Mean' coordinates.
type LegendreExponentialFamily x =
    ( PotentialCoordinates x ~ Natural, Legendre x, ExponentialFamily x
    , Transition (PotentialCoordinates x) (Dual (PotentialCoordinates x)) x )

-- | When additionally, the (negative) entropy and its derivative of the given
-- 'ExponentialFamily' may be computed in closed-form, then we refer to it as a
-- 'DuallyFlatExponentialFamily'.
--
-- Note that the negative entropy is the 'dualPotential' of the 'DuallyFlat' class,
-- and its derivative maps 'Mean' coordinates to 'Natural' coordinates.
type DuallyFlatExponentialFamily x =
    ( LegendreExponentialFamily x, DuallyFlat x
    , Transition (Dual (PotentialCoordinates x)) (PotentialCoordinates x) x )

-- | The relative entropy \(D(P \parallel Q)\), also known as the KL-divergence.
-- This is simply the 'canonicalDivergence' with its arguments flipped.
relativeEntropy :: DuallyFlatExponentialFamily x => Mean # x -> Natural # x -> Double
relativeEntropy = flip canonicalDivergence

-- | A function for computing the cross-entropy, which is the relative entropy
-- plus the entropy of the first distribution.
crossEntropy :: DuallyFlatExponentialFamily x => Mean # x -> Natural # x ->
    Double
crossEntropy mp nq = potential nq - (mp <.> nq)

-- | The differential of the relative entropy with respect to the 'Natural' parameters of
-- the second argument.
relativeEntropyDifferential :: LegendreExponentialFamily x => Mean # x -> Natural # x -> Mean # x
relativeEntropyDifferential mp nq = transition nq - mp

-- | Monte Carlo estimate of the differential of the relative entropy with
-- respect to the 'Natural' parameters of the second argument, based on samples from
-- the two distributions.
stochasticRelativeEntropyDifferential
    :: ExponentialFamily x
    => Sample x -- ^ True Samples
    -> Sample x -- ^ Model Samples
    -> Mean # x -- ^ Differential Estimate
stochasticRelativeEntropyDifferential pxs qxs =
    averageSufficientStatistic qxs - averageSufficientStatistic pxs

-- | Estimate of the differential of relative entropy with respect to the
-- 'Natural' parameters of the first argument, based a 'Sample' from the first
-- argument and the unnormalized log-density of the second.
stochasticInformationProjectionDifferential
    :: ExponentialFamily x
    => Natural # x -- ^ Model Distribution
    -> Sample x -- ^ Model Samples
    -> (SamplePoint x -> Double) -- ^ Unnormalized log-density of target distribution
    -> Mean # x -- ^ Differential Estimate
stochasticInformationProjectionDifferential px xs f =
    let mxs = sufficientStatistic <$> xs
        mys = (\x -> sufficientStatistic x <.> px - f x) <$> xs
        ln = fromIntegral $ length xs
        mxht = ln /> sum mxs
        myht = sum mys / ln
     in (ln - 1) /> sum [ (my - myht) .> (mx - mxht) | (mx,my) <- zip mxs mys ]

-- | The density of an exponential family distribution that has an exact
-- expression for the log-partition function.
exponentialFamilyLogDensities
    :: (ExponentialFamily x, Legendre x, PotentialCoordinates x ~ Natural) => Natural # x -> Sample x -> [Double]
exponentialFamilyLogDensities p xs = subtract (potential p) <$> unnormalizedLogDensities p xs

-- | The unnormalized log-density of an arbitrary exponential family distribution.
unnormalizedLogDensities :: forall x . ExponentialFamily x => Natural # x -> Sample x -> [Double]
unnormalizedLogDensities p xs =
    zipWith (+) (dotMap p $ sufficientStatistic <$> xs) (logBaseMeasure (Proxy @ x) <$> xs)

-- | 'logLikelihood' for a 'LegendreExponentialFamily'.
exponentialFamilyLogLikelihood
    :: forall x . LegendreExponentialFamily x
    => Sample x -> Natural # x -> Double
exponentialFamilyLogLikelihood xs nq =
    let mp = averageSufficientStatistic xs
        bm = average $ logBaseMeasure (Proxy :: Proxy x) <$> xs
     in -potential nq + (mp <.> nq) + bm

-- | 'logLikelihoodDifferential' for a 'LegendreExponentialFamily'.
exponentialFamilyLogLikelihoodDifferential
    :: LegendreExponentialFamily x
    => Sample x -> Natural # x -> Mean # x
exponentialFamilyLogLikelihoodDifferential xs nq =
    let mp = averageSufficientStatistic xs
     in mp - transition nq


--- Internal ---


replicatedlogBaseMeasure0 :: (ExponentialFamily x, Storable (SamplePoint x), KnownNat k)
                       => Proxy x -> Proxy (Replicated k x) -> S.Vector k (SamplePoint x) -> Double
replicatedlogBaseMeasure0 prxym _ xs = S.sum $ S.map (logBaseMeasure prxym) xs

pairlogBaseMeasure
    :: (ExponentialFamily x, ExponentialFamily y)
    => Proxy x
    -> Proxy y
    -> Proxy (x,y)
    -> SamplePoint (x,y)
    -> Double
pairlogBaseMeasure prxym prxyn _ (xm,xn) =
     logBaseMeasure prxym xm + logBaseMeasure prxyn xn


--- Instances ---


-- Replicated --

instance Transition Natural Natural x where
    transition = id

instance Transition Mean Mean x where
    transition = id

instance Transition Source Source x where
    transition = id

instance (ExponentialFamily x, Storable (SamplePoint x), KnownNat k)
  => ExponentialFamily (Replicated k x) where
    sufficientStatistic xs = joinReplicated $ S.map sufficientStatistic xs
    logBaseMeasure = replicatedlogBaseMeasure0 Proxy

-- Sum --

instance (ExponentialFamily x, ExponentialFamily y) => ExponentialFamily (x,y) where
    sufficientStatistic (xm,xn) =
         join (sufficientStatistic xm) (sufficientStatistic xn)
    logBaseMeasure = pairlogBaseMeasure Proxy Proxy


-- Source Coordinates --

instance Primal Source where
    type Dual Source = Source

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Implementation of Conway-Maxwell Poisson distributions (CoMPoisson).
-- (<https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9876.2005.00474.x>) CoMPoisson distributions generalize Poisson distributions with
-- a shape parameter that can concentrate or disperse the underlying Poisson
-- distribution.
module Goal.Probability.Distributions.CoMPoisson
    (
    -- * CoMPoisson
      CoMPoisson
    , CoMShape
    -- ** Numerics
    , comPoissonLogPartitionSum
    , comPoissonExpectations
    ) where

-- Package --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import System.Random.MWC.Probability


--- Analysis ---

--- CoMPoisson Distribution ---

-- | A type for storing the shape of a 'CoMPoisson' distribution.
data CoMShape

-- | The 'Manifold' of 'CoMPoisson' distributions. The 'Source' coordinates of the
-- 'CoMPoisson' are the mode $\mu$ and the "pseudo-precision" parameter $\nu$, such that $\mu / \nu$ is approximately the variance of the distribution.
type CoMPoisson = LocationShape Poisson CoMShape

-- | Approximates the log-partition function of the given CoMPoisson
-- distribution up to the specified precision.
comPoissonLogPartitionSum :: Double -> Natural # CoMPoisson -> Double
{-# INLINE comPoissonLogPartitionSum #-}
comPoissonLogPartitionSum eps np =
    let (tht1,tht2) = S.toPair $ coordinates np
     in fst $ comPoissonLogPartitionSum0 eps tht1 tht2

-- | Approximates the expectations of functions given the natural parameters of
-- a CoM-Poisson distribution.
comPoissonExpectations
    :: KnownNat n
    => Double
    -> (Int -> S.Vector n Double)
    -> Natural # CoMPoisson
    -> S.Vector n Double
{-# INLINE comPoissonExpectations #-}
comPoissonExpectations eps f np =
    let (tht1,tht2) = S.toPair $ coordinates np
        (lgprt,ln) = comPoissonLogPartitionSum0 eps tht1 tht2
        js = [0..ln]
        dns = exp . subtract lgprt <$> unnormalizedLogDensities np js
     in sum $ zipWith S.scale dns (f <$> js)

-- | Approximates the mean mparameters of a CoM-Poisson distribution.
comPoissonMeans :: Double -> Natural # CoMPoisson -> Mean # CoMPoisson
{-# INLINE comPoissonMeans #-}
comPoissonMeans eps cp =
    let ss :: Int -> Mean # CoMPoisson
        ss = sufficientStatistic
     in Point $ comPoissonExpectations eps (coordinates . ss) cp


--- Internal ---


comPoissonSequence :: Double -> Double -> [Double]
comPoissonSequence tht1 tht2 =
    [ tht1 * fromIntegral j + logFactorial j *tht2 | (j :: Int) <- [0..] ]

comPoissonLogPartitionSum0 :: Double -> Double -> Double -> (Double, Int)
{-# INLINE comPoissonLogPartitionSum0 #-}
comPoissonLogPartitionSum0 eps tht1 tht2 =
    let md = floor $ comPoissonSmoothMode tht1 tht2
        (hdsqs,tlsqs) = splitAt md $ comPoissonSequence tht1 tht2
        mx = tht1 * fromIntegral md + logFactorial md *tht2
        ehdsqs = exp . subtract mx <$> hdsqs
        etlsqs = exp . subtract mx <$> tlsqs
        sqs' = ehdsqs ++ takeWhile (> eps) etlsqs
     in ((+ mx) . log1p . subtract 1 $ sum sqs' , length sqs')

comPoissonSmoothMode :: Double -> Double -> Double
comPoissonSmoothMode tht1 tht2 = exp (tht1/negate tht2)

--comPoissonApproximateMean :: Double -> Double -> Double
--comPoissonApproximateMean mu nu =
--    mu + 1/(2*nu) - 0.5
--
--comPoissonApproximateVariance :: Double -> Double -> Double
--comPoissonApproximateVariance mu nu = mu / nu

overDispersedEnvelope :: Double -> Double -> Double -> Double
overDispersedEnvelope p mu nu =
    let mnm1 = 1 - p
        flrd = max 0 . floor $ mu / (mnm1**recip nu)
        nmr = mu**(nu * fromIntegral flrd)
        dmr = (mnm1^flrd) * (factorial flrd ** nu)
     in recip p * nmr / dmr

underDispersedEnvelope :: Double -> Double -> Double
underDispersedEnvelope mu nu =
    let fmu = floor mu
     in (mu ^ fmu / factorial fmu)** (nu - 1)

sampleOverDispersed :: Double -> Double -> Double -> Double -> Random r Int
sampleOverDispersed p bnd0 mu nu = do
    u0 <- uniform
    let y' = max 0 . floor $ logBase (1 - p) u0
        nmr = (mu^y' / factorial y')**nu
        dmr = bnd0 * (1-p)^y' * p
        alph = nmr/dmr
    u <- uniform
    if isNaN alph
       then error "NaN in sampling CoMPoisson: Parameters out of bounds"
       else if u <= alph
       then return y'
       else sampleOverDispersed p bnd0 mu nu

sampleUnderDispersed :: Double -> Double -> Double -> Random r Int
sampleUnderDispersed bnd0 mu nu = do
    let psn :: Source # Poisson
        psn = Point $ S.singleton mu
    y' <- samplePoint psn
    let alph0 = mu^y' / factorial y'
        alph = alph0**nu / (bnd0*alph0)
    u <- uniform
    if u <= alph
       then return y'
    else sampleUnderDispersed bnd0 mu nu

sampleCoMPoisson :: Int -> Double -> Double -> Random r [Int]
sampleCoMPoisson n mu nu
  | nu >= 1 =
      let bnd0 = underDispersedEnvelope mu nu
       in replicateM n $ sampleUnderDispersed bnd0 mu nu
  | otherwise =
      let p = 2*nu / (2*mu*nu + 1 + nu)
          bnd0 = overDispersedEnvelope p mu nu
       in replicateM n $ sampleOverDispersed p bnd0 mu nu


-- Instances --


instance ExponentialFamily CoMPoisson where
    sufficientStatistic k = fromTuple (fromIntegral k, logFactorial k)
    logBaseMeasure _ _ = 0

type instance PotentialCoordinates CoMPoisson = Natural

instance Legendre CoMPoisson where
    potential =
         comPoissonLogPartitionSum 1e-16

instance AbsolutelyContinuous Natural CoMPoisson where
    logDensities = exponentialFamilyLogDensities

instance Transition Source Natural CoMPoisson where
    transition p =
        let (mu,nu) = S.toPair $ coordinates p
         in fromTuple (nu * log mu, -nu)

instance Transition Natural Source CoMPoisson where
    transition p =
        let (tht1,tht2) = S.toPair $ coordinates p
         in fromTuple (exp (-tht1/tht2), -tht2)

instance (Transition c Source CoMPoisson) => Generative c CoMPoisson where
    sample n p = do
        let (mu,nu) = S.toPair . coordinates $ toSource p
         in sampleCoMPoisson n mu nu

instance Transition Natural Mean CoMPoisson where
    transition = comPoissonMeans 1e-16

instance Transition Source Mean CoMPoisson where
    transition = toMean . toNatural

instance LogLikelihood Natural CoMPoisson Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance Manifold CoMShape where
    type Dimension CoMShape = 1

instance Statistical CoMShape where
    type SamplePoint CoMShape = Int


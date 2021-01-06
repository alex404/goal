{-# LANGUAGE UndecidableInstances #-}

-- | Core types, classes, and functions for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical
    ( -- * Random
      Random
    , Statistical (SamplePoint)
    , Sample
    , realize
    -- * Initializiation
    , initialize
    , uniformInitialize
    , uniformInitialize'
    -- * Properties of Distributions
    , Generative (sample,samplePoint)
    , AbsolutelyContinuous (logDensity,logDensities,density,densities)
    , Discrete (Cardinality,sampleSpace)
    , pointSampleSpace
    , expectation
    -- ** Maximum Likelihood Estimation
    , MaximumLikelihood (mle)
    , LogLikelihood (logLikelihood,logLikelihoodDifferential)
    ) where


--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

-- Qualified --

import qualified System.Random.MWC.Probability as P
import qualified Control.Monad.ST as ST

import Foreign.Storable


--- Probability Theory ---


-- | A 'Manifold' is 'Statistical' if it is a set of probability distributions
-- over a particular sample space, where the sample space is a set of the
-- specified 'SamplePoint's.
class Manifold x => Statistical x where
    type SamplePoint x :: Type

-- | A 'Sample' is a list of 'SamplePoint's.
type Sample x = [SamplePoint x]

-- | A random variable.
type Random s = P.Prob (ST.ST s)

-- | Turn a random variable into an IO action.
realize :: Random s a -> IO a
realize = P.withSystemRandom . P.sample

-- | Probability distributions for which the sample space is countable. This
-- affords brute force computation of expectations.
class KnownNat (Cardinality x) => Discrete x where
    type Cardinality x :: Nat
    sampleSpace :: Proxy x -> Sample x

-- | Convenience function for getting the sample space of a 'Discrete'
-- probability distribution.
pointSampleSpace :: forall c x . Discrete x => c # x -> Sample x
pointSampleSpace _ = sampleSpace (Proxy :: Proxy x)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by @mwc-random@.
class Statistical x => Generative c x where
    samplePoint :: Point c x -> Random r (SamplePoint x)
    samplePoint = fmap head . sample 1
    sample :: Int -> Point c x -> Random r (Sample x)
    sample n = replicateM n . samplePoint


-- | The distributions \(P \in \mathcal M\) in a 'Statistical' 'Manifold'
-- \(\mathcal M\) are 'AbsolutelyContinuous' if there is a reference measure
-- \(\mu\) and a function \(p\) such that
-- \(P(A) = \int_A p d\mu\). We refer to \(p(x)\) as the 'density' of the
-- probability distribution.
class Statistical x => AbsolutelyContinuous c x where
    logDensity :: Point c x -> SamplePoint x -> Double
    logDensity p = head . logDensities p . (:[])
    logDensities :: Point c x -> Sample x -> [Double]
    logDensities p = map (logDensity p)

    density :: Point c x -> SamplePoint x -> Double
    density p = exp . head . logDensities p . (:[])
    densities :: Point c x -> Sample x -> [Double]
    densities p = map (exp . logDensity p)

-- | 'expectation' computes the brute force expected value of a 'Finite' set
-- given an appropriate 'density'.
expectation
    :: forall c x . (AbsolutelyContinuous c x, Discrete x)
    => Point c x
    -> (SamplePoint x -> Double)
    -> Double
expectation p f =
    let xs = sampleSpace (Proxy :: Proxy x)
     in sum $ zipWith (*) (f <$> xs) (densities p xs)

-- Maximum Likelihood Estimation

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical x => MaximumLikelihood c x where
    mle :: Sample x -> c # x

-- | Average log-likelihood and the differential for gradient ascent.
class Manifold x => LogLikelihood c x s where
    logLikelihood :: [s] -> c # x -> Double
    --logLikelihood xs p = average $ log <$> densities p xs
    logLikelihoodDifferential :: [s] -> c # x -> c #* x


--- Construction ---


-- | Generates a random point on the target 'Manifold' by generating random
-- samples from the given distribution.
initialize
    :: (Manifold x, Generative d y, SamplePoint y ~ Double)
    => d # y
    -> Random r (c # x)
initialize q = Point <$> S.replicateM (samplePoint q)

-- | Generates an initial point on the target 'Manifold' by generating uniform
-- samples from the given vector of bounds.
uniformInitialize' :: Manifold x => B.Vector (Dimension x) (Double,Double) -> Random r (Point c x)
uniformInitialize' bnds =
    Point . G.convert <$> mapM P.uniformR bnds

-- | Generates an initial point on the target 'Manifold' by generating uniform
-- samples from the given vector of bounds.
uniformInitialize :: Manifold x => (Double,Double) -> Random r (Point c x)
uniformInitialize bnds =
    Point <$> S.replicateM (P.uniformR bnds)


--- Instances ---


-- Replicated --

instance (Statistical x, KnownNat k, Storable (SamplePoint x))
  => Statistical (Replicated k x) where
    type SamplePoint (Replicated k x) = S.Vector k (SamplePoint x)

instance (KnownNat k, Generative c x, Storable (SamplePoint x))
  => Generative c (Replicated k x) where
    samplePoint = S.mapM samplePoint . splitReplicated

instance (KnownNat k, Storable (SamplePoint x), AbsolutelyContinuous c x)
  => AbsolutelyContinuous c (Replicated k x) where
    density cxs = S.product . S.zipWith density (splitReplicated cxs)

instance (KnownNat k, LogLikelihood c x s, Storable s)
  => LogLikelihood c (Replicated k x) (S.Vector k s) where
    logLikelihood cxs ps = S.sum . S.imap subLogLikelihood $ splitReplicated ps
        where subLogLikelihood fn = logLikelihood (flip S.index fn <$> cxs)
    logLikelihoodDifferential cxs ps =
        joinReplicated . S.imap subLogLikelihoodDifferential $ splitReplicated ps
            where subLogLikelihoodDifferential fn =
                    logLikelihoodDifferential (flip S.index fn <$> cxs)


-- Pair --


instance (Statistical x) => Statistical [x] where
    type SamplePoint [x] = [SamplePoint x]

instance (Statistical x, Statistical y)
  => Statistical (x,y) where
    type SamplePoint (x,y) = (SamplePoint x, SamplePoint y)

instance (Generative c x, Generative c y) => Generative c (x,y) where
    samplePoint pmn = do
        let (pm,pn) = splitPair pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c x, AbsolutelyContinuous c y)
  => AbsolutelyContinuous c (x,y) where
    density pmn (xm,xn) =
        let (pm,pn) = splitPair pmn
         in density pm xm * density pn xn


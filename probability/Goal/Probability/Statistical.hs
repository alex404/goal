{-# LANGUAGE UndecidableInstances,TypeFamilies,DataKinds,FlexibleContexts,MultiParamTypeClasses,FlexibleInstances #-}
-- | Here we provide the basic types and classes for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical (
    -- * Random
      Random
    , realize
    -- * Statistical Manifolds
    , Statistical (SamplePoint)
    , Sample
    -- * Construction
    , initialize
    -- * Properties of Distributions
    , Generative (samplePoint)
    , sample
    , AbsolutelyContinuous (density)
    , MaximumLikelihood (mle)
    , Discrete (Cardinality,sampleSpace)
    , expectation
    -- * Model Selection
    , akaikesInformationCriterion
    , bayesianInformationCriterion
    -- * Goodness of Fit
    -- , pseudoRSquared
    ) where


--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

-- Qualified --

import qualified System.Random.MWC.Probability as P
import qualified Control.Monad.ST as ST


--- Test Bed ---


--- Probability Theory ---

-- | A random variable.
type Random s a = P.Prob (ST.ST s) a

-- | Turn a random variable into an IO action.
realize :: Random s a -> IO a
realize = P.withSystemRandom . P.sample

-- | A 'Statistical' 'Manifold' is a 'Manifold' of probability distributions,
-- which all have in common a particular 'SampleSpace'.
class Manifold m => Statistical m where
    type SamplePoint m :: *

-- | A 'Vector' of 'SamplePoint's.
type Sample k m = B.Vector k (SamplePoint m)

class (KnownNat (Cardinality m), Statistical m) => Discrete m where
    type Cardinality m :: Nat
    sampleSpace :: Proxy m -> Sample (Cardinality m) m

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by MWC Monad.
class Statistical m => Generative c m where
    samplePoint :: Point c m -> Random r (SamplePoint m)

-- | Generates a 'Vector' of 'samplePoint's.
sample :: (Generative c m, KnownNat k) => Point c m -> Random r (Sample k m)
{-# INLINE sample #-}
sample = B.replicateM . samplePoint

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Statistical m => AbsolutelyContinuous c m where
    density :: Point c m -> SamplePoint m -> Double

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
expectation
    :: forall c m. (AbsolutelyContinuous c m, Discrete m)
    => Point c m
    -> (SamplePoint m -> Double)
    -> Double
expectation p f =
     B.sum $ (\x -> f x * density p x) <$> sampleSpace (Proxy :: Proxy m)

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical m => MaximumLikelihood c m where
    mle :: (KnownNat k, 1 <= k) => Sample k m -> Point c m


--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'dimension' m
-- samples from the given distribution.
initialize :: (Manifold m, Generative d n, SamplePoint n ~ Double) => d # n -> Random r (Point c m)
initialize q = fromBoxed <$> sample q


--- Model Evaluation ---


akaikesInformationCriterion
    :: forall c m k . (AbsolutelyContinuous c m, KnownNat k)
    => c # m
    -> Sample k m
    -> Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension m))
     in 2 * fromIntegral d - 2 * sum (log . density p <$> xs)

bayesianInformationCriterion
    :: forall c m k . (AbsolutelyContinuous c m, KnownNat k)
    => c # m
    -> Sample k m
    -> Double
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension m))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum (log . density p <$> xs)


--- Instances ---


instance (KnownNat k, Statistical m) => Statistical (Replicated k m) where
    type SamplePoint (Replicated k m) = B.Vector k (SamplePoint m)

instance (KnownNat k, Statistical m, Generative c m) => Generative c (Replicated k m) where
    {-# INLINE samplePoint #-}
    samplePoint = B.mapM samplePoint . G.convert . splitReplicated


instance (KnownNat k, Statistical m, AbsolutelyContinuous c m) => AbsolutelyContinuous c (Replicated k m) where
    {-# INLINE density #-}
    density ps xs = B.product $ B.zipWith density (G.convert $ splitReplicated ps) xs

instance (Statistical m, Statistical n) => Statistical (Sum m n) where
    type SamplePoint (Sum m n) = (SamplePoint m, SamplePoint n)

instance (Generative c m, Generative c n) => Generative c (Sum m n) where
    {-# INLINE samplePoint #-}
    samplePoint pmn = do
        let (pm,pn) = splitSum pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c m, AbsolutelyContinuous c n) => AbsolutelyContinuous c (Sum m n) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitSum pmn
         in density pm xm * density pn xn


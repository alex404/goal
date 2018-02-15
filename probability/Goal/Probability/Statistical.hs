{-# LANGUAGE TypeFamilies,DataKinds,FlexibleContexts,MultiParamTypeClasses,FlexibleInstances #-}
-- | Here we provide the basic types and classes for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical (
    -- * Statistical Manifolds
      Statistical (Sample)
    , Discrete (Cardinality,sampleSpace)
    -- ** Properties of Distributions
    , Generative (generate)
    , AbsolutelyContinuous (density)
    , expectation
    , MaximumLikelihood (mle)
    -- ** Construction
    , initialize
    -- * Random
    , Random
    , realize
    ) where


--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry

-- Qualified --

import qualified System.Random.MWC.Monad as R


--- Test Bed ---


--- Probability Theory ---

type Random s a = R.RandST s a

realize :: Random s a -> IO a
realize = R.runWithSystemRandom

-- | A 'Statistical' 'Manifold' is a 'Manifold' of probability distributions,
-- which all have in common a particular 'SampleSpace'.
class Manifold m => Statistical m where
    type Sample m :: *

class (KnownNat (Cardinality m), Statistical m) => Discrete m where
    type Cardinality m :: Nat
    sampleSpace :: Proxy m -> Vector (Cardinality m) (Sample m)

-- | A distribution is 'Generative' if we can 'generate' samples from it. Generation is
-- powered by MWC Monad.
class Statistical m => Generative c m where
    generate :: RealFloat x => Point c m x -> Random r (Sample m)

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Statistical m => AbsolutelyContinuous c m where
    density :: RealFloat x => Point c m x -> Sample m -> x

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
expectation0 :: (AbsolutelyContinuous c m, Discrete m, RealFloat x) => Proxy m -> Point c m x -> (Sample m -> x) -> x
expectation0 prxy p f =
    let xs = sampleSpace prxy
     in sum $ zipWithV (*) (f <$> xs) (density p <$> xs)

expectation :: (AbsolutelyContinuous c m, Discrete m, RealFloat x) => Point c m x -> (Sample m -> x) -> x
expectation = expectation0 Proxy


-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical m => MaximumLikelihood c m where
    mle :: (RealFloat x, Traversable f) => f (Sample m) -> Point c m x


--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'dimension' m
-- samples from the given distribution.
initialize :: (Manifold m, Generative d n, Sample n ~ Double, RealFloat x) => Point d n x -> Random r (Point c m x)
initialize q = do
    c0s <- replicateMV $ generate q
    return . Point $ realToFrac <$> c0s

instance (KnownNat k, Statistical m) => Statistical (Replicated k m) where
    type Sample (Replicated k m) = Vector k (Sample m)

instance (KnownNat k, Statistical m, Generative c m) => Generative c (Replicated k m) where
    {-# INLINE generate #-}
    generate = sequence . mapReplicated generate

instance (KnownNat k, Statistical m, AbsolutelyContinuous c m) => AbsolutelyContinuous c (Replicated k m) where
    {-# INLINE density #-}
    density ps xs = product $ zipWithV ($) (mapReplicated density ps) xs

instance (Statistical m, Statistical n) => Statistical (Sum m n) where
    type Sample (Sum m n) = (Sample m, Sample n)

instance (Generative c m, Generative c n) => Generative c (Sum m n) where
    {-# INLINE generate #-}
    generate pmn = do
        let (pm,pn) = splitSum pmn
        xm <- generate pm
        xn <- generate pn
        return (xm,xn)

instance (AbsolutelyContinuous c m, AbsolutelyContinuous c n) => AbsolutelyContinuous c (Sum m n) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitSum pmn
         in density pm xm * density pn xn


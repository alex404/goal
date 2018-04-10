{-# LANGUAGE UndecidableInstances,TypeFamilies,DataKinds,FlexibleContexts,MultiParamTypeClasses,FlexibleInstances #-}
-- | Here we provide the basic types and classes for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical (
    -- * Statistical Manifolds
      Statistical (SampleDimension)
    , Sample
    -- , Discrete (Cardinality,sampleSpace)
    -- ** Properties of Distributions
    , Generative (sample)
    , AbsolutelyContinuous (density)
    --, expectation
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

import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified System.Random.MWC.Probability as P
import qualified Control.Monad.ST as ST


--- Test Bed ---


--- Probability Theory ---

type Random s a = P.Prob (ST.ST s) a

type Sample m = S.Vector (SampleDimension m) Double

realize :: Random s a -> IO a
realize = P.withSystemRandom . P.sample

-- | A 'Statistical' 'Manifold' is a 'Manifold' of probability distributions,
-- which all have in common a particular 'SampleSpace'.
class (Manifold m, KnownNat (SampleDimension m)) => Statistical m where
    type SampleDimension m :: Nat

--class (KnownNat (Cardinality m), Statistical m) => Discrete m where
--    type Cardinality m :: Nat
--    sampleSpace :: Proxy m -> S.Vector (Cardinality m) (Sample m)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by MWC Monad.
class Statistical m => Generative c m where
    sample :: Point c m -> Random r (Sample m)

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Statistical m => AbsolutelyContinuous c m where
    density :: Point c m -> Sample m -> Double

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
--expectation
--    :: forall c m. (AbsolutelyContinuous c m, Discrete m)
--    => Point c m
--    -> (Sample m -> Double)
--    -> Double
--expectation p f =
--    let xs = sampleSpace (Proxy :: Proxy m)
--     in S.sum $ S.map (\x -> f x * density p x) xs

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical m => MaximumLikelihood c m where
    mle :: (KnownNat k, 1 <= k) => S.Vector k (Sample m) -> Point c m


--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'dimension' m
-- samples from the given distribution.
initialize :: (Manifold m, Generative d n, SampleDimension n ~ 1) => d # n -> Random r (Point c m)
initialize q = do
    c0s <- S.replicateM $ S.head <$> sample q
    return $ Point c0s

instance (KnownNat k, Statistical m) => Statistical (Replicated k m) where
    type SampleDimension (Replicated k m) = k * SampleDimension m

instance (KnownNat k, Statistical m, Generative c m) => Generative c (Replicated k m) where
    {-# INLINE sample #-}
    sample = fmap S.concat . S.mapM sample . splitReplicated

instance (KnownNat k, Statistical m, AbsolutelyContinuous c m) => AbsolutelyContinuous c (Replicated k m) where
    {-# INLINE density #-}
    density ps xs = S.product $ S.zipWith density (splitReplicated ps) $ S.breakEvery xs

instance (Statistical m, Statistical n) => Statistical (Sum m n) where
    type SampleDimension (Sum m n) = SampleDimension m + SampleDimension n

instance (Generative c m, Generative c n) => Generative c (Sum m n) where
    {-# INLINE sample #-}
    sample pmn = do
        let (pm,pn) = splitSum pmn
        xm <- sample pm
        xn <- sample pn
        return $ xm S.++ xn

instance (AbsolutelyContinuous c m, AbsolutelyContinuous c n) => AbsolutelyContinuous c (Sum m n) where
    {-# INLINE density #-}
    density pmn xmn =
        let (xm,xn) = S.splitAt xmn
            (pm,pn) = splitSum pmn
         in density pm xm * density pn xn


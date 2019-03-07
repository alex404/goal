{-# LANGUAGE
    UndecidableInstances,
    TypeFamilies,
    DataKinds,
    FlexibleContexts,
    MultiParamTypeClasses,
    RankNTypes,
    TypeOperators,
    ScopedTypeVariables,
    FlexibleInstances #-}
-- | Here we provide the basic types and classes for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical
    ( -- * Random
      Random
    , realize
    -- * Statistical Manifolds
    , Statistical (SamplePoint)
    , Sample
    , SamplePoints
    -- * Construction
    , initialize
    , uniformInitialize
    -- * Properties of Distributions
    , Generative (samplePoint)
    , sample
    , AbsolutelyContinuous (density)
    , MaximumLikelihood (mle)
    , Discrete (Cardinality,sampleSpace)
    , pointSampleSpace
    , expectation
    -- * Model Selection
    , akaikesInformationCriterion
    , bayesianInformationCriterion
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


--- Probability Theory ---


-- | A random variable.
type Random s = P.Prob (ST.ST s)

-- | Turn a random variable into an IO action.
realize :: Random s a -> IO a
{-# INLINE realize #-}
realize = P.withSystemRandom . P.sample

-- | A 'Statistical' 'Manifold' is a 'Manifold' of probability distributions,
-- which all have in common a particular 'SampleSpace'.
class Manifold x => Statistical x where
    type SamplePoint x :: *

-- | A 'Vector' of 'SamplePoint's.
type Sample x = [SamplePoint x]

-- | 'SamplePoint' mapped over an 'HList'.
type family SamplePoints (xs :: [*]) where
    SamplePoints '[] = '[]
    SamplePoints (x : xs) = SamplePoint x : SamplePoints xs

-- | Probability distributions for which the sample space is countable. This
-- affords brute force computation of expectations.
class (KnownNat (Cardinality x), Statistical x) => Discrete x where
    type Cardinality x :: Nat
    sampleSpace :: Proxy x -> Sample x

pointSampleSpace :: forall c x . Discrete x => c # x -> Sample x
pointSampleSpace _ = sampleSpace (Proxy :: Proxy x)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by MWC Monad.
class Statistical x => Generative c x where
    samplePoint :: Point c x -> Random r (SamplePoint x)

-- | Generates a 'Vector' of 'samplePoint's.
sample :: Generative c x => Int -> Point c x -> Random r (Sample x)
{-# INLINE sample #-}
sample k = replicateM k . samplePoint

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Statistical x => AbsolutelyContinuous c x where
    density :: Point c x -> SamplePoint x -> Double

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
expectation
    :: forall c x . (AbsolutelyContinuous c x, Discrete x)
    => Point c x
    -> (SamplePoint x -> Double)
    -> Double
{-# INLINE expectation #-}
expectation p f =
     sum $ (\x -> f x * density p x) <$> sampleSpace (Proxy :: Proxy x)

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical x => MaximumLikelihood c x where
    mle :: Sample x -> Point c x


--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'Dimension' m
-- samples from the given distribution.
initialize
    :: (Manifold x, Generative d y, SamplePoint y ~ Double)
    => d # y
    -> Random r (c # x)
initialize q = fromBoxed <$> B.replicateM (samplePoint q)

-- | Generates an initial point on the 'Manifold' m by generating uniform samples from the given vector of bounds.
uniformInitialize :: Manifold x => B.Vector (Dimension x) (Double,Double) -> Random r (Point c x)
uniformInitialize bnds =
    Point . G.convert <$> mapM P.uniformR bnds


--- Model Evaluation ---


-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion
    :: forall c x . AbsolutelyContinuous c x
    => c # x
    -> Sample x
    -> Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
     in 2 * fromIntegral d - 2 * sum (log . density p <$> xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion
    :: forall c x . AbsolutelyContinuous c x
    => c # x
    -> Sample x
    -> Double
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum (log . density p <$> xs)


--- Instances ---


instance (KnownNat k, Statistical x) => Statistical (Replicated k x) where
    type SamplePoint (Replicated k x) = B.Vector k (SamplePoint x)

instance (KnownNat k, Statistical x, Generative c x) => Generative c (Replicated k x) where
    {-# INLINE samplePoint #-}
    samplePoint = B.mapM samplePoint . G.convert . splitReplicated


instance (KnownNat k, Statistical x, AbsolutelyContinuous c x) => AbsolutelyContinuous c (Replicated k x) where
    {-# INLINE density #-}
    density ps xs = B.product $ B.zipWith density (G.convert $ splitReplicated ps) xs

instance Manifold (Sum xs) => Statistical (Sum xs) where
    type SamplePoint (Sum xs) = HList (SamplePoints xs)

instance Generative c (Sum '[]) where
    {-# INLINE samplePoint #-}
    samplePoint _ = return Null

instance (Generative c x, Generative c (Sum xs)) => Generative c (Sum (x : xs)) where
    {-# INLINE samplePoint #-}
    samplePoint pms = do
        let (pm,pms') = splitSum pms
        xm <- samplePoint pm
        xms <- samplePoint pms'
        return $ xm :+: xms

instance AbsolutelyContinuous c (Sum '[]) where
    {-# INLINE density #-}
    density _ _ = 1

instance (AbsolutelyContinuous c x, AbsolutelyContinuous c (Sum xs))
  => AbsolutelyContinuous c (Sum (x : xs)) where
    {-# INLINE density #-}
    density pms (xm :+: xms) =
        let (pm,pms') = splitSum pms
         in density pm xm * density pms' xms

instance (Statistical x, Statistical y) => Statistical (x,y) where
    type SamplePoint (x,y) = (SamplePoint x, SamplePoint y)

instance (Generative c x, Generative c y) => Generative c (x,y) where
    {-# INLINE samplePoint #-}
    samplePoint pmn = do
        let (pm,pn) = splitPair pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c x, AbsolutelyContinuous c y) => AbsolutelyContinuous c (x,y) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitPair pmn
         in density pm xm * density pn xn


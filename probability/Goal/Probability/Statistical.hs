-- | Here we provide the basic types and classes for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical
    ( -- * Random
      Random
    , realize
    -- * Construction
    , initialize
    , uniformInitialize
    -- * Properties of Distributions
    , Generative (sample,samplePoint)
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
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

-- Qualified --

import qualified System.Random.MWC.Probability as P
import qualified Control.Monad.ST as ST

import Foreign.Storable

--- Probability Theory ---


-- | A random variable.
type Random s = P.Prob (ST.ST s)

-- | Turn a random variable into an IO action.
realize :: Random s a -> IO a
{-# INLINE realize #-}
realize = P.withSystemRandom . P.sample

-- | Probability distributions for which the sample space is countable. This
-- affords brute force computation of expectations.
class KnownNat (Cardinality x) => Discrete x s where
    type Cardinality x :: Nat
    sampleSpace :: Proxy x -> [s]

pointSampleSpace :: forall c x s . Discrete x s => c # x -> [s]
pointSampleSpace _ = sampleSpace (Proxy :: Proxy x)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by MWC Monad.
class Manifold x => Generative c x s where
    samplePoint :: Point c x -> Random r s
    samplePoint = fmap head . sample 1
    sample :: Int -> Point c x -> Random r [s]
    sample n = replicateM n . samplePoint

---- | Generates a 'Vector' of 'sample's.
--sample :: Generative c x => Int -> Point c x -> Random r (Sample x)
--{-# INLINE sample #-}
--sample k = replicateM k . sample

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Manifold x => AbsolutelyContinuous c x s where
    density :: Point c x -> s -> Double
    density p = head . densities p . (:[])
    densities :: Point c x -> [s] -> [Double]
    densities p = map (density p)

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
expectation
    :: forall c x s . (AbsolutelyContinuous c x s, Discrete x s)
    => Point c x
    -> (s -> Double)
    -> Double
{-# INLINE expectation #-}
expectation p f =
    let xs = sampleSpace (Proxy :: Proxy x)
     in sum $ zipWith (*) (f <$> xs) (densities p xs)

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class MaximumLikelihood c x s where
    mle :: [s] -> Point c x


--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'Dimension' m
-- samples from the given distribution.
initialize
    :: (Manifold x, Generative d y Double)
    => d # y
    -> Random r (c # x)
initialize q = Point <$> S.replicateM (samplePoint q)

-- | Generates an initial point on the 'Manifold' m by generating uniform samples from the given vector of bounds.
uniformInitialize :: Manifold x => B.Vector (Dimension x) (Double,Double) -> Random r (Point c x)
uniformInitialize bnds =
    Point . G.convert <$> mapM P.uniformR bnds


--- Model Evaluation ---


-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion
    :: forall c x s . (Manifold x, AbsolutelyContinuous c x s)
    => c # x
    -> [s]
    -> Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
     in 2 * fromIntegral d - 2 * sum (log <$> densities p xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion
    :: forall c x s . (AbsolutelyContinuous c x s, Manifold x)
    => c # x
    -> [s]
    -> Double
bayesianInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension x))
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum (log <$> densities p xs)


--- Instances ---


instance (KnownNat k, Storable s, Generative c x s)
  => Generative c (Replicated k x) (S.Vector k s) where
    {-# INLINE samplePoint #-}
    samplePoint = S.mapM samplePoint . splitReplicated


instance (KnownNat k, Storable s, AbsolutelyContinuous c x s)
  => AbsolutelyContinuous c (Replicated k x) (S.Vector k s) where
    {-# INLINE density #-}
    density cxs = S.product . S.zipWith density (splitReplicated cxs)

instance Generative c (Sum '[]) (HList '[]) where
    {-# INLINE samplePoint #-}
    samplePoint _ = return Null

instance (Generative c x s, Generative c (Sum xs) (HList ss))
  => Generative c (Sum (x : xs)) (HList (s : ss)) where
    {-# INLINE samplePoint #-}
    samplePoint pms = do
        let (pm,pms') = splitSum pms
        xm <- samplePoint pm
        xms <- samplePoint pms'
        return $ xm :+: xms

instance AbsolutelyContinuous c (Sum '[]) (HList '[]) where
    {-# INLINE density #-}
    density _ _ = 1

instance (AbsolutelyContinuous c x s, AbsolutelyContinuous c (Sum xs) (HList ss))
  => AbsolutelyContinuous c (Sum (x : xs)) (HList (s : ss)) where
    {-# INLINE density #-}
    density pms (xm :+: xms) =
        let (pm,pms') = splitSum pms
         in density pm xm * density pms' xms

instance (Generative c x s, Generative c y t)
  => Generative c (x,y) (s,t) where
    {-# INLINE samplePoint #-}
    samplePoint pmn = do
        let (pm,pn) = splitPair pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c x s, AbsolutelyContinuous c y t)
  => AbsolutelyContinuous c (x,y) (s,t) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitPair pmn
         in density pm xm * density pn xn


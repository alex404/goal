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
    , type (~>)
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
    -- * Sample Estimates
    , estimateMeanVariance
    , estimateFanoFactor
    , estimateCoefficientOfVariation
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
import qualified Statistics.Sample as S
import qualified Data.Vector.Storable as VS


--- Probability Theory ---


-- | A random variable.
type Random s a = P.Prob (ST.ST s) a

type (~>) s a = Random s a
infixl 1 ~>

-- | Turn a random variable into an IO action.
realize :: Random s a -> IO a
{-# INLINE realize #-}
realize = P.withSystemRandom . P.sample

-- | A 'Statistical' 'Manifold' is a 'Manifold' of probability distributions,
-- which all have in common a particular 'SampleSpace'.
class Manifold m => Statistical m where
    type SamplePoint m :: *

-- | A 'Vector' of 'SamplePoint's.
type Sample m = [SamplePoint m]

-- | 'SamplePoint' mapped over an 'HList'.
type family SamplePoints (ms :: [*]) where
    SamplePoints '[] = '[]
    SamplePoints (m : ms) = SamplePoint m : SamplePoints ms

-- | Probability distributions for which the sample space is countable. This
-- affords brute force computation of expectations.
class (KnownNat (Cardinality m), Statistical m) => Discrete m where
    type Cardinality m :: Nat
    sampleSpace :: Proxy m -> Sample m

pointSampleSpace :: forall c m . Discrete m => c # m -> Sample m
pointSampleSpace _ = sampleSpace (Proxy :: Proxy m)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by MWC Monad.
class Statistical m => Generative c m where
    samplePoint :: Point c m -> Random r (SamplePoint m)

-- | Generates a 'Vector' of 'samplePoint's.
sample :: Generative c m => Int -> Point c m -> Random r (Sample m)
{-# INLINE sample #-}
sample k = replicateM k . samplePoint

-- | If a distribution is 'AbsolutelyContinuous' with respect to a reference
-- measure on its 'SampleSpace', then we may define the 'density' of a
-- probability distribution as the Radon-Nikodym derivative of the probability
-- measure with respect to the base measure.
class Statistical m => AbsolutelyContinuous c m where
    density :: Point c m -> SamplePoint m -> Double

-- | 'expectation' computes the brute force expected value of a 'Finite' set given an appropriate 'density'.
expectation
    :: forall c m . (AbsolutelyContinuous c m, Discrete m)
    => Point c m
    -> (SamplePoint m -> Double)
    -> Double
expectation p f =
     sum $ (\x -> f x * density p x) <$> sampleSpace (Proxy :: Proxy m)

-- | 'mle' computes the 'MaximumLikelihood' estimator.
class Statistical m => MaximumLikelihood c m where
    mle :: Sample m -> Point c m

-- | Estimate the mean and variance of a sample (with Bessel's correction)
estimateMeanVariance
    :: Traversable f
    => f Double
    -> (Double,Double)
{-# INLINE estimateMeanVariance #-}
estimateMeanVariance xs = S.meanVarianceUnb . VS.fromList $ toList xs

---- | Estimate the mean and variance of a sample (with Bessel's correction)
--estimateMeanVariance
--    :: (Traversable f, Real x)
--    => f x
--    -> (Double,Double)
--{-# INLINE estimateMeanVariance #-}
--estimateMeanVariance xs0 =
--    let xs = realToFrac <$> xs0
--        xht = average xs
--        x2s = square . subtract xht <$> xs
--     in (xht, uncurry (/) $ foldr (\e (s,c) -> (e+s,c+1)) (0,-1) x2s)

-- | Estimate the Fano Factor of a sample.
estimateFanoFactor
    :: Traversable f
    => f Double
    -> Double
{-# INLINE estimateFanoFactor #-}
estimateFanoFactor xs =
    let (mu,vr) = estimateMeanVariance xs
     in vr / mu

estimateCoefficientOfVariation :: Traversable f => f Double -> Double
{-# INLINE estimateCoefficientOfVariation #-}
estimateCoefficientOfVariation zs =
    let (mu,vr) = estimateMeanVariance zs
     in sqrt vr / mu

--- Construction ---


-- | Generates an initial point on the 'Manifold' m by generating 'Dimension' m
-- samples from the given distribution.
initialize
    :: (Manifold m, Generative d n, SamplePoint n ~ Double)
    => d # n
    -> Random r (Point c m)
initialize q = fromBoxed <$> B.replicateM (samplePoint q)

-- | Generates an initial point on the 'Manifold' m by generating uniform samples from the given vector of bounds.
uniformInitialize :: Manifold m => B.Vector (Dimension m) (Double,Double) -> Random r (Point c m)
uniformInitialize bnds =
    Point . G.convert <$> mapM P.uniformR bnds


--- Model Evaluation ---


-- | Calculate the AIC for a given model and sample.
akaikesInformationCriterion
    :: forall c m . AbsolutelyContinuous c m
    => c # m
    -> Sample m
    -> Double
akaikesInformationCriterion p xs =
    let d = natVal (Proxy :: Proxy (Dimension m))
     in 2 * fromIntegral d - 2 * sum (log . density p <$> xs)

-- | Calculate the BIC for a given model and sample.
bayesianInformationCriterion
    :: forall c m . AbsolutelyContinuous c m
    => c # m
    -> Sample m
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

instance Manifold (Sum ms) => Statistical (Sum ms) where
    type SamplePoint (Sum ms) = HList (SamplePoints ms)

instance Generative c (Sum '[]) where
    {-# INLINE samplePoint #-}
    samplePoint _ = return Null

instance (Generative c m, Generative c (Sum ms)) => Generative c (Sum (m : ms)) where
    {-# INLINE samplePoint #-}
    samplePoint pms = do
        let (pm,pms') = splitSum pms
        xm <- samplePoint pm
        xms <- samplePoint pms'
        return $ xm :+: xms

instance AbsolutelyContinuous c (Sum '[]) where
    {-# INLINE density #-}
    density _ _ = 1

instance (AbsolutelyContinuous c m, AbsolutelyContinuous c (Sum ms))
  => AbsolutelyContinuous c (Sum (m : ms)) where
    {-# INLINE density #-}
    density pms (xm :+: xms) =
        let (pm,pms') = splitSum pms
         in density pm xm * density pms' xms

instance (Statistical m, Statistical n) => Statistical (m,n) where
    type SamplePoint (m,n) = (SamplePoint m, SamplePoint n)

instance (Generative c m, Generative c n) => Generative c (m,n) where
    {-# INLINE samplePoint #-}
    samplePoint pmn = do
        let (pm,pn) = splitPair pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c m, AbsolutelyContinuous c n) => AbsolutelyContinuous c (m,n) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitPair pmn
         in density pm xm * density pn xn


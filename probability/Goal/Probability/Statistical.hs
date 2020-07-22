{-# LANGUAGE UndecidableInstances #-}

-- | Core types, classes, and functions for working with manifolds of
-- probability distributions.
module Goal.Probability.Statistical
    ( -- * Random
      Random
    , Statistical (SamplePoint)
    , Sample
    , SamplePoints
    , realize
    , observableSample
    -- * Initializiation
    , initialize
    , uniformInitialize
    , uniformInitialize'
    -- * Properties of Distributions
    , Generative (sample,samplePoint)
    , AbsolutelyContinuous (density,densities)
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
{-# INLINE realize #-}
realize = P.withSystemRandom . P.sample

-- | Probability distributions for which the sample space is countable. This
-- affords brute force computation of expectations.
class KnownNat (Cardinality x) => Discrete x where
    type Cardinality x :: Nat
    sampleSpace :: Proxy x -> Sample x

-- | Convenience function for getting the sample space of a 'Discrete'
-- probability distribution.
pointSampleSpace :: forall c x . Discrete x => c # x -> Sample x
{-# INLINE pointSampleSpace #-}
pointSampleSpace _ = sampleSpace (Proxy :: Proxy x)

-- | A distribution is 'Generative' if we can 'sample' from it. Generation is
-- powered by @mwc-random@.
class Statistical x => Generative c x where
    samplePoint :: Point c x -> Random r (SamplePoint x)
    {-# INLINE samplePoint #-}
    samplePoint = fmap head . sample 1
    sample :: Int -> Point c x -> Random r (Sample x)
    {-# INLINE sample #-}
    sample n = replicateM n . samplePoint

-- | A 'SamplePoint' construction for 'HList's.
type family SamplePoints (xs :: [Type]) where
    SamplePoints '[] = '[]
    SamplePoints (x : xs) = SamplePoint x : SamplePoints xs


type family HHead as where
    HHead (HList (a ': as)) = '[a]

type family Head as where
    Head (HList (a ': as)) = a

type Observation x = Head (SamplePoint x)

type Observations x = [Observation x]


observableSample
    :: (Generative c x, SamplePoint x ~ HList (a : as))
    => Int -> c # x -> Random r (Observations x)
{-# INLINE observableSample #-}
observableSample nsmps p = map hHead <$> sample nsmps p


-- | The distributions \(P \in \mathcal M\) in a 'Statistical' 'Manifold'
-- \(\mathcal M\) are 'AbsolutelyContinuous' if there is a reference measure
-- \(\mu\) and a function \(p\) such that
-- \(P(A) = \int_A p d\mu\). We refer to \(p(x)\) as the 'density' of the
-- probability distribution.
class Statistical x => AbsolutelyContinuous c x where
    density :: Point c x -> SamplePoint x -> Double
    {-# INLINE density #-}
    density p = head . densities p . (:[])
    densities :: Point c x -> Sample x -> [Double]
    {-# INLINE densities #-}
    densities p = map (density p)

-- | 'expectation' computes the brute force expected value of a 'Finite' set
-- given an appropriate 'density'.
expectation
    :: forall c x . (AbsolutelyContinuous c x, Discrete x)
    => Point c x
    -> (SamplePoint x -> Double)
    -> Double
{-# INLINE expectation #-}
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
{-# INLINE initialize #-}
initialize q = Point <$> S.replicateM (samplePoint q)

-- | Generates an initial point on the target 'Manifold' by generating uniform
-- samples from the given vector of bounds.
uniformInitialize' :: Manifold x => B.Vector (Dimension x) (Double,Double) -> Random r (Point c x)
{-# INLINE uniformInitialize' #-}
uniformInitialize' bnds =
    Point . G.convert <$> mapM P.uniformR bnds

-- | Generates an initial point on the target 'Manifold' by generating uniform
-- samples from the given vector of bounds.
uniformInitialize :: Manifold x => (Double,Double) -> Random r (Point c x)
{-# INLINE uniformInitialize #-}
uniformInitialize bnds =
    Point <$> S.replicateM (P.uniformR bnds)



--- Instances ---


-- Replicated --

instance (Statistical x, KnownNat k, Storable (SamplePoint x))
  => Statistical (Replicated k x) where
    type SamplePoint (Replicated k x) = S.Vector k (SamplePoint x)

instance (KnownNat k, Generative c x, Storable (SamplePoint x))
  => Generative c (Replicated k x) where
    {-# INLINE samplePoint #-}
    samplePoint = S.mapM samplePoint . splitReplicated

instance (KnownNat k, Storable (SamplePoint x), AbsolutelyContinuous c x)
  => AbsolutelyContinuous c (Replicated k x) where
    {-# INLINE density #-}
    density cxs = S.product . S.zipWith density (splitReplicated cxs)

instance (KnownNat k, LogLikelihood c x s, Storable s)
  => LogLikelihood c (Replicated k x) (S.Vector k s) where
    {-# INLINE logLikelihood #-}
    logLikelihood cxs ps = S.sum . S.imap subLogLikelihood $ splitReplicated ps
        where subLogLikelihood fn = logLikelihood (flip S.index fn <$> cxs)
    {-# INLINE logLikelihoodDifferential #-}
    logLikelihoodDifferential cxs ps =
        joinReplicated . S.imap subLogLikelihoodDifferential $ splitReplicated ps
            where subLogLikelihoodDifferential fn = logLikelihoodDifferential (flip S.index fn <$> cxs)

-- Sum --


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

-- Pair --

instance (Statistical x, Statistical y)
  => Statistical (x,y) where
    type SamplePoint (x,y) = (SamplePoint x, SamplePoint y)


instance (Generative c x, Generative c y) => Generative c (x,y) where
    {-# INLINE samplePoint #-}
    samplePoint pmn = do
        let (pm,pn) = splitPair pmn
        xm <- samplePoint pm
        xn <- samplePoint pn
        return (xm,xn)

instance (AbsolutelyContinuous c x, AbsolutelyContinuous c y)
  => AbsolutelyContinuous c (x,y) where
    {-# INLINE density #-}
    density pmn (xm,xn) =
        let (pm,pn) = splitPair pmn
         in density pm xm * density pn xn


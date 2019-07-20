{-# LANGUAGE UndecidableInstances,TypeApplications #-}
-- | Definitions for working with exponential families.
module Goal.Probability.ExponentialFamily
    ( -- * Exponential Families
    ExponentialFamily (sufficientStatistic, baseMeasure)
    , LegendreExponentialFamily
    , DuallyFlatExponentialFamily
    , sufficientStatisticT
    , exponentialFamilyDensity
    , unnormalizedDensity
    , unnormalizedLogDensity
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
    , crossEntropyDifferential
    , stochasticCrossEntropyDifferential
    , stochasticInformationProjectionDifferential
    -- ** Conditional Distributions
    , (>.>*)
    , (>$>*)
    , (*<.<)
    , (*<$<)
    -- *** Maximum Likelihood
    , exponentialFamilyLogLikelihood
    , exponentialFamilyLogLikelihoodDifferential
    , conditionalLogLikelihood
    , conditionalLogLikelihoodDifferential
    ) where

--- Imports ---


-- Package --

import Goal.Probability.Statistical

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Foreign.Storable

--- Exponential Families ---

-- Source Chart --

-- | A parameterization which represents the standard or typical parameterization of
-- the given manifold, e.g. the 'Poisson' rate or 'Normal' mean and standard deviation.
data Source


-- | A parameterization in terms of the natural coordinates of an exponential family.
data Natural

-- | A representation in terms of the mean sufficient statistics of an exponential family.
data Mean

instance Primal Natural where
    type Dual Natural = Mean

instance Primal Mean where
    type Dual Mean = Natural

-- | Expresses an exponential family distribution in natural coordinates.
toNatural :: (Transition c Natural x) => c # x -> Natural # x
{-# INLINE toNatural #-}
toNatural = transition

-- | Expresses an exponential family distribution in mean coordinates.
toMean :: (Transition c Mean x) => c # x -> Mean # x
{-# INLINE toMean #-}
toMean = transition

-- | Expresses an exponential family distribution in source coordinates.
toSource :: (Transition c Source x) => c # x -> Source # x
{-# INLINE toSource #-}
toSource = transition


-- | A 'Statistical' 'Manifold' is a member of the 'ExponentialFamily' if we can
-- specify a 'sufficientStatistic' of fixed length, and a 'baseMeasure' which
-- fully determines the normalization of the distribution in exponential family
-- form.
--
-- 'ExponentialFamily' distributions theoretically have a 'Riemannian' geometry
-- given by the Fisher information metric, given rise to the 'Dual' coordinate
-- system given by the 'Natural' and 'Mean' coordinates. However, not all
-- distributions (e.g. the von Mises distribution) allow for closed form
-- expressions of the relevant structures, and so we define a distinct class for
-- this purpose.
class Statistical x => ExponentialFamily x where
    sufficientStatistic :: SamplePoint x -> Mean # x
    baseMeasure :: Proxy x -> SamplePoint x -> Double

-- | When the 'Riemannian' properties of the given 'ExponentialFamily' may be
-- computed in closed-form, then we refer to it as a
-- 'ClosedFormExponentialFamily'.
type LegendreExponentialFamily x =
    ( PotentialCoordinates x ~ Natural, Legendre x, ExponentialFamily x
    , Transition (PotentialCoordinates x) (Dual (PotentialCoordinates x)) x )

type DuallyFlatExponentialFamily x =
    ( LegendreExponentialFamily x, DuallyFlat x
    , Transition (Dual (PotentialCoordinates x)) (PotentialCoordinates x) x )


-- | The sufficient statistic of a traversable set of iid random variables.
sufficientStatisticT
    :: (ExponentialFamily x, Traversable f)
    => f (SamplePoint x) -> Mean # x
{-# INLINE sufficientStatisticT #-}
sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs

-- | A function for computing the relative entropy, also known as the KL-divergence.
relativeEntropy :: (DuallyFlatExponentialFamily x) => Mean # x -> Natural # x -> Double
{-# INLINE relativeEntropy #-}
relativeEntropy = flip canonicalDivergence

-- | A function for computing the cross-entropy, which is the relative entropy plus the entropy of the first distribution.
crossEntropy :: DuallyFlatExponentialFamily x => Mean # x -> Natural # x -> Double
{-# INLINE crossEntropy #-}
crossEntropy mp nq = potential nq - (mp <.> nq)

-- | The differential of the cross-entropy with respect to the parameters of the second argument.
crossEntropyDifferential :: (LegendreExponentialFamily x) => Mean # x -> Natural # x -> Mean # x
{-# INLINE crossEntropyDifferential #-}
crossEntropyDifferential mp nq = transition nq <-> mp

-- | The differential of the cross-entropy with respect to the parameters of the
-- second argument, based only on samples from the two distributions.
stochasticCrossEntropyDifferential
    :: ExponentialFamily x
    => Sample x -- ^ True Samples
    -> Sample x -- ^ Model Samples
    -> Mean # x -- ^ Differential Estimate
{-# INLINE stochasticCrossEntropyDifferential #-}
stochasticCrossEntropyDifferential pxs qxs =
    sufficientStatisticT qxs <-> sufficientStatisticT pxs

-- | The differential of the dual relative entropy.
stochasticInformationProjectionDifferential
    :: ExponentialFamily x
    => Natural # x -- ^ Model Distribution
    -> Sample x -- ^ Model Samples
    -> (SamplePoint x -> Double) -- ^ Unnormalized log-density of target distribution
    -> Mean # x -- ^ Differential Estimate
{-# INLINE stochasticInformationProjectionDifferential #-}
stochasticInformationProjectionDifferential px xs f =
    let mxs = sufficientStatistic <$> xs
        mys = (\x -> sufficientStatistic x <.> px - f x) <$> xs
        ln = fromIntegral $ length xs
        mxht = ln /> foldr1 (<+>) mxs
        myht = sum mys / ln
        cvr = (ln - 1) /> foldr1 (<+>) [ (my - myht) .> (mx <-> mxht) | (mx,my) <- zip mxs mys ]
     in cvr

-- | The exact density of an exponential family distribution with an exact
-- expression for the log-partition function.
exponentialFamilyDensity
    :: (ExponentialFamily x, LegendreExponentialFamily x)
    => Natural # x -> SamplePoint x -> Double
{-# INLINE exponentialFamilyDensity #-}
exponentialFamilyDensity p x = unnormalizedDensity p x * (exp . negate $ potential p)

-- | The unnormalized density of an arbitrary exponential family distribution.
unnormalizedDensity :: forall x . ExponentialFamily x => Natural # x -> SamplePoint x -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy @ x) x

-- | The unnormalized log-density of an arbitrary exponential family distribution.
unnormalizedLogDensity :: forall x . ExponentialFamily x => Natural # x -> SamplePoint x -> Double
unnormalizedLogDensity p x =
    p <.> sufficientStatistic x  + log (baseMeasure (Proxy @ x) x)


-- | An approximate cross-entropy based on samples from the first argument, and
-- an exact expression for the second argument.
exponentialFamilyLogLikelihood
    :: forall x . (ExponentialFamily x, LegendreExponentialFamily x)
    => Sample x -> Natural # x -> Double
{-# INLINE exponentialFamilyLogLikelihood #-}
exponentialFamilyLogLikelihood xs nq =
    let mp = sufficientStatisticT xs
        bm = average $ log . baseMeasure (Proxy :: Proxy x) <$> xs
     in -potential nq + (mp <.> nq) + bm

-- | An approximate cross-entropy differential based on samples from the first argument, and
-- an exact expression for differentiated distribution.
exponentialFamilyLogLikelihoodDifferential
    :: (ExponentialFamily x, LegendreExponentialFamily x)
    => Sample x -> Natural # x -> Mean # x
{-# INLINE exponentialFamilyLogLikelihoodDifferential #-}
exponentialFamilyLogLikelihoodDifferential xs nq =
    let mp = sufficientStatisticT xs
     in mp <-> transition nq


-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
conditionalLogLikelihood
    :: ( Map Mean Natural f y x, ExponentialFamily x, LogLikelihood Natural y t)
    => [(t,SamplePoint x)] -- ^ (Output,Input) sample
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE conditionalLogLikelihood #-}
conditionalLogLikelihood yxs f =
    let (ys,xs) = unzip yxs
     in average . zipWith logLikelihood ((:[]) <$> ys) $ f >$>* xs

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates (this is
-- primarily of internal use).
conditionalLogLikelihoodDifferential
    :: ( Propagate Mean Natural f y x, LogLikelihood Natural y t
       , LegendreExponentialFamily y, ExponentialFamily x )
    => [(t,SamplePoint x)] -- ^ Output/Input mean distributions
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
{-# INLINE conditionalLogLikelihoodDifferential #-}
conditionalLogLikelihoodDifferential yxs f =
    let (ys,xs) = unzip yxs
        (df,yhts) = propagate mys (sufficientStatistic <$> xs) f
        mys = zipWith logLikelihoodDifferential ((:[]) <$> ys) yhts
     in df


--- Conditional Distributions ---


-- | Applies the given conditional distribution to a samplePoint.
(>.>*) :: (Map Mean c f y x, ExponentialFamily x)
       => Function Mean c # f y x
       -> SamplePoint x
       -> c # y
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application on samples.
(>$>*) :: (Map Mean c f y x, ExponentialFamily x)
       => Function Mean c # f y x
       -> Sample x
       -> [c # y]
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the transpose conditional distribution to a samplePoint.
(*<.<) :: (Map Mean Natural f x y, Bilinear f y x, ExponentialFamily y)
       => SamplePoint y
       -> Natural #> f y x
       -> Natural # x
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped application on samples.
(*<$<) :: (Map Mean Natural f x y, Bilinear f y x, ExponentialFamily y)
       => Sample y
       -> Natural #> f y x
       -> [Natural # x]
{-# INLINE (*<$<) #-}
(*<$<) xs p = (sufficientStatistic <$> xs) <$< p

infix 8 *<.<
infix 8 *<$<


--- Model Selection ---



--- Internal ---


replicatedBaseMeasure0 :: (ExponentialFamily x, Storable (SamplePoint x), KnownNat k)
                       => Proxy x -> Proxy (Replicated k x) -> S.Vector k (SamplePoint x) -> Double
{-# INLINE replicatedBaseMeasure0  #-}
replicatedBaseMeasure0 prxym _ xs = S.product $ S.map (baseMeasure prxym) xs

sumBaseMeasure
    :: (ExponentialFamily x, ExponentialFamily (Sum xs))
    => Proxy x
    -> Proxy (Sum xs)
    -> Proxy (Sum (x : xs))
    -> SamplePoint (Sum (x : xs))
    -> Double
{-# INLINE sumBaseMeasure #-}
sumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

pairBaseMeasure
    :: (ExponentialFamily x, ExponentialFamily y)
    => Proxy x
    -> Proxy y
    -> Proxy (x,y)
    -> SamplePoint (x,y)
    -> Double
{-# INLINE pairBaseMeasure #-}
pairBaseMeasure prxym prxyn _ (xm,xn) =
     baseMeasure prxym xm * baseMeasure prxyn xn


--- Instances ---


-- Replicated --

instance Transition Natural Natural x where
    {-# INLINE transition #-}
    transition = id

instance Transition Mean Mean x where
    {-# INLINE transition #-}
    transition = id

instance Transition Source Source x where
    {-# INLINE transition #-}
    transition = id

instance (ExponentialFamily x, Storable (SamplePoint x), KnownNat k)
  => ExponentialFamily (Replicated k x) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = joinReplicated $ S.map sufficientStatistic xs
    {-# INLINE baseMeasure #-}
    baseMeasure = replicatedBaseMeasure0 Proxy

-- Sum --

instance ExponentialFamily (Sum '[]) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic _ = zero
    {-# INLINE baseMeasure #-}
    baseMeasure _ _ = 1

instance (ExponentialFamily x, ExponentialFamily (Sum xs)) => ExponentialFamily (Sum (x : xs)) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm :+: xms) =
         joinSum (sufficientStatistic xm) (sufficientStatistic xms)
    {-# INLINE baseMeasure #-}
    baseMeasure = sumBaseMeasure Proxy Proxy

instance (ExponentialFamily x, ExponentialFamily y) => ExponentialFamily (x,y) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm,xn) =
         joinPair (sufficientStatistic xm) (sufficientStatistic xn)
    {-# INLINE baseMeasure #-}
    baseMeasure = pairBaseMeasure Proxy Proxy

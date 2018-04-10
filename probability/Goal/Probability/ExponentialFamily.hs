-- | Definitions for working with exponential families.
module Goal.Probability.ExponentialFamily
    ( -- * Exponential Families
    ExponentialFamily (sufficientStatistic, baseMeasure)
    , ClosedFormExponentialFamily
    , sufficientStatisticT
    , unnormalizedDensity
    , exponentialFamilyDensity
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
    , crossEntropyDifferential
    , stochasticCrossEntropy
    , stochasticCrossEntropyDifferential
    , stochasticConditionalCrossEntropy
    , stochasticConditionalCrossEntropyDifferential
    , estimateStochasticCrossEntropyDifferential
    -- ** Conditional Distributions
    , (>.>*)
    , (>$>*)
    , (*<.<)
    , (*<$<)
    ) where

--- Imports ---


-- Package --

import Goal.Probability.Statistical

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

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

toNatural :: (Transition c Natural m) => Point c m -> Point Natural m
toNatural = transition

toMean :: (Transition c Mean m) => Point c m -> Point Mean m
toMean = transition

toSource :: (Transition c Source m) => Point c m -> Point Source m
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
class Statistical m => ExponentialFamily m where
    sufficientStatistic :: Sample m -> Point Mean m
    baseMeasure :: Proxy m -> Sample m -> Double

-- | When the 'Riemannian' properties of the given 'ExponentialFamily' may be
-- computed in closed-form, then we refer to it as a
-- 'ClosedFormExponentialFamily'.
type ClosedFormExponentialFamily m =
    ( ExponentialFamily m, Legendre Natural m, Legendre Mean m, AbsolutelyContinuous Natural m
    , Transition Natural Mean m, Transition Mean Natural m )

-- | The sufficient statistic of N iid random variables.
sufficientStatisticR
    :: (ExponentialFamily m, KnownNat k, 1 <= k)
    => S.Vector k (Sample m) -> Point Mean m
{-# INLINE sufficientStatisticT #-}
sufficientStatisticR xs = averagePoint (S.map sufficientStatistic xs)

-- | A function for computing the relative entropy, also known as the KL-divergence.
relativeEntropy
    :: (ClosedFormExponentialFamily m, Transition c Mean m, Transition d Natural m)
    => Point c m -> Point d m -> Double
{-# INLINE relativeEntropy #-}
relativeEntropy p q = divergence (toMean p) (toNatural q)

-- | A function for computing the relative entropy, also known as the KL-divergence.
crossEntropy
    :: (ClosedFormExponentialFamily m, Transition c Mean m, Transition d Natural m)
    => Point c m -> Point d m -> Double
{-# INLINE crossEntropy #-}
crossEntropy p q =
    let mp = toMean p
        nq = toNatural q
     in potential nq - (mp <.> nq)

-- | A function for computing the relative entropy, also known as the KL-divergence.
crossEntropyDifferential
    :: (ClosedFormExponentialFamily m, Transition c Mean m, Transition d Natural m)
    => Point c m -> Point d m -> CotangentVector Natural m
{-# INLINE crossEntropyDifferential #-}
crossEntropyDifferential p q =
    let mp = primalIsomorphism $ toMean p
        nq = toNatural q
     in potentialDifferential nq <-> mp

-- | A function for computing the relative entropy, also known as the KL-divergence.
estimateStochasticCrossEntropyDifferential0
    :: (KnownNat k, 1 <= k, ExponentialFamily m)
    => S.Vector k (Sample m) -- ^ True Samples
    -> S.Vector k (Sample m) -- ^ Model Samples
    -> Point Mean m -- ^ Differential Estimate
{-# INLINE estimateStochasticCrossEntropyDifferential0 #-}
estimateStochasticCrossEntropyDifferential0 pxs qxs =
    sufficientStatisticT qxs <-> sufficientStatisticT pxs

-- | A function for computing the relative entropy, also known as the KL-divergence.
estimateStochasticCrossEntropyDifferential1
    :: (ExponentialFamily m)
    => Point Mean m -- ^ Differential Estimate
    -> CotangentVector Natural m -- ^ Differential Estimate
{-# INLINE estimateStochasticCrossEntropyDifferential1 #-}
estimateStochasticCrossEntropyDifferential1 = Point . coordinates

-- | A function for computing the relative entropy, also known as the KL-divergence.
estimateStochasticCrossEntropyDifferential
    :: (KnownNat k, 1 <= k, ExponentialFamily m)
    => S.Vector k (Sample m) -- ^ True Samples
    -> S.Vector k (Sample m) -- ^ Model Samples
    -> CotangentVector Natural m -- ^ Differential Estimate
{-# INLINE estimateStochasticCrossEntropyDifferential #-}
estimateStochasticCrossEntropyDifferential pxs qxs =
    estimateStochasticCrossEntropyDifferential1 $ estimateStochasticCrossEntropyDifferential0 pxs qxs

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticCrossEntropy
    :: (KnownNat k, 1 <= k, ClosedFormExponentialFamily m)
    => S.Vector k (Sample m) -> Point Natural m -> Double
{-# INLINE stochasticCrossEntropy #-}
stochasticCrossEntropy xs nq =
    let mp = sufficientStatisticT xs
     in potential nq - (mp <.> nq)

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticCrossEntropyDifferential
    :: (KnownNat k, 1 <= k, ClosedFormExponentialFamily m)
    => S.Vector k (Sample m) -> Point Natural m -> CotangentVector Natural m
{-# INLINE stochasticCrossEntropyDifferential #-}
stochasticCrossEntropyDifferential xs nq =
    let mp = sufficientStatisticT xs
     in potentialDifferential nq <-> primalIsomorphism mp

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticConditionalCrossEntropy
    :: (KnownNat k, Apply Mean Natural f, ExponentialFamily (Domain f), ClosedFormExponentialFamily (Codomain f))
    => S.Vector k (Sample (Domain f))
    -> S.Vector k (Sample (Codomain f))
    -> Point (Mean ~> Natural) f
    -> Double
{-# INLINE stochasticConditionalCrossEntropy #-}
stochasticConditionalCrossEntropy xs ys f =
    S.average. S.zipWith stochasticCrossEntropy (S.map S.singleton ys) . splitReplicated $ f >$>* xs

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticConditionalCrossEntropyDifferential
    :: ( Propagate Mean Natural f
       , ExponentialFamily (Domain f)
       , ClosedFormExponentialFamily (Codomain f)
       , KnownNat k, 1 <= k)
      => S.Vector k (Sample (Domain f))
      -> S.Vector k (Sample (Codomain f))
      -> Mean ~> Natural # f
      -> CotangentVector (Mean ~> Natural) f
{-# INLINE stochasticConditionalCrossEntropyDifferential #-}
stochasticConditionalCrossEntropyDifferential xs ys f =
    let (df,yhts) = propagate mys (joinReplicated $ S.map sufficientStatistic xs) f
        mys = joinReplicated $ S.zipWith differentiator ys (splitReplicated yhts)
     in primalIsomorphism df
        where differentiator y yht =
                  dualIsomorphism $ stochasticCrossEntropyDifferential (S.singleton y) yht

unnormalizedDensity :: forall m. ExponentialFamily m => Point Natural m -> Sample m -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy :: Proxy m) x

exponentialFamilyDensity
    :: (ExponentialFamily m, Legendre Natural m)
    => Point Natural m -> Sample m -> Double
{-# INLINE exponentialFamilyDensity #-}
exponentialFamilyDensity p x = unnormalizedDensity p x * (exp . negate $ potential p)

replicatedBaseMeasure0 :: (ExponentialFamily m, Storable (Sample m), KnownNat k)
                       => Proxy m -> Proxy (Replicated k m) -> Sample (Replicated k m) -> Double
{-# INLINE replicatedBaseMeasure0  #-}
replicatedBaseMeasure0 prxym _ xs = S.sum . S.map (baseMeasure prxym) $ S.breakEvery xs

sumBaseMeasure0
    :: (ExponentialFamily m, ExponentialFamily n)
    => Proxy m
    -> Proxy n
    -> Proxy (Sum m n)
    -> Sample (Sum m n)
    -> Double
{-# INLINE sumBaseMeasure0  #-}
sumBaseMeasure0 prxym prxyn _ xmn =
    let (xm,xn) = S.splitAt xmn
     in baseMeasure prxym xm * baseMeasure prxyn xn


--- Conditional Distributions ---


-- | Applies the given conditional distribution to a sample from the 'SampleSpace' of the 'Domain'.
(>.>*) :: (Apply Mean c m, ExponentialFamily (Domain m))
       => Point (Mean ~> c) m
       -> Sample (Domain m)
       -> Point c (Codomain m)
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application on samples.
(>$>*) :: (Apply Mean c m, ExponentialFamily (Domain m), KnownNat k)
       => Point (Mean ~> c) m
       -> S.Vector k (Sample (Domain m))
       -> Point c (Replicated k (Codomain m))
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> joinReplicated (S.map sufficientStatistic xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the given conditional distribution to a sample from the 'SampleSpace' of the 'Domain'.
(*<.<) :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f))
       => Sample (Codomain f)
       -> Point (Function Mean Natural) f
       -> Point Natural (Domain f)
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped application on samples.
(*<$<) :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f), KnownNat k)
       => S.Vector k (Sample (Codomain f))
       -> Point (Function Mean Natural) f
       -> Point Natural (Replicated k (Domain f))
{-# INLINE (*<$<) #-}
(*<$<) xs p = joinReplicated (S.map sufficientStatistic xs) <$< p

infix 8 *<.<
infix 8 *<$<



--- Instances ---


-- Transitions --

instance Legendre Mean m => Transition Mean Natural m where
    {-# INLINE transition #-}
    transition = dualTransition

instance Legendre Natural m => Transition Natural Mean m where
    {-# INLINE transition #-}
    transition = dualTransition


-- Replicated --

instance (Legendre Natural m, Riemannian Natural m, KnownNat k) => Riemannian Natural (Replicated k m) where
    --{-# INLINE metric #-}
    --metric = hessian potential
    {-# INLINE flat #-}
    flat = replicatedJoinTangentPair . S.map flat . replicatedSplitTangentPair
    {-# INLINE sharp #-}
    sharp = replicatedJoinTangentPair . S.map sharp . replicatedSplitTangentPair

instance (ExponentialFamily m, ExponentialFamily n) => ExponentialFamily (Sum m n) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xmn =
        let (xm,xn) = S.splitAt xmn
         in joinSum (sufficientStatistic xm) (sufficientStatistic xn)
    {-# INLINE baseMeasure #-}
    baseMeasure = sumBaseMeasure0 Proxy Proxy

instance (ExponentialFamily m, KnownNat k) => ExponentialFamily (Replicated k m) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = joinReplicated $ S.map sufficientStatistic $ S.breakEvery xs
    {-# INLINE baseMeasure #-}
    baseMeasure = replicatedBaseMeasure0 Proxy

instance (Manifold m, Manifold n, Transition Natural Source m, Transition Natural Source n) => Transition Natural Source (Sum m n) where
    {-# INLINE transition #-}
    transition pmn =
        let (pm,pn) = splitSum pmn
         in joinSum (transition pm) (transition pn)

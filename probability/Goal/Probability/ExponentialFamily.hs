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
    -- ** Entropies and their Differentials
    , relativeEntropy
    , crossEntropy
    , crossEntropyDifferential
    , stochasticCrossEntropy
    , stochasticCrossEntropyDifferential
    , estimateStochasticCrossEntropyDifferential
    -- *** Conditional
    , stochasticConditionalCrossEntropy
    , stochasticConditionalCrossEntropyDifferential
    , stochasticConditionalCrossEntropyDifferential0
    -- *** Dual
    --, stochasticDualRelativeEntropyDifferential
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

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Boxed as B

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
{-# INLINE toNatural #-}
toNatural = transition

toMean :: (Transition c Mean m) => Point c m -> Point Mean m
{-# INLINE toMean #-}
toMean = transition

toSource :: (Transition c Source m) => Point c m -> Point Source m
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
class Statistical m => ExponentialFamily m where
    sufficientStatistic :: SamplePoint m -> Point Mean m
    baseMeasure :: Proxy m -> SamplePoint m -> Double

-- | When the 'Riemannian' properties of the given 'ExponentialFamily' may be
-- computed in closed-form, then we refer to it as a
-- 'ClosedFormExponentialFamily'.
type ClosedFormExponentialFamily m =
    ( ExponentialFamily m, Legendre Natural m, Legendre Mean m, AbsolutelyContinuous Natural m
    , Transition Natural Mean m, Transition Mean Natural m )

-- | The sufficient statistic of N iid random variables.
sufficientStatisticT
    :: (ExponentialFamily m, KnownNat k, 1 <= k)
    => Sample k m -> Point Mean m
{-# INLINE sufficientStatisticT #-}
sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs

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
    => Sample k m -- ^ True Samples
    -> Sample k m -- ^ Model Samples
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
    => Sample k m -- ^ True Samples
    -> Sample k m -- ^ Model Samples
    -> CotangentVector Natural m -- ^ Differential Estimate
{-# INLINE estimateStochasticCrossEntropyDifferential #-}
estimateStochasticCrossEntropyDifferential pxs qxs =
    estimateStochasticCrossEntropyDifferential1 $ estimateStochasticCrossEntropyDifferential0 pxs qxs

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticCrossEntropy
    :: (KnownNat k, 1 <= k, ClosedFormExponentialFamily m)
    => Sample k m -> Point Natural m -> Double
{-# INLINE stochasticCrossEntropy #-}
stochasticCrossEntropy xs nq =
    let mp = sufficientStatisticT xs
     in potential nq - (mp <.> nq)

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticCrossEntropyDifferential
    :: (KnownNat k, 1 <= k, ClosedFormExponentialFamily m)
    => Sample k m -> Point Natural m -> CotangentVector Natural m
{-# INLINE stochasticCrossEntropyDifferential #-}
stochasticCrossEntropyDifferential xs nq =
    let mp = sufficientStatisticT xs
     in potentialDifferential nq <-> primalIsomorphism mp

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticConditionalCrossEntropy
    :: (KnownNat k, Apply Mean Natural f, ExponentialFamily (Domain f), ClosedFormExponentialFamily (Codomain f))
    => Sample k (Domain f)
    -> Sample k (Codomain f)
    -> Point (Mean ~> Natural) f
    -> Double
{-# INLINE stochasticConditionalCrossEntropy #-}
stochasticConditionalCrossEntropy xs ys f =
    average . G.zipWith stochasticCrossEntropy (B.singleton <$> ys) . G.convert . splitReplicated $ f >$>* xs

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticConditionalCrossEntropyDifferential0
    :: ( Propagate Mean Natural f
       , ExponentialFamily (Domain f)
       , ClosedFormExponentialFamily (Codomain f)
       , KnownNat k, 1 <= k)
      => Point Mean (Replicated k (Domain f))
      -> Point Mean (Replicated k (Codomain f))
      -> Mean ~> Natural # f
      -> CotangentVector (Mean ~> Natural) f
{-# INLINE stochasticConditionalCrossEntropyDifferential0 #-}
stochasticConditionalCrossEntropyDifferential0 xs ys f =
    let (df,yhts) = propagate mys xs f
        mys = dualIsomorphism $ potentialDifferential yhts <-> primalIsomorphism ys
     in primalIsomorphism df

-- | A function for computing the relative entropy, also known as the KL-divergence.
stochasticConditionalCrossEntropyDifferential
    :: ( Propagate Mean Natural f
       , ExponentialFamily (Domain f)
       , ClosedFormExponentialFamily (Codomain f)
       , KnownNat k, 1 <= k)
      => Sample k (Domain f)
      -> Sample k (Codomain f)
      -> Mean ~> Natural # f
      -> CotangentVector (Mean ~> Natural) f
{-# INLINE stochasticConditionalCrossEntropyDifferential #-}
stochasticConditionalCrossEntropyDifferential xs ys =
    stochasticConditionalCrossEntropyDifferential0
        (joinBoxedReplicated $ sufficientStatistic <$> xs)
        (joinBoxedReplicated $ sufficientStatistic <$> ys)

unnormalizedDensity :: forall m. ExponentialFamily m => Point Natural m -> SamplePoint m -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy :: Proxy m) x

exponentialFamilyDensity
    :: (ExponentialFamily m, Legendre Natural m)
    => Point Natural m -> SamplePoint m -> Double
{-# INLINE exponentialFamilyDensity #-}
exponentialFamilyDensity p x = unnormalizedDensity p x * (exp . negate $ potential p)

replicatedBaseMeasure0 :: (ExponentialFamily m, KnownNat k)
                       => Proxy m -> Proxy (Replicated k m) -> SamplePoint (Replicated k m) -> Double
{-# INLINE replicatedBaseMeasure0  #-}
replicatedBaseMeasure0 prxym _ xs = sum $ baseMeasure prxym <$> xs

sumBaseMeasure0
    :: (ExponentialFamily m, ExponentialFamily n)
    => Proxy m
    -> Proxy n
    -> Proxy (Sum m n)
    -> SamplePoint (Sum m n)
    -> Double
{-# INLINE sumBaseMeasure0 #-}
sumBaseMeasure0 prxym prxyn _ (xm,xn) =
     baseMeasure prxym xm * baseMeasure prxyn xn

---- | Estimates the dual relative entropy differential of an exponential family
---- distribution from another distribution, given the exponential family
---- distribution and unnormalized, log density.
--stochasticDualRelativeEntropyDifferential
--    :: (ClosedFormExponentialFamily x, Generative Natural x)
--    => Int
--    -> Natural # x
--    -> (Domain x -> Double)
--    -> RandST s (Differential # Tangent Natural x)
--stochasticDualRelativeEntropyDifferential n rx ldns = do
--    let (Harmonium xm _) = manifold hrm
--        (nx,nz,imtx) = splitHarmonium hrm
--        covariate x = ssx x <.> rx + ldns x
--    xs <- sample rx
--    let cv0 = average (\x -> covariate x .> sufficientStatistic x) <$> xs
--        cv1 = average [covariate x | x <- xs ] .> averagePoint [ ssx x | x <- xs ]
--        cv = cv0 <-> cv1
--    return . fromCoordinates (Tangent rx) $ coordinates cv


--- Conditional Distributions ---


-- | Applies the given conditional distribution to a sample from the 'SampleSpace' of the 'Domain'.
(>.>*) :: (Apply Mean c m, ExponentialFamily (Domain m))
       => Point (Mean ~> c) m
       -> SamplePoint (Domain m)
       -> Point c (Codomain m)
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application on samples.
(>$>*) :: (Apply Mean c m, ExponentialFamily (Domain m), KnownNat k)
       => Point (Mean ~> c) m
       -> Sample k (Domain m)
       -> Point c (Replicated k (Codomain m))
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> joinBoxedReplicated (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the given conditional distribution to a sample from the 'SampleSpace' of the 'Domain'.
(*<.<) :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f))
       => SamplePoint (Codomain f)
       -> Point (Function Mean Natural) f
       -> Point Natural (Domain f)
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped application on samples.
(*<$<) :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f), KnownNat k)
       => Sample k (Codomain f)
       -> Point (Function Mean Natural) f
       -> Point Natural (Replicated k (Domain f))
{-# INLINE (*<$<) #-}
(*<$<) xs p = joinBoxedReplicated (sufficientStatistic <$> xs) <$< p

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


instance (ExponentialFamily m, ExponentialFamily n) => ExponentialFamily (Sum m n) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm,xn) =
         joinSum (sufficientStatistic xm) (sufficientStatistic xn)
    {-# INLINE baseMeasure #-}
    baseMeasure = sumBaseMeasure0 Proxy Proxy

instance (ExponentialFamily m, KnownNat k) => ExponentialFamily (Replicated k m) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = joinBoxedReplicated $ sufficientStatistic <$> xs
    {-# INLINE baseMeasure #-}
    baseMeasure = replicatedBaseMeasure0 Proxy

instance (Manifold m, Manifold n, Transition Natural Source m, Transition Natural Source n) => Transition Natural Source (Sum m n) where
    {-# INLINE transition #-}
    transition pmn =
        let (pm,pn) = splitSum pmn
         in joinSum (transition pm) (transition pn)

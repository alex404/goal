{-# LANGUAGE
    TypeOperators,
    DataKinds,
    TypeFamilies,
    ConstraintKinds,
    FlexibleContexts,
    FlexibleInstances,
    MultiParamTypeClasses,
    ScopedTypeVariables,
    RankNTypes
#-}
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
    , stochasticCrossEntropyDifferential'
    , stochasticInformationProjectionDifferential
    -- *** Conditional Entropies
    , stochasticConditionalCrossEntropy
    , stochasticConditionalCrossEntropyDifferential
    , stochasticConditionalCrossEntropyDifferential0
    -- ** Conditional Distributions
    , (>.>*)
    , (>$>*)
    , (*<.<)
    , (*<$<)
    -- ** Model Selection
    , conditionalAkaikesInformationCriterion
    , conditionalBayesianInformationCriterion
    ) where

--- Imports ---


-- Package --

import Goal.Probability.Statistical

import Goal.Core
import Goal.Geometry

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
toNatural :: (Transition c Natural m) => Point c m -> Point Natural m
{-# INLINE toNatural #-}
toNatural = transition

-- | Expresses an exponential family distribution in mean coordinates.
toMean :: (Transition c Mean m) => Point c m -> Point Mean m
{-# INLINE toMean #-}
toMean = transition

-- | Expresses an exponential family distribution in source coordinates.
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
type ClosedFormExponentialFamily m = (ExponentialFamily m, Legendre Natural m, Legendre Mean m)

-- -- | The sufficient statistic of N iid random variables.
-- sufficientStatisticT
--     :: (ExponentialFamily m, KnownNat k, 1 <= k)
--     => Sample m -> Point Mean m
-- {-# INLINE sufficientStatisticT #-}
-- sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs
--
-- | The sufficient statistic of a traversable set of iid random variables.
sufficientStatisticT
    :: (ExponentialFamily m, Traversable f)
    => f (SamplePoint m) -> Point Mean m
{-# INLINE sufficientStatisticT #-}
sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs

-- | A function for computing the relative entropy, also known as the KL-divergence.
relativeEntropy :: ClosedFormExponentialFamily m => Mean # m -> Natural # m -> Double
{-# INLINE relativeEntropy #-}
relativeEntropy = divergence

-- | A function for computing the cross-entropy, which is the relative entropy plus the entropy of the first distribution.
crossEntropy :: Legendre Natural m => Mean # m -> Natural # m -> Double
{-# INLINE crossEntropy #-}
crossEntropy mp nq = potential nq - (mp <.> nq)

-- | The differential of the cross-entropy with respect to the parameters of the second argument.
crossEntropyDifferential :: Legendre Natural m => Mean # m -> Natural # m -> CotangentVector Natural m
{-# INLINE crossEntropyDifferential #-}
crossEntropyDifferential mp nq =
     potentialDifferential nq <-> primalIsomorphism mp

-- | An approximate cross-entropy based on samples from the first argument, and
-- an exact expression for the second argument.
stochasticCrossEntropy
    :: forall m . (ExponentialFamily m, Legendre Natural m)
    => Sample m -> Point Natural m -> Double
{-# INLINE stochasticCrossEntropy #-}
stochasticCrossEntropy xs nq =
    let mp = sufficientStatisticT xs
        bm = average $ log . baseMeasure (Proxy :: Proxy m) <$> xs
     in potential nq - (mp <.> nq) - bm

-- | An approximate cross-entropy differential based on samples from the first argument, and
-- an exact expression for differentiated distribution.
stochasticCrossEntropyDifferential
    :: (ExponentialFamily m, Legendre Natural m)
    => Sample m -> Point Natural m -> CotangentVector Natural m
{-# INLINE stochasticCrossEntropyDifferential #-}
stochasticCrossEntropyDifferential xs nq =
    let mp = sufficientStatisticT xs
     in potentialDifferential nq <-> primalIsomorphism mp

-- | The differential of the cross-entropy with respect to the parameters of the
-- second argument, based only on samples from the two distributions.
stochasticCrossEntropyDifferential'
    :: ExponentialFamily m
    => Sample m -- ^ True Samples
    -> Sample m -- ^ Model Samples
    -> CotangentVector Natural m -- ^ Differential Estimate
{-# INLINE stochasticCrossEntropyDifferential' #-}
stochasticCrossEntropyDifferential' pxs qxs =
    primalIsomorphism $ sufficientStatisticT qxs <-> sufficientStatisticT pxs


-- | The differential of the dual relative entropy.
stochasticInformationProjectionDifferential
    :: ExponentialFamily m
    => Natural # m -- ^ Model Distribution
    -> Sample m -- ^ Model Samples
    -> (SamplePoint m -> Double) -- ^ Unnormalized log-density of target distribution
    -> CotangentVector Natural m -- ^ Differential Estimate
{-# INLINE stochasticInformationProjectionDifferential #-}
stochasticInformationProjectionDifferential px xs f =
    let mxs = sufficientStatistic <$> xs
        mys = (\x -> sufficientStatistic x <.> px - f x) <$> xs
        ln = fromIntegral $ length xs
        mxht = ln /> foldr1 (<+>) mxs
        myht = sum mys / ln
        cvr = (ln - 1) /> foldr1 (<+>) [ (my - myht) .> (mx <-> mxht) | (mx,my) <- zip mxs mys ]
     in primalIsomorphism cvr

-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
stochasticConditionalCrossEntropy
    :: (Map Mean Natural f m n, ExponentialFamily n, ExponentialFamily m, Legendre Natural m)
    => Sample n -- ^ Input sample
    -> Sample m -- ^ Output sample
    -> Mean #> Natural # f m n -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE stochasticConditionalCrossEntropy #-}
stochasticConditionalCrossEntropy xs ys f =
    average . zipWith stochasticCrossEntropy ((:[]) <$> ys) $ f >$>* xs

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates (this is
-- primarily of internal use).
stochasticConditionalCrossEntropyDifferential0
    :: (Propagate Mean Natural f m n, ExponentialFamily n, Legendre Natural m)
    => [Mean # n] -- ^ Input mean distributions
    -> [Mean # m] -- ^ Output mean distributions
    -> Mean #> Natural # f m n -- ^ Function
    -> CotangentVector (Mean #> Natural) (f m n) -- ^ Differential
{-# INLINE stochasticConditionalCrossEntropyDifferential0 #-}
stochasticConditionalCrossEntropyDifferential0 xs ys f =
    let (df,yhts) = propagate mys xs f
        mys = dualIsomorphism <$> zipWith (<->) (potentialDifferential <$> yhts) (primalIsomorphism <$> ys)
     in primalIsomorphism df

-- | The stochastic conditional cross-entropy differential.
stochasticConditionalCrossEntropyDifferential
    :: (Propagate Mean Natural f m n, ExponentialFamily n, ExponentialFamily m, Legendre Natural m)
      => Sample n -- ^ Input Sample
      -> Sample m -- ^ Output sample
      -> Mean #> Natural # f m n -- ^ Parametric Function
      -> CotangentVector (Mean #> Natural) (f m n) -- ^ Function differential
{-# INLINE stochasticConditionalCrossEntropyDifferential #-}
stochasticConditionalCrossEntropyDifferential xs ys =
    stochasticConditionalCrossEntropyDifferential0 (sufficientStatistic <$> xs) (sufficientStatistic <$> ys)

-- | The unnormalized density of an arbitrary exponential family distribution.
unnormalizedDensity :: forall m. ExponentialFamily m => Point Natural m -> SamplePoint m -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy :: Proxy m) x

-- | The exact density of an exponential family distribution with an exact
-- expression for the log-partition function.
exponentialFamilyDensity
    :: (ExponentialFamily m, Legendre Natural m)
    => Point Natural m -> SamplePoint m -> Double
{-# INLINE exponentialFamilyDensity #-}
exponentialFamilyDensity p x = unnormalizedDensity p x * (exp . negate $ potential p)


--- Conditional Distributions ---


-- | Applies the given conditional distribution to a samplePoint.
(>.>*) :: (Map Mean c f m n, ExponentialFamily n)
       => Mean #> c # f m n
       -> SamplePoint n
       -> c # m
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application on samples.
(>$>*) :: (Map Mean c f m n, ExponentialFamily n)
       => Mean #> c # f m n
       -> Sample n
       -> [c # m]
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the transpose conditional distribution to a samplePoint.
(*<.<) :: (Bilinear f m n, ExponentialFamily m)
       => SamplePoint m
       -> Mean #> Natural # f m n
       -> Natural # n
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped application on samples.
(*<$<) :: (Bilinear f m n, ExponentialFamily m)
       => Sample m
       -> Mean #> Natural # f m n
       -> [Natural # n]
{-# INLINE (*<$<) #-}
(*<$<) xs p = (sufficientStatistic <$> xs) <$< p

infix 8 *<.<
infix 8 *<$<


--- Model Selection ---


-- | Calculate the conditional AIC for a given model and sample.
conditionalAkaikesInformationCriterion
    :: forall d f m n
    . (AbsolutelyContinuous d m, ExponentialFamily n, Map Mean d f m n)
    => Mean #> d # f m n
    -> Sample n
    -> Sample m
    -> Double
conditionalAkaikesInformationCriterion f xs ys =
    let d = natVal (Proxy :: Proxy (Dimension m))
        yhts = f >$>* xs
     in 2 * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]

---- | Calculate the conditional BIC for a given model and sample.
conditionalBayesianInformationCriterion
    :: forall d f m n
    . (AbsolutelyContinuous d m, ExponentialFamily n, Map Mean d f m n)
    => Mean #> d # f m n
    -> Sample n
    -> Sample m
    -> Double
conditionalBayesianInformationCriterion f xs ys =
    let d = natVal (Proxy :: Proxy (Dimension m))
        yhts = f >$>* xs
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]


--- Internal ---


replicatedBaseMeasure0 :: (ExponentialFamily m, KnownNat k)
                       => Proxy m -> Proxy (Replicated k m) -> SamplePoint (Replicated k m) -> Double
{-# INLINE replicatedBaseMeasure0  #-}
replicatedBaseMeasure0 prxym _ xs = sum $ baseMeasure prxym <$> xs

sumBaseMeasure
    :: (ExponentialFamily m, ExponentialFamily (Sum ms))
    => Proxy m
    -> Proxy (Sum ms)
    -> Proxy (Sum (m : ms))
    -> SamplePoint (Sum (m : ms))
    -> Double
{-# INLINE sumBaseMeasure #-}
sumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

pairBaseMeasure
    :: (ExponentialFamily m, ExponentialFamily n)
    => Proxy m
    -> Proxy n
    -> Proxy (m,n)
    -> SamplePoint (m,n)
    -> Double
{-# INLINE pairBaseMeasure #-}
pairBaseMeasure prxym prxyn _ (xm,xn) =
     baseMeasure prxym xm * baseMeasure prxyn xn


--- Instances ---


-- Replicated --

instance Transition Natural Natural m where
    {-# INLINE transition #-}
    transition = id

instance Transition Mean Mean m where
    {-# INLINE transition #-}
    transition = id

instance Transition Source Source m where
    {-# INLINE transition #-}
    transition = id

instance (ExponentialFamily m, KnownNat k) => ExponentialFamily (Replicated k m) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = joinBoxedReplicated $ sufficientStatistic <$> xs
    {-# INLINE baseMeasure #-}
    baseMeasure = replicatedBaseMeasure0 Proxy

-- Sum --

instance ExponentialFamily (Sum '[]) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic _ = zero
    {-# INLINE baseMeasure #-}
    baseMeasure _ _ = 1

instance (ExponentialFamily m, ExponentialFamily (Sum ms)) => ExponentialFamily (Sum (m : ms)) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm :+: xms) =
         joinSum (sufficientStatistic xm) (sufficientStatistic xms)
    {-# INLINE baseMeasure #-}
    baseMeasure = sumBaseMeasure Proxy Proxy

instance (ExponentialFamily m, ExponentialFamily n) => ExponentialFamily (m,n) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm,xn) =
         joinPair (sufficientStatistic xm) (sufficientStatistic xn)
    {-# INLINE baseMeasure #-}
    baseMeasure = pairBaseMeasure Proxy Proxy

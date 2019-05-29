{-# LANGUAGE TypeApplications #-}
-- | Definitions for working with exponential families.
module Goal.Probability.ExponentialFamily
    ( -- * Exponential Families
    ExponentialFamily (sufficientStatistic, baseMeasure)
    , ClosedFormExponentialFamily
    , sufficientStatisticT
    , unnormalizedDensity
    , unnormalizedLogDensity
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
class Manifold x => ExponentialFamily x s where
    sufficientStatistic :: s -> Mean # x
    baseMeasure :: Proxy x -> s -> Double

-- | When the 'Riemannian' properties of the given 'ExponentialFamily' may be
-- computed in closed-form, then we refer to it as a
-- 'ClosedFormExponentialFamily'.
type ClosedFormExponentialFamily x s = (ExponentialFamily x s, Legendre Natural x, Legendre Mean x)

-- -- | The sufficient statistic of N iid random variables.
-- sufficientStatisticT
--     :: (ExponentialFamily x, KnownNat k, 1 <= k)
--     => Sample x -> Mean # x
-- {-# INLINE sufficientStatisticT #-}
-- sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs
--
-- | The sufficient statistic of a traversable set of iid random variables.
sufficientStatisticT
    :: (ExponentialFamily x s, Traversable f)
    => f s -> Mean # x
{-# INLINE sufficientStatisticT #-}
sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (<+>) $ sufficientStatistic <$> xs

-- | A function for computing the relative entropy, also known as the KL-divergence.
relativeEntropy :: (Legendre Natural x, Legendre Mean x) => Mean # x -> Natural # x -> Double
{-# INLINE relativeEntropy #-}
relativeEntropy = divergence

-- | A function for computing the cross-entropy, which is the relative entropy plus the entropy of the first distribution.
crossEntropy :: Legendre Natural x => Mean # x -> Natural # x -> Double
{-# INLINE crossEntropy #-}
crossEntropy mp nq = potential nq - (mp <.> nq)

-- | The differential of the cross-entropy with respect to the parameters of the second argument.
crossEntropyDifferential :: Legendre Natural x => Mean # x -> Natural # x -> Mean # x
{-# INLINE crossEntropyDifferential #-}
crossEntropyDifferential mp nq =
     potentialDifferential nq <-> mp

-- | An approximate cross-entropy based on samples from the first argument, and
-- an exact expression for the second argument.
stochasticCrossEntropy
    :: forall x s . (ExponentialFamily x s, Legendre Natural x)
    => [s] -> Natural # x -> Double
{-# INLINE stochasticCrossEntropy #-}
stochasticCrossEntropy xs nq =
    let mp = sufficientStatisticT xs
        bm = average $ log . baseMeasure (Proxy :: Proxy x) <$> xs
     in potential nq - (mp <.> nq) - bm

-- | An approximate cross-entropy differential based on samples from the first argument, and
-- an exact expression for differentiated distribution.
stochasticCrossEntropyDifferential
    :: (ExponentialFamily x s, Legendre Natural x)
    => [s] -> Natural # x -> Mean # x
{-# INLINE stochasticCrossEntropyDifferential #-}
stochasticCrossEntropyDifferential xs nq =
    let mp = sufficientStatisticT xs
     in potentialDifferential nq <-> mp

-- | The differential of the cross-entropy with respect to the parameters of the
-- second argument, based only on samples from the two distributions.
stochasticCrossEntropyDifferential'
    :: ExponentialFamily x s
    => [s] -- ^ True Samples
    -> [s] -- ^ Model Samples
    -> Mean # x -- ^ Differential Estimate
{-# INLINE stochasticCrossEntropyDifferential' #-}
stochasticCrossEntropyDifferential' pxs qxs =
    sufficientStatisticT qxs <-> sufficientStatisticT pxs

-- | The differential of the dual relative entropy.
stochasticInformationProjectionDifferential
    :: ExponentialFamily x s
    => Natural # x -- ^ Model Distribution
    -> [s] -- ^ Model Samples
    -> (s -> Double) -- ^ Unnormalized log-density of target distribution
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

-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
stochasticConditionalCrossEntropy
    :: (Map Mean Natural f y x, ExponentialFamily x s, ExponentialFamily y t, Legendre Natural y)
    => [(t,s)] -- ^ (Output,Input) sample
    -> Mean #> Natural # f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE stochasticConditionalCrossEntropy #-}
stochasticConditionalCrossEntropy yxs f =
    let (ys,xs) = unzip yxs
     in average . zipWith stochasticCrossEntropy ((:[]) <$> ys) $ f >$>* xs

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates (this is
-- primarily of internal use).
stochasticConditionalCrossEntropyDifferential0
    :: (Propagate Mean Natural f y x, Manifold x, Manifold y, Legendre Natural y)
    => [(Mean # y, Mean # x)] -- ^ Output/Input mean distributions
    -> Mean #> Natural # f y x -- ^ Function
    -> Mean #> Natural #* f y x -- ^ Differential
{-# INLINE stochasticConditionalCrossEntropyDifferential0 #-}
stochasticConditionalCrossEntropyDifferential0 yxs f =
    let (ys,xs) = unzip yxs
        (df,yhts) = propagate mys xs f
        mys = zipWith (<->) (potentialDifferential <$> yhts) ys
     in df

-- | The stochastic conditional cross-entropy differential.
stochasticConditionalCrossEntropyDifferential
    :: (Propagate Mean Natural f y x, ExponentialFamily y t, ExponentialFamily x s, Legendre Natural y)
      => [(t,s)] -- ^ Output/Input Sample
      -> Mean #> Natural # f y x -- ^ Parametric Function
      -> Mean #> Natural #* f y x -- ^ Function differential
{-# INLINE stochasticConditionalCrossEntropyDifferential #-}
stochasticConditionalCrossEntropyDifferential yxs =
    let (ys,xs) = unzip yxs
     in stochasticConditionalCrossEntropyDifferential0 $ zip (sufficientStatistic <$> ys) (sufficientStatistic <$> xs)

-- | The unnormalized density of an arbitrary exponential family distribution.
unnormalizedDensity :: forall x s . ExponentialFamily x s => Natural # x -> s -> Double
unnormalizedDensity p x =
    exp (p <.> sufficientStatistic x) * baseMeasure (Proxy @ x) x

-- | The unnormalized log-density of an arbitrary exponential family distribution.
unnormalizedLogDensity :: forall x s . ExponentialFamily x s => Natural # x -> s -> Double
unnormalizedLogDensity p x =
    p <.> sufficientStatistic x  + log (baseMeasure (Proxy @ x) x)

-- | The exact density of an exponential family distribution with an exact
-- expression for the log-partition function.
exponentialFamilyDensity
    :: (ExponentialFamily x s, Legendre Natural x)
    => Natural # x -> s -> Double
{-# INLINE exponentialFamilyDensity #-}
exponentialFamilyDensity p x = unnormalizedDensity p x * (exp . negate $ potential p)


--- Conditional Distributions ---


-- | Applies the given conditional distribution to a samplePoint.
(>.>*) :: (Map Mean c f y x, ExponentialFamily x s)
       => Mean #> c # f y x
       -> s
       -> c # y
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application on samples.
(>$>*) :: (Map Mean c f y x, ExponentialFamily x s)
       => Mean #> c # f y x
       -> [s]
       -> [c # y]
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the transpose conditional distribution to a samplePoint.
(*<.<) :: (Bilinear f y x, ExponentialFamily y t)
       => t
       -> Mean #> Natural # f y x
       -> Natural # x
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped application on samples.
(*<$<) :: (Bilinear f y x, ExponentialFamily y t)
       => [t]
       -> Mean #> Natural # f y x
       -> [Natural # x]
{-# INLINE (*<$<) #-}
(*<$<) xs p = (sufficientStatistic <$> xs) <$< p

infix 8 *<.<
infix 8 *<$<


--- Model Selection ---


-- | Calculate the conditional AIC for a given model and sample.
conditionalAkaikesInformationCriterion
    :: forall d f x y s t
    . (AbsolutelyContinuous d y t, ExponentialFamily x s, Map Mean d f y x)
    => Mean #> d # f y x
    -> [(t,s)]
    -> Double
conditionalAkaikesInformationCriterion f yxs =
    let (ys,xs) = unzip yxs
        d = natVal (Proxy :: Proxy (Dimension y))
        yhts = f >$>* xs
     in 2 * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]

---- | Calculate the conditional BIC for a given model and sample.
conditionalBayesianInformationCriterion
    :: forall d f x y t s
    . (AbsolutelyContinuous d y t, ExponentialFamily x s, Map Mean d f y x)
    => Mean #> d # f y x
    -> [(t,s)]
    -> Double
conditionalBayesianInformationCriterion f yxs =
    let (ys,xs) = unzip yxs
        d = natVal (Proxy :: Proxy (Dimension y))
        yhts = f >$>* xs
        n = length xs
     in log (fromIntegral n) * fromIntegral d - 2 * sum
         [ log $ density yht y | (y,yht) <- zip ys yhts ]


--- Internal ---


replicatedBaseMeasure0 :: (ExponentialFamily x s, Storable s, KnownNat k)
                       => Proxy x -> Proxy (Replicated k x) -> S.Vector k s -> Double
{-# INLINE replicatedBaseMeasure0  #-}
replicatedBaseMeasure0 prxym _ xs = S.product $ S.map (baseMeasure prxym) xs

sumBaseMeasure
    :: (ExponentialFamily x s, ExponentialFamily (Sum xs) (HList ss))
    => Proxy x
    -> Proxy (Sum xs)
    -> Proxy (Sum (x : xs))
    -> HList (s : ss)
    -> Double
{-# INLINE sumBaseMeasure #-}
sumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

pairBaseMeasure
    :: (ExponentialFamily x s, ExponentialFamily y t)
    => Proxy x
    -> Proxy y
    -> Proxy (x,y)
    -> (s,t)
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

instance (ExponentialFamily x s, Storable s, KnownNat k)
  => ExponentialFamily (Replicated k x) (S.Vector k s) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = joinReplicated $ S.map sufficientStatistic xs
    {-# INLINE baseMeasure #-}
    baseMeasure = replicatedBaseMeasure0 Proxy

-- Sum --

instance ExponentialFamily (Sum '[]) (HList '[]) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic _ = zero
    {-# INLINE baseMeasure #-}
    baseMeasure _ _ = 1

instance (ExponentialFamily x s, ExponentialFamily (Sum xs) (HList ss))
  => ExponentialFamily (Sum (x : xs)) (HList (s : ss)) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm :+: xms) =
         joinSum (sufficientStatistic xm) (sufficientStatistic xms)
    {-# INLINE baseMeasure #-}
    baseMeasure = sumBaseMeasure Proxy Proxy

instance (ExponentialFamily x s, ExponentialFamily y t) => ExponentialFamily (x,y) (s,t) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic (xm,xn) =
         joinPair (sufficientStatistic xm) (sufficientStatistic xn)
    {-# INLINE baseMeasure #-}
    baseMeasure = pairBaseMeasure Proxy Proxy

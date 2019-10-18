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
    , relativeEntropyDifferential
    , stochasticRelativeEntropyDifferential
    , stochasticInformationProjectionDifferential
    -- ** Conditional Distributions
    , (>.>*)
    , (>$>*)
    , (*<.<)
    , (*<$<)
    , conditionalLogLikelihood
    , conditionalLogLikelihoodDifferential
    , sortedConditionalLogLikelihood
    , sortedConditionalLogLikelihoodDifferential
    -- *** Maximum Likelihood Instances
    , exponentialFamilyLogLikelihood
    , exponentialFamilyLogLikelihoodDifferential
    ) where

--- Imports ---


-- Package --

import Goal.Probability.Statistical

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import qualified Data.Map.Strict as M
import Foreign.Storable
import Data.Tuple

--- Exponential Families ---

-- Source Chart --

-- | A parameterization which represents the standard or typical parameterization of
-- the given manifold, e.g. the 'Poisson' rate or 'Normal' mean and standard deviation.
data Source


-- | A parameterization in terms of the natural parameters of an exponential family.
data Natural

-- | A parameterization in terms of the mean 'sufficientStatistic' of an exponential family.
data Mean

instance Primal Natural where
    type Dual Natural = Mean

instance Primal Mean where
    type Dual Mean = Natural

-- | Expresses an exponential family distribution in 'Natural' coordinates.
toNatural :: (Transition c Natural x) => c # x -> Natural # x
{-# INLINE toNatural #-}
toNatural = transition

-- | Expresses an exponential family distribution in 'Mean' coordinates.
toMean :: (Transition c Mean x) => c # x -> Mean # x
{-# INLINE toMean #-}
toMean = transition

-- | Expresses an exponential family distribution in 'Source' coordinates.
toSource :: (Transition c Source x) => c # x -> Source # x
{-# INLINE toSource #-}
toSource = transition


-- | An 'ExponentialFamily' is a 'Statistical' 'Manifold' \( \mathcal M \)
-- determined by a fixed-length 'sufficientStatistic' \(s_i\) and a
-- 'baseMeasure' \(\mu\). Each distribution \(P \in \mathcal M\) may then be
-- identified with 'Natural' parameters \(\theta_i\) such that
-- \(p(x) \propto e^{\sum_{i=1}^n \theta_i s_i(x)}\mu(x)\).  'ExponentialFamily'
-- distributions theoretically have a 'Riemannian' geometry, with 'metric'
-- 'Tensor' given by the Fisher information metric. However, not all
-- distributions (e.g. the von Mises distribution) afford closed-form
-- expressions for all the relevant structures.
class Statistical x => ExponentialFamily x where
    sufficientStatistic :: SamplePoint x -> Mean # x
    baseMeasure :: Proxy x -> SamplePoint x -> Double

-- | When the log-partition function and its derivative of the given
-- 'ExponentialFamily' may be computed in closed-form, then we refer to it as a
-- 'LegendreExponentialFamily'.
--
-- Note that the log-partition function is the 'potential' of the 'Legendre'
-- class, and its derivative maps 'Natural' coordinates to 'Mean' coordinates.
type LegendreExponentialFamily x =
    ( PotentialCoordinates x ~ Natural, Legendre x, ExponentialFamily x
    , Transition (PotentialCoordinates x) (Dual (PotentialCoordinates x)) x )

-- | When additionally, the (negative) entropy and its derivative of the given
-- 'ExponentialFamily' may be computed in closed-form, then we refer to it as a
-- 'DuallyFlatExponentialFamily'.
--
-- Note that the negative entropy is the 'dualPotential' of the 'DuallyFlat' class,
-- and its derivative maps 'Mean' coordinates to 'Natural'' coordinates.
type DuallyFlatExponentialFamily x =
    ( LegendreExponentialFamily x, DuallyFlat x
    , Transition (Dual (PotentialCoordinates x)) (PotentialCoordinates x) x )

-- | The sufficient statistic of a traversable set of iid random variables.
sufficientStatisticT
    :: (ExponentialFamily x, Traversable f)
    => f (SamplePoint x) -> Mean # x
{-# INLINE sufficientStatisticT #-}
sufficientStatisticT xs = (fromIntegral (length xs) />) . foldr1 (+) $ sufficientStatistic <$> xs

-- | The relative entropy \(D(P \parallel Q)\), also known as the KL-divergence.
-- This is simply the 'canonicalDivergence' with its arguments flipped.
relativeEntropy :: DuallyFlatExponentialFamily x => Mean # x -> Natural # x -> Double
{-# INLINE relativeEntropy #-}
relativeEntropy = flip canonicalDivergence

-- | A function for computing the cross-entropy, which is the relative entropy
-- plus the entropy of the first distribution.
crossEntropy :: DuallyFlatExponentialFamily x => Mean # x -> Natural # x ->
    Double
{-# INLINE crossEntropy #-}
crossEntropy mp nq = potential nq - (mp <.> nq)

-- | The differential of the relative entropy with respect to the 'Natural' parameters of
-- the second argument.
relativeEntropyDifferential :: LegendreExponentialFamily x => Mean # x -> Natural # x -> Mean # x
{-# INLINE relativeEntropyDifferential #-}
relativeEntropyDifferential mp nq = transition nq - mp

-- | Monte Carlo estimate of the differential of the relative entropy with
-- respect to the 'Natural' parameters of the second argument, based on samples from
-- the two distributions.
stochasticRelativeEntropyDifferential
    :: ExponentialFamily x
    => Sample x -- ^ True Samples
    -> Sample x -- ^ Model Samples
    -> Mean # x -- ^ Differential Estimate
{-# INLINE stochasticRelativeEntropyDifferential #-}
stochasticRelativeEntropyDifferential pxs qxs =
    sufficientStatisticT qxs - sufficientStatisticT pxs

-- | Estimate of the differential of relative entropy with respect to the
-- 'Natural' parameters of the first argument, based a 'Sample' from the first
-- argument and the unnormalized log-density of the second.
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
        mxht = ln /> foldr1 (+) mxs
        myht = sum mys / ln
     in (ln - 1) /> foldr1 (+) [ (my - myht) .> (mx - mxht) | (mx,my) <- zip mxs mys ]

-- | The density of an exponential family distribution that has an exact
-- expression for the log-partition function.
exponentialFamilyDensity
    :: LegendreExponentialFamily x => Natural # x -> SamplePoint x -> Double
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

-- | 'logLikelihood' for a 'LegendreExponentialFamily'.
exponentialFamilyLogLikelihood
    :: forall x . LegendreExponentialFamily x
    => Sample x -> Natural # x -> Double
{-# INLINE exponentialFamilyLogLikelihood #-}
exponentialFamilyLogLikelihood xs nq =
    let mp = sufficientStatisticT xs
        bm = average $ log . baseMeasure (Proxy :: Proxy x) <$> xs
     in -potential nq + (mp <.> nq) + bm

-- | 'logLikelihoodDifferential' for a 'LegendreExponentialFamily'.
exponentialFamilyLogLikelihoodDifferential
    :: LegendreExponentialFamily x
    => Sample x -> Natural # x -> Mean # x
{-# INLINE exponentialFamilyLogLikelihoodDifferential #-}
exponentialFamilyLogLikelihoodDifferential xs nq =
    let mp = sufficientStatisticT xs
     in mp - transition nq


--- Conditional Distributions ---


dependantLogLikelihood
    :: (LogLikelihood d y s, Map Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Double
{-# INLINE dependantLogLikelihood #-}
dependantLogLikelihood ysxs chrm =
    let (yss,xs) = unzip ysxs
     in average . zipWith logLikelihood yss $ chrm >$> xs

dependantLogLikelihoodDifferential
    :: (LogLikelihood d y s, Propagate Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Function Mean d #* f y x
{-# INLINE dependantLogLikelihoodDifferential #-}
dependantLogLikelihoodDifferential ysxs chrm =
    let (yss,xs) = unzip ysxs
        (df,yhts) = propagate mys xs chrm
        mys = zipWith logLikelihoodDifferential yss yhts
     in df

-- | The conditional 'logLikelihood' for a conditional distribution.
conditionalLogLikelihood
    :: (ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t)
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE conditionalLogLikelihood #-}
conditionalLogLikelihood yxs f =
    let ysxs = [ ([y],sufficientStatistic x) | (y,x) <- yxs ]
     in dependantLogLikelihood ysxs f

-- | The conditional 'logLikelihoodDifferential' for a conditional distribution.
conditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t, Propagate Mean Natural f y x )
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
{-# INLINE conditionalLogLikelihoodDifferential #-}
conditionalLogLikelihoodDifferential yxs f =
    let ysxs = [ ([y],sufficientStatistic x) | (y,x) <- yxs ]
     in dependantLogLikelihoodDifferential ysxs f

-- | The conditional 'logLikelihood' for a conditional distribution, where
-- redundant conditions/inputs are combined. This can dramatically increase performance when
-- the number of distinct conditions/inputs is small.
sortedConditionalLogLikelihood
    :: ( ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t )
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE sortedConditionalLogLikelihood #-}
sortedConditionalLogLikelihood yxs f =
    let ysxs = map swap . M.toList
            $ M.fromListWith (++) [(sufficientStatistic x, [y]) | (y, x) <- yxs]
     in dependantLogLikelihood ysxs f

-- | The conditional 'logLikelihoodDifferential', where redundant conditions are
-- combined. This can dramatically increase performance when the number of
-- distinct conditions is small.
sortedConditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t
       , Propagate Mean Natural f y x, Ord (SamplePoint x) )
    => [(t, SamplePoint x)] -- ^ Output/Input pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
{-# INLINE sortedConditionalLogLikelihoodDifferential #-}
sortedConditionalLogLikelihoodDifferential yxs f =
    let ysxs = map swap . M.toList
            $ M.fromListWith (++) [(sufficientStatistic x, [y]) | (y, x) <- yxs]
     in dependantLogLikelihoodDifferential ysxs f


-- | Evalutes the given conditional distribution at a 'SamplePoint'.
(>.>*) :: (Map Mean c f y x, ExponentialFamily x)
       => Function Mean c # f y x
       -> SamplePoint x
       -> c # y
{-# INLINE (>.>*) #-}
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application of conditional distributions on a 'Sample'.
(>$>*) :: (Map Mean c f y x, ExponentialFamily x)
       => Function Mean c # f y x
       -> Sample x
       -> [c # y]
{-# INLINE (>$>*) #-}
(>$>*) p xs = p >$> (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the transpose of a 'Bilinear' 'Map' to a 'SamplePoint'.
(*<.<) :: (Map Mean Natural f x y, Bilinear f y x, ExponentialFamily y)
       => SamplePoint y
       -> Natural #> f y x
       -> Natural # x
{-# INLINE (*<.<) #-}
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped transpose application on a 'Sample'.
(*<$<) :: (Map Mean Natural f x y, Bilinear f y x, ExponentialFamily y)
       => Sample y
       -> Natural #> f y x
       -> [Natural # x]
{-# INLINE (*<$<) #-}
(*<$<) xs p = (sufficientStatistic <$> xs) <$< p

infix 8 *<.<
infix 8 *<$<


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
    sufficientStatistic _ = 0
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

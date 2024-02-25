{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | Various instances of statistical manifolds, with a focus on exponential
families. In the documentation we use \(X\) to indicate a random variable
with the distribution being documented.
-}
module Goal.Probability.Distributions.Gaussian (
    -- * Manifolds
    StandardNormal,
    CovarianceMatrix,
    KnownCovariance,
    MultivariateNormal,
    Normal,
    FullNormal,
    DiagonalNormal,
    IsotropicNormal,

    -- * Construction
    splitNaturalNormal,
    joinNaturalNormal,
    standardNormal,

    -- * Analysis
    bivariateNormalConfidenceEllipse,
    multivariateNormalCorrelations,

    -- * Linear Models
    LinearModel,
    FullLinearModel,
    FactorAnalysis,
    PrincipleComponentAnalysis,
) where

--- Imports ---

--- Goal

import Goal.Core
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.Statistical

import Goal.Geometry

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

import System.Random.MWC.Distributions qualified as R

--- Misc

import Control.Monad (replicateM)
import Data.Proxy (Proxy (..))

-- Normal Distribution --

{- | The Mean of a normal distribution. When used as a distribution itself, it
is a Normal distribution with unit variance.
-}
data StandardNormal (n :: Nat)

-- | The variance of a normal distribution.
type CovarianceMatrix t n = Linear t (StandardNormal n) (StandardNormal n)

-- | Synonym for known positive definite covariance matrices.
type KnownCovariance t n = KnownLinear t (StandardNormal n) (StandardNormal n)

--- Multivariate Normal ---

{- | The 'Manifold' of 'MultivariateNormal' distributions. The 'Source'
coordinates are the (vector) mean and the covariance matrix. For the
coordinates of a multivariate normal distribution, the elements of the mean
come first, and then the elements of the covariance matrix in row major
order.

Note that we only store the lower triangular elements of the covariance
matrix, to better reflect the true dimension of a MultivariateNormal
Manifold. In short, be careful when using 'join' and 'split' to access the
values of the Covariance matrix, and consider using the specific instances
for MVNs.
-}
type MultivariateNormal t (n :: Nat) = LocationShape (StandardNormal n) (CovarianceMatrix t n)

type Normal = MultivariateNormal L.PositiveDefinite 1
type FullNormal n = MultivariateNormal L.PositiveDefinite n
type DiagonalNormal n = MultivariateNormal L.Diagonal n
type IsotropicNormal n = MultivariateNormal L.Scale n

-- | Linear models are linear functions with additive Guassian noise.
type LinearModel f n k = Affine L.Full (StandardNormal n) (MultivariateNormal f n) (StandardNormal k)

type FullLinearModel n k = LinearModel L.PositiveDefinite n k
type FactorAnalysis n k = LinearModel L.Diagonal n k
type PrincipleComponentAnalysis n k = LinearModel L.Scale n k

splitNaturalNormal ::
    (KnownCovariance t n) =>
    Natural # MultivariateNormal t n ->
    (Natural # StandardNormal n, Natural # CovarianceMatrix t n)
splitNaturalNormal mvn =
    let (mu, sgma) = split mvn
     in (mu, precisionPreCorrection sgma)

joinNaturalNormal ::
    (KnownCovariance t n) =>
    Natural # StandardNormal n ->
    Natural # CovarianceMatrix t n ->
    Natural # MultivariateNormal t n
joinNaturalNormal mu sgma =
    join mu $ precisionPostCorrection sgma

bivariateNormalConfidenceEllipse ::
    ( KnownCovariance L.PositiveDefinite 2
    , KnownCovariance t 2
    ) =>
    Int ->
    Double ->
    Source # MultivariateNormal t 2 ->
    [(Double, Double)]
bivariateNormalConfidenceEllipse nstps sds mvn =
    let (mu, sgma) = split mvn
        pd :: Source # CovarianceMatrix L.PositiveDefinite 2
        pd = fromTensor $ toTensor sgma
        mrt = sds .> choleskyDecomposition pd
        xs = range 0 (2 * pi) nstps
        sxs = [fromTuple (cos x, sin x) | x <- xs]
     in S.toPair . coordinates . (mu +) <$> mrt >$> sxs

-- | Create a standard normal distribution in a variety of forms
standardNormal ::
    forall c t n.
    (KnownCovariance t n, Transition Source c (MultivariateNormal t n)) =>
    c # MultivariateNormal t n
standardNormal =
    let sgm0 :: Source # CovarianceMatrix t n
        sgm0 = identity
     in transition $ join 0 sgm0

-- | Computes the correlation matrix of a 'MultivariateNormal' distribution.
multivariateNormalCorrelations ::
    forall t n.
    (KnownCovariance t n) =>
    Source # MultivariateNormal t n ->
    Source # Tensor (StandardNormal n) (StandardNormal n)
multivariateNormalCorrelations mvn =
    let cvrs = toTensor . snd $ split mvn
        diag :: Source # Diagonal (StandardNormal n)
        diag = fromTensor cvrs
        sds = breakManifold $ sqrt diag
        sdmtx = sds >.< sds
     in cvrs / sdmtx

--- Internal ---

standardNormalLogBaseMeasure ::
    forall n.
    (KnownNat n) =>
    Proxy (StandardNormal n) ->
    S.Vector n Double ->
    Double
standardNormalLogBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n / 2 * log (2 * pi)

multivariateNormalLogBaseMeasure ::
    forall f n.
    (KnownNat n) =>
    Proxy (MultivariateNormal f n) ->
    S.Vector n Double ->
    Double
multivariateNormalLogBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n / 2 * log (2 * pi)

-- | Halve the off diagonal elements of a triangular matrix.
preCorrect :: (KnownNat n) => S.Vector n Double -> S.Vector n Double
preCorrect trng = S.triangularMapDiagonal (* 2) $ trng / 2

-- | Double the off diagonal elements of a triangular matrix.
postCorrect :: (KnownNat n) => S.Vector n Double -> S.Vector n Double
postCorrect trng = S.triangularMapDiagonal (/ 2) $ trng * 2

precisionPreCorrection0 :: forall t n. (KnownNat n) => L.Linear t n n -> L.Linear t n n
{-# INLINE precisionPreCorrection0 #-}
precisionPreCorrection0 f@(L.PositiveDefiniteLinear _) =
    L.PositiveDefiniteLinear . preCorrect $ L.toVector f
precisionPreCorrection0 m = m

precisionPostCorrection0 :: forall t n. (KnownNat n) => L.Linear t n n -> L.Linear t n n
{-# INLINE precisionPostCorrection0 #-}
precisionPostCorrection0 f@(L.PositiveDefiniteLinear _) =
    L.PositiveDefiniteLinear . postCorrect $ L.toVector f
precisionPostCorrection0 m = m

precisionPreCorrection ::
    (KnownCovariance t n) =>
    Natural # CovarianceMatrix t n ->
    Natural # CovarianceMatrix t n
precisionPreCorrection = Point . L.toVector . precisionPreCorrection0 . useLinear

precisionPostCorrection ::
    (KnownCovariance t n) =>
    Natural # CovarianceMatrix t n ->
    Natural # CovarianceMatrix t n
precisionPostCorrection = Point . L.toVector . precisionPostCorrection0 . useLinear

{- | samples a multivariateNormal by way of a covariance matrix i.e. by taking
the square root.
-}
sampleFullNormal ::
    (KnownCovariance L.PositiveDefinite n) =>
    Int ->
    Source # FullNormal n ->
    Random [S.Vector n Double]
sampleFullNormal n p = do
    let (mu, sgma) = split p
        rtsgma = choleskyDecomposition sgma
    x0s <- replicateM n . S.replicateM $ Random (R.normal 0 1)
    return $ coordinates . (mu +) <$> rtsgma >$> (Point <$> x0s)

sampleDiagonalNormal ::
    (KnownCovariance L.Diagonal n) =>
    Int ->
    Source # DiagonalNormal n ->
    Random [S.Vector n Double]
sampleDiagonalNormal n p = do
    let (mu, sgma) = split p
        rtsgma = sqrt sgma
    x0s <- replicateM n . S.replicateM $ Random (R.normal 0 1)
    return $ coordinates . (mu +) <$> rtsgma >$> (Point <$> x0s)

sampleScaleNormal ::
    (KnownCovariance L.Scale n) =>
    Int ->
    Source # IsotropicNormal n ->
    Random [S.Vector n Double]
sampleScaleNormal n p = do
    let (mu, sgma) = split p
        rtsgma = sqrt sgma
    x0s <- replicateM n . S.replicateM $ Random (R.normal 0 1)
    return $ coordinates . (mu +) <$> rtsgma >$> (Point <$> x0s)

--- Instances ---

--- Standard Normal ---

type instance PotentialCoordinates (StandardNormal n) = Natural

instance (KnownNat n) => Manifold (StandardNormal n) where
    type Dimension (StandardNormal n) = n

instance (KnownNat n) => Statistical (StandardNormal n) where
    type SamplePoint (StandardNormal n) = S.Vector n Double

instance (KnownNat n) => ExponentialFamily (StandardNormal n) where
    sufficientStatistic = Point
    logBaseMeasure = standardNormalLogBaseMeasure

instance Transition Natural Mean (StandardNormal n) where
    transition = breakChart

instance (KnownNat n) => Legendre (StandardNormal n) where
    potential p = 0.5 * (p <.> toMean p)

instance
    ( ExponentialFamily (StandardNormal n)
    , Transition Natural Mean (StandardNormal n)
    , Legendre (StandardNormal n)
    ) =>
    LogLikelihood Natural (StandardNormal n) (S.Vector n Double)
    where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

--- MultivariateNormal ---

--- Transition Functions

type instance PotentialCoordinates (MultivariateNormal t n) = Natural

instance
    (KnownCovariance t n) =>
    Transition Source Natural (MultivariateNormal t n)
    where
    transition p =
        let (mu, sgma) = split p
            invsgma = inverse sgma
         in joinNaturalNormal (breakChart $ invsgma >.> mu) . breakChart $ (-2) /> invsgma

instance
    (KnownCovariance t n) =>
    Transition Natural Source (MultivariateNormal t n)
    where
    transition p =
        let (nmu, nsgma) = splitNaturalNormal p
            insgma = inverse $ (-2) .> nsgma
         in join (breakChart $ insgma >.> nmu) $ breakChart insgma

instance
    (KnownCovariance t n) =>
    Transition Source Mean (MultivariateNormal t n)
    where
    transition mvn =
        let (mu, sgma) = split mvn
            mmvn = join mu $ sgma + (mu >.< mu)
         in breakChart mmvn

instance
    (KnownCovariance t n) =>
    Transition Mean Source (MultivariateNormal t n)
    where
    transition mmvn =
        let (mu, msgma) = split mmvn
            mvn = join mu $ msgma - (mu >.< mu)
         in breakChart mvn

instance
    ( KnownCovariance t n
    , Transition Source Mean (MultivariateNormal t n)
    ) =>
    Transition Natural Mean (MultivariateNormal t n)
    where
    transition = toMean . toSource

instance
    ( KnownCovariance t n
    , Transition Mean Source (MultivariateNormal t n)
    ) =>
    Transition Mean Natural (MultivariateNormal t n)
    where
    transition = toNatural . toSource

--- Basic Instances

instance
    (KnownCovariance t n) =>
    AbsolutelyContinuous Source (MultivariateNormal t n)
    where
    densities mvn xs =
        let (mu, sgma) = split mvn
            n = fromIntegral $ natValInt (Proxy @n)
            scl = (2 * pi) ** (-n / 2) * determinant sgma ** (-1 / 2)
            dffs = [Point $ x - coordinates mu | x <- xs]
            expvals = zipWith (<.>) dffs $ inverse sgma >$> dffs
         in (scl *) . exp . negate . (/ 2) <$> expvals

instance
    (KnownCovariance L.Scale n, Transition c Source (IsotropicNormal n)) =>
    Generative c (IsotropicNormal n)
    where
    sample n = sampleScaleNormal n . toSource

instance
    (KnownCovariance L.Diagonal n, Transition c Source (DiagonalNormal n)) =>
    Generative c (DiagonalNormal n)
    where
    sample n = sampleDiagonalNormal n . toSource

instance
    (KnownCovariance L.PositiveDefinite n, Transition c Source (FullNormal n)) =>
    Generative c (FullNormal n)
    where
    sample n = sampleFullNormal n . toSource

--- Exponential Family Instances

instance
    (KnownCovariance t n) =>
    ExponentialFamily (MultivariateNormal t n)
    where
    sufficientStatistic x =
        let mx = sufficientStatistic x
         in join mx $ mx >.< mx
    averageSufficientStatistic xs =
        let mxs = sufficientStatistic <$> xs
         in join (average mxs) $ mxs >$< mxs
    logBaseMeasure = multivariateNormalLogBaseMeasure

instance (KnownCovariance t n) => Legendre (MultivariateNormal t n) where
    potential p =
        let (nmu, nsgma) = splitNaturalNormal p
            (insgma, lndt, _) = inverseLogDeterminant . negate $ 2 * nsgma
         in 0.5 * (nmu <.> (insgma >.> nmu)) - 0.5 * lndt

instance
    ( KnownCovariance t n
    , Transition Mean Source (MultivariateNormal t n)
    ) =>
    DuallyFlat (MultivariateNormal t n)
    where
    dualPotential p =
        let sgma = snd . split $ toSource p
            n = natValInt (Proxy @n)
            lndet0 = log $ determinant sgma
            lndet = fromIntegral n * log (2 * pi * exp 1) + lndet0
         in -0.5 * lndet

instance
    ( ExponentialFamily (MultivariateNormal t n)
    , Transition Natural Mean (MultivariateNormal t n)
    , Legendre (MultivariateNormal t n)
    ) =>
    LogLikelihood Natural (MultivariateNormal t n) (S.Vector n Double)
    where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance
    ( ExponentialFamily (MultivariateNormal t n)
    , Transition Natural Mean (MultivariateNormal t n)
    , Legendre (MultivariateNormal t n)
    ) =>
    AbsolutelyContinuous Natural (MultivariateNormal t n)
    where
    logDensities = exponentialFamilyLogDensities

instance
    ( ExponentialFamily (MultivariateNormal t n)
    , Transition Mean c (MultivariateNormal t n)
    ) =>
    MaximumLikelihood c (MultivariateNormal t n)
    where
    mle = transition . averageSufficientStatistic

--- Linear Models ---

instance
    (KnownCovariance t n, KnownNat n, KnownNat k) =>
    Transition Natural Source (LinearModel t n k)
    where
    transition nlm =
        let (nmvn, nmtx) = split nlm
            smvn = toSource nmvn
            sgma = snd $ split smvn
            smtx = sgma <#> breakChart nmtx
         in join smvn smtx

instance
    (KnownCovariance t n, KnownNat n, KnownNat k) =>
    Transition Source Natural (LinearModel t n k)
    where
    transition sfa =
        let (smvn, smtx) = split sfa
            nmvn = toNatural smvn
            nvr = snd $ split nmvn
            nmtx = -2 .> breakChart (breakChart nvr <#> smtx)
         in join nmvn nmtx

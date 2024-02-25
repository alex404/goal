{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | Various instances of statistical manifolds, with a focus on exponential
families. In the documentation we use \(X\) to indicate a random variable
with the distribution being documented.
-}
module Goal.Probability.Distributions.Boltzmann where

--     -- * Manifolds
--     StandardNormal,
--     CovarianceMatrix,
--     KnownCovariance,
--     MultivariateNormal,
--     Normal,
--     FullNormal,
--     DiagonalNormal,
--     IsotropicNormal,
--
--     -- * Construction
--     splitNaturalNormal,
--     joinNaturalNormal,
--     standardNormal,
--
--     -- * Analysis
--     bivariateNormalConfidenceEllipse,
--     multivariateNormalCorrelations,
--
--     -- * Linear Models
--     LinearModel,
--     FullLinearModel,
--     FactorAnalysis,
--     PrincipleComponentAnalysis,
-- ) where
--

--- Imports ---

--- Goal

import Goal.Core
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.Statistical

import Goal.Geometry

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Control.Monad (foldM, replicateM)
import Data.Finite (Finite, natToFinite)
import Data.Maybe (fromJust)
import Data.Proxy (Proxy (..))
import Data.Vector qualified as V
import System.Random.MWC.Distributions qualified as R

-- Boltzmann Distribution

-- | The Boltzmann distribution is a probability distribution over the set of all \(2^n\) possible states of a system with \(n\) binary degrees of freedom.
data Boltzmann (n :: Nat)

-- | The variance of a normal distribution.
type InteractionMatrix n = Symmetric (Replicated n Bernoulli)

--- Functions

-- | Convert a Boltzmann distribution to an interaction matrix. Note that the diagonal includes the self-interaction/bias terms.
boltzmannToInteractionMatrix :: (KnownNat n) => Natural # Boltzmann n -> Natural # InteractionMatrix n
boltzmannToInteractionMatrix = preCorrection . breakManifold

-- | Convert an interaction matrix to a Boltzmann distribution. Note that the diagonal includes the self-interaction/bias terms.
interactionMatrixToBoltzmann :: (KnownNat n) => Natural # InteractionMatrix n -> Natural # Boltzmann n
interactionMatrixToBoltzmann = breakManifold . postCorrection

boltzmannBiases :: (KnownNat n) => Natural # Boltzmann n -> Natural # Replicated n Bernoulli
boltzmannBiases = Point . S.triangularTakeDiagonal . coordinates

-- | The Gibbs sampling algorithm for a Boltzmann distribution.
gibbsBoltzmann :: forall n. (1 <= n, KnownNat n) => Int -> Natural # Boltzmann n -> Random (S.Vector n Bool)
gibbsBoltzmann ncycs bltz = do
    let prr :: Natural # Replicated n Bernoulli
        prr = 0
    bls <- samplePoint prr
    iterateM' ncycs (cycleBoltzmann bltz) bls

-- | The probability of a single unit being active in a Boltzmann distribution.
unitDistribution :: (KnownNat n) => Natural # Boltzmann n -> S.Vector n Bool -> Finite n -> Natural # Bernoulli
unitDistribution bltz bls idx =
    let blstru = bls S.// [(idx, True)]
        blsfls = bls S.// [(idx, False)]
     in singleton $ bltz <.> (sufficientStatistic blstru - sufficientStatistic blsfls)

--- Internal ---

-- | A single step of the Gibbs sampling algorithm for a Boltzmann distribution.
stepBoltzmann ::
    (KnownNat n) =>
    Natural # Boltzmann n ->
    S.Vector n Bool ->
    Finite n ->
    Random (S.Vector n Bool)
stepBoltzmann bltz bls idx = do
    let brn = unitDistribution bltz bls idx
    bl <- samplePoint brn
    return $ bls S.// [(idx, bl)]

-- | A single step of the Gibbs sampling algorithm for a Boltzmann distribution.
cycleBoltzmann ::
    forall n.
    (KnownNat n, 1 <= n) =>
    Natural # Boltzmann n ->
    S.Vector n Bool ->
    Random (S.Vector n Bool)
cycleBoltzmann bltz bls = do
    let n1 = natToFinite (Proxy @(n - 1))
    -- idxs :: V.Vector Int <- Random (R.uniformPermutation n)
    let idxs0 = V.fromList [0 .. n1]
    idxs <- Random (R.uniformShuffle idxs0)
    foldM (stepBoltzmann bltz) bls idxs

--- Functions

-- | Halve the off diagonal elements of a triangular matrix.
preCorrect :: (KnownNat n) => S.Vector n Double -> S.Vector n Double
preCorrect trng = S.triangularMapDiagonal (* 2) $ trng / 2

-- | Double the off diagonal elements of a triangular matrix.
postCorrect :: (KnownNat n) => S.Vector n Double -> S.Vector n Double
postCorrect trng = S.triangularMapDiagonal (/ 2) $ trng * 2

preCorrection0 :: forall n. (KnownNat n) => L.Linear L.Symmetric n n -> L.Linear L.Symmetric n n
{-# INLINE preCorrection0 #-}
preCorrection0 f = L.SymmetricLinear . preCorrect $ L.toVector f

postCorrection0 :: forall n. (KnownNat n) => L.Linear L.Symmetric n n -> L.Linear L.Symmetric n n
{-# INLINE postCorrection0 #-}
postCorrection0 f =
    L.SymmetricLinear . postCorrect $ L.toVector f

preCorrection ::
    (KnownNat n) =>
    Natural # InteractionMatrix n ->
    Natural # InteractionMatrix n
preCorrection = Point . L.toVector . preCorrection0 . useLinear

postCorrection ::
    (KnownNat n) =>
    Natural # InteractionMatrix n ->
    Natural # InteractionMatrix n
postCorrection = Point . L.toVector . postCorrection0 . useLinear

--- Instances

type instance PotentialCoordinates (Boltzmann n) = Natural

instance (KnownNat n) => Manifold (Boltzmann n) where
    type Dimension (Boltzmann n) = Triangular n

instance (KnownNat n) => Statistical (Boltzmann n) where
    type SamplePoint (Boltzmann n) = S.Vector n Bool

instance (KnownNat n) => Discrete (Boltzmann n) where
    type Cardinality (Boltzmann n) = 2 ^ n
    sampleSpace _ =
        fromJust . S.fromList
            <$> replicateM (natValInt @n Proxy) [False, True]

instance (KnownNat n) => ExponentialFamily (Boltzmann n) where
    sufficientStatistic bls =
        let bls' = S.map (fromIntegral . fromEnum) bls
         in Point . S.lowerTriangular $ S.outerProduct bls' bls'
    logBaseMeasure _ _ = 0

instance (KnownNat n) => Legendre (Boltzmann n) where
    potential p =
        let blss = pointSampleSpace p
         in logSumExp . dotMap p $ sufficientStatistic <$> blss

instance
    ( ExponentialFamily (Boltzmann n)
    , Transition Natural Mean (Boltzmann n)
    , Legendre (Boltzmann n)
    ) =>
    AbsolutelyContinuous Natural (Boltzmann n)
    where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n) => Transition Natural Mean (Boltzmann n) where
    transition bltz =
        let blss = pointSampleSpace bltz
            blss' = S.map (fromIntegral . fromEnum) <$> blss
            prbs = densities bltz blss
         in -- in weightedAverage . zip prbs $ sufficientStatistic <$> blss
            Point . S.lowerTriangular . S.weightedAverageOuterProduct $ zip3 prbs blss' blss'

instance
    ( ExponentialFamily (Boltzmann n)
    , Transition Natural Mean (Boltzmann n)
    , Legendre (Boltzmann n)
    ) =>
    LogLikelihood Natural (Boltzmann n) (S.Vector n Bool)
    where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

--- Instances 2

type Boltzmann2 (n :: Nat) = LocationShape (Replicated n Bernoulli) (Symmetric (Replicated n Bernoulli))

type instance PotentialCoordinates (Boltzmann2 n) = Natural

instance (KnownNat n) => Statistical (Boltzmann2 n) where
    type SamplePoint (Boltzmann2 n) = S.Vector n Bool

instance (KnownNat n) => Discrete (Boltzmann2 n) where
    type Cardinality (Boltzmann2 n) = 2 ^ n
    sampleSpace _ =
        fromJust . S.fromList
            <$> replicateM (natValInt @n Proxy) [False, True]

instance (KnownNat n) => ExponentialFamily (Boltzmann2 n) where
    sufficientStatistic bls =
        let bls' = S.map (fromIntegral . fromEnum) bls
         in join (Point bls') . Point . S.lowerTriangular $ S.outerProduct bls' bls'
    logBaseMeasure _ _ = 0

-- instance (KnownNat n) => Legendre (Boltzmann n) where
--     potential p =
--         let blss = pointSampleSpace p
--          in logSumExp . dotMap p $ sufficientStatistic <$> blss
--
-- instance
--     ( ExponentialFamily (Boltzmann n)
--     , Transition Natural Mean (Boltzmann n)
--     , Legendre (Boltzmann n)
--     ) =>
--     AbsolutelyContinuous Natural (Boltzmann n)
--     where
--     logDensities = exponentialFamilyLogDensities
--
-- instance (KnownNat n) => Transition Natural Mean (Boltzmann n) where
--     transition bltz =
--         let blss = pointSampleSpace bltz
--             blss' = S.map (fromIntegral . fromEnum) <$> blss
--             prbs = densities bltz blss
--          in -- in weightedAverage . zip prbs $ sufficientStatistic <$> blss
--             Point . S.lowerTriangular . S.weightedAverageOuterProduct $ zip3 prbs blss' blss'
--
-- instance
--     ( ExponentialFamily (Boltzmann n)
--     , Transition Natural Mean (Boltzmann n)
--     , Legendre (Boltzmann n)
--     ) =>
--     LogLikelihood Natural (Boltzmann n) (S.Vector n Bool)
--     where
--     logLikelihood = exponentialFamilyLogLikelihood
--     logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

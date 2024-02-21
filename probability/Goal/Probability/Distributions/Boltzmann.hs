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

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

import System.Random.MWC.Distributions qualified as R

--- Misc

import Control.Monad (replicateM)
import Data.Maybe (fromJust)
import Data.Proxy (Proxy (..))

-- Normal Distribution --

{- | The Mean of a normal distribution. When used as a distribution itself, it
is a Normal distribution with unit variance.
-}
data Boltzmann (n :: Nat)

--- Internal ---

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
         in Point . S.lowerTriangular . S.weightedAverageOuterProduct $
                zip3 prbs blss' blss'

instance
    ( ExponentialFamily (Boltzmann n)
    , Transition Natural Mean (Boltzmann n)
    , Legendre (Boltzmann n)
    ) =>
    LogLikelihood Natural (Boltzmann n) (S.Vector n Bool)
    where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

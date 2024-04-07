{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | Various instances of statistical manifolds, with a focus on exponential
families. In the documentation we use \(X\) to indicate a random variable
with the distribution being documented.
-}
module Goal.Probability.Distributions.Boltzmann (
    Boltzmann,
    KnownBoltzmann,
    InteractionMatrix (..),
    gibbsBoltzmann,
    unitDistribution,
) where

--- Imports ---

--- Goal

import Goal.Core
import Goal.Probability.Distributions
import Goal.Probability.Distributions.Gaussian
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

-- | The variance of a normal distribution.
newtype InteractionMatrix n = InteractionMatrix (Symmetric (Replicated (n - 1) Bernoulli))

deriving instance (KnownNat n, 1 <= n) => Manifold (InteractionMatrix n)

type instance CovarianceRep (InteractionMatrix n) = L.Symmetric

{- | The Boltzmann distribution is a probability distribution over the set of all \(2^n\) possible states of a system with \(n\) binary degrees of freedom.
type Boltzmann (n :: Nat) = MomentParameters (Replicated n Bernoulli) (Symmetric (Replicated n Bernoulli))
-}
type Boltzmann (n :: Nat) = MomentParameters (Replicated n Bernoulli) (InteractionMatrix n)

type KnownBoltzmann t k =
    ( GeneralizedGaussian Natural (Replicated k Bernoulli) (InteractionMatrix k)
    , KnownNat k
    , 1 <= k
    )

--- Functions

-- | The Gibbs sampling algorithm for a Boltzmann distribution.
gibbsBoltzmann :: forall n. (1 <= n, KnownNat n) => Int -> Natural # Boltzmann n -> Random (S.Vector n Bool)
{-# INLINE gibbsBoltzmann #-}
gibbsBoltzmann ncycs bltz = do
    let prr :: Natural # Replicated n Bernoulli
        prr = 0
    bls <- samplePoint prr
    iterateM' ncycs (cycleBoltzmann bltz) bls

-- | The probability of a single unit being active in a Boltzmann distribution.
unitDistribution :: (1 <= n, KnownNat n) => Natural # Boltzmann n -> S.Vector n Bool -> Finite n -> Natural # Bernoulli
{-# INLINE unitDistribution #-}
unitDistribution bltz bls idx =
    let blstru = bls S.// [(idx, True)]
        blsfls = bls S.// [(idx, False)]
     in singleton $ bltz <.> (sufficientStatistic blstru - sufficientStatistic blsfls)

--- Internal ---

-- | A single step of the Gibbs sampling algorithm for a Boltzmann distribution.
stepBoltzmann ::
    (KnownNat n, 1 <= n) =>
    Natural # Boltzmann n ->
    S.Vector n Bool ->
    Finite n ->
    Random (S.Vector n Bool)
{-# INLINE stepBoltzmann #-}
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
{-# INLINE cycleBoltzmann #-}
cycleBoltzmann bltz bls = do
    let n1 = natToFinite (Proxy @(n - 1))
    -- idxs :: V.Vector Int <- Random (R.uniformPermutation n)
    let idxs0 = V.fromList [0 .. n1]
    idxs <- Random (R.uniformShuffle idxs0)
    foldM (stepBoltzmann bltz) bls idxs

--- Instances

type instance PotentialCoordinates (Boltzmann n) = Natural

instance (KnownNat n) => Discrete (Boltzmann n) where
    type Cardinality (Boltzmann n) = 2 ^ n
    sampleSpace _ =
        fromJust . S.fromList
            <$> replicateM (natValInt @n Proxy) [False, True]

instance (KnownNat n, 1 <= n) => ExponentialFamily (Boltzmann n) where
    sufficientStatistic bls =
        let bls' = S.map (fromIntegral . fromEnum) bls
            intmtx = S.lowerTriangular $ S.outerProduct bls' bls'
            (mu, cvr) = S.triangularSplitDiagonal intmtx
         in join (Point mu) (Point cvr)
    logBaseMeasure _ _ = 0

instance (KnownNat n, 1 <= n) => Legendre (Boltzmann n) where
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

instance (KnownNat n, 1 <= n) => Transition Natural Mean (Boltzmann n) where
    {-# INLINE transition #-}
    transition bltz =
        let blss = pointSampleSpace bltz
            blss' = S.map (fromIntegral . fromEnum) <$> blss
            prbs = densities bltz blss
            intmtx = S.lowerTriangular . S.weightedAverageOuterProduct $ zip3 prbs blss' blss'
            (mu, cvr) = S.triangularSplitDiagonal intmtx
         in join (Point mu) (Point cvr)

instance
    ( ExponentialFamily (Boltzmann n)
    , Transition Natural Mean (Boltzmann n)
    , Legendre (Boltzmann n)
    ) =>
    LogLikelihood Natural (Boltzmann n) (S.Vector n Bool)
    where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance (KnownNat n, 1 <= n) => Generative Natural (Boltzmann n) where
    {-# INLINE sample #-}
    sample n bltz = do
        x0 <- gibbsBoltzmann brnn bltz
        xs <- iterateM (n - 1) (iterateM' skp (cycleBoltzmann bltz)) x0
        return $ x0 : xs
      where
        brnn = 1000
        skp = 10

instance
    (KnownNat n, 1 <= n) =>
    GeneralizedGaussian Mean (Replicated n Bernoulli) (InteractionMatrix n)
    where
    splitGaussian bltz =
        let (mu, symm) = split bltz
            symm' = coordinates symm
            cvr0 = S.triangularJoinDiagonal (coordinates mu) symm'
         in (0, Point cvr0)
    joinGaussian mu0 cvr0 =
        let (mu', symm) = S.triangularSplitDiagonal $ coordinates cvr0
            mu = mu0 + Point mu'
         in join mu (Point symm)

instance
    (KnownNat n, 1 <= n) =>
    GeneralizedGaussian Natural (Replicated n Bernoulli) (InteractionMatrix n)
    where
    splitGaussian bltz =
        let (mu, symm) = split bltz
            symm' = coordinates symm / 2
            cvr0 = S.triangularJoinDiagonal (coordinates mu) symm'
         in (0, Point cvr0)
    joinGaussian mu0 cvr0 =
        let (mu', symm') = S.triangularSplitDiagonal $ coordinates cvr0
            mu = mu0 + Point mu'
            symm = 2 * Point symm'
         in join mu symm

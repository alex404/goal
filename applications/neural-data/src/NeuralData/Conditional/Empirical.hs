{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeOperators
    #-}

module NeuralData.Conditional.Empirical
    ( -- * Empirical Analysis
      empiricalTuningCurves
    , subsampleEmpiricalTuningCurves
    , empiricalPosterior0
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import qualified Numeric.Sum as NS

-- Qualified --

import qualified Data.Map as M


--- Inference ---


-- Under the assumption of a flat prior
unnormalizedEmpiricalPPCLogPosterior0
    :: KnownNat k
    => Bool
    -> Mean # Neurons k
    -> Response k
    -> Double
unnormalizedEmpiricalPPCLogPosterior0 True mz z =
     toNatural mz <.> sufficientStatistic z - NS.sum NS.kbn (listCoordinates mz)
unnormalizedEmpiricalPPCLogPosterior0 False mz z =
     toNatural mz <.> sufficientStatistic z

-- Under the assumption of a flat prior
empiricalPosterior0
    :: (Ord s, KnownNat k)
    => Bool -> M.Map s (Mean # Neurons k) -> Response k -> M.Map s Double
empiricalPosterior0 nrmb xzmp z =
    let uldns = flip (unnormalizedEmpiricalPPCLogPosterior0 nrmb) z <$> xzmp
        avg = NS.sum NS.kbn uldns / fromIntegral (length uldns)
        udns = exp . subtract avg <$> uldns
        nrm = traceGiven $ NS.sum NS.kbn udns
     in (/nrm) <$> udns

empiricalTuningCurves :: (Ord s, KnownNat k) => M.Map s [Response k] -> M.Map s (Mean # Neurons k)
empiricalTuningCurves zxmp = sufficientStatisticT <$> zxmp

subsampleEmpiricalTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subsampleEmpiricalTuningCurves nzxmp idxs =
     Point . flip S.backpermute idxs . coordinates <$> nzxmp

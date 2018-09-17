{-# LANGUAGE GADTs,ScopedTypeVariables,DataKinds,TypeOperators #-}

module NeuralData
    ( -- * Types
    NeuralModel
    , Neurons
    , Response
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    , getPopulationSize
    -- * Processing
    , stimulusResponseMap
    , empiricalTuningCurves
    , empiricalPPCPosterior0
    -- * Subsampling
    , generateIndices
    , subSampleTuningCurves
    , subSampleResponses
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Numeric.Sum as NS

-- Qualified --

import qualified Data.Map as M


--- Types ---


type NeuralModel s k = Harmonium Tensor (Replicated k Poisson) s
type Neurons k = Replicated k Poisson
type Response k = SamplePoint (Neurons k)


--- IO ---

-- Under the assumption of a flat prior
unnormalizedEmpiricalPPCLogPosterior0
    :: KnownNat k
    => Bool
    -> Mean # Neurons k
    -> Response k
    -> Double
unnormalizedEmpiricalPPCLogPosterior0 True mz z =
     dualTransition mz <.> sufficientStatistic z - NS.sum NS.kbn (listCoordinates mz)
unnormalizedEmpiricalPPCLogPosterior0 False mz z =
     dualTransition mz <.> sufficientStatistic z

-- Under the assumption of a flat prior
empiricalPPCPosterior0
    :: (Ord s, KnownNat k)
    => Bool -> M.Map s (Mean # Neurons k) -> Response k -> M.Map s Double
empiricalPPCPosterior0 nrmb xzmp z =
    let uldns = flip (unnormalizedEmpiricalPPCLogPosterior0 nrmb) z <$> xzmp
        avg = NS.sum NS.kbn uldns / fromIntegral (length uldns)
        udns = exp . subtract avg <$> uldns
        nrm = traceGiven $ NS.sum NS.kbn udns
     in (/nrm) <$> udns

getNeuralData :: Read s => Collection -> Dataset -> IO [([Int], s)]
getNeuralData (Collection clcstr) (Dataset dststr) =
    read <$> readFile (clcstr ++ "/data/" ++ dststr ++ ".dat")

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . B.fromList <$> ks) ss

getPopulationSize :: [([Int], s)] -> Int
getPopulationSize = length . fst . head



--- Processing ---


stimulusResponseMap :: Ord s => [(Response k, s)] -> M.Map s [Response k]
stimulusResponseMap zxs = M.fromListWith (++) [(x,[z]) | (z,x) <- zxs]

empiricalTuningCurves :: (Ord s, KnownNat k) => M.Map s [Response k] -> M.Map s (Mean # Neurons k)
empiricalTuningCurves zxmp = sufficientStatisticT <$> zxmp


--- Subsampling ---


subSampleTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subSampleTuningCurves nzxmp idxs =
     Point . flip S.backpermute idxs . coordinates <$> nzxmp

subSampleResponses
    :: (Ord s, KnownNat k)
    => M.Map s [Response k]
    -> B.Vector k' Int
    -> M.Map s [Response k']
subSampleResponses zxmp idxs =
     map (`B.backpermute` idxs) <$> zxmp

generateIndices
    :: forall k k' r . (KnownNat k, KnownNat k', k' <= k)
    => Proxy k
    -> Random r (B.Vector k' Int)
generateIndices _ = do
    let idxs :: B.Vector k Int
        idxs = B.generate finiteInt
    subsampleVector idxs

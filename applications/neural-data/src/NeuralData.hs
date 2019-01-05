{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeApplications,
    TypeOperators
    #-}

module NeuralData
    ( -- * Types
      prjnm
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    -- * General
    , meanSDInliers
    -- * Subsampling
    , generateIndices
    , subSampleResponses
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B

-- Qualified --

import qualified Data.Map as M


--- Types ---


--- CSV ---

prjnm :: String
prjnm = "neural-data"

--- Inference ---


getNeuralData :: Read s => String -> String -> IO (NatNumber,[([Int], s)])
getNeuralData expnm dst = read <$> goalReadDataset (Experiment prjnm expnm) dst

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . B.fromList <$> ks) ss


--- Analysis ---


meanSDInliers :: [Double] -> (Double,Double)
{-# INLINE meanSDInliers #-}
meanSDInliers xs =
    let (mu,vr) = estimateMeanVariance xs
        xs' = filter (\x -> square (x-mu) < 4*vr) xs
        (mu',vr') = estimateMeanVariance xs'
     in (mu',sqrt vr')

subSampleResponses
    :: (Ord s, KnownNat k, KnownNat m)
    => M.Map s [Response (k+m)]
    -> B.Vector k Int
    -> M.Map s [Response k]
subSampleResponses zxmp idxs =
     map (`B.backpermute` idxs) <$> zxmp

generateIndices
    :: forall k m r . (KnownNat k, KnownNat m)
    => Proxy (k + m)
    -> Random r (B.Vector k Int)
{-# INLINE generateIndices #-}
generateIndices _ = do
    let idxs :: B.Vector (k + m) Int
        idxs = B.generate finiteInt
    subsampleVector idxs

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Flags

import Criterion.Main qualified as C
import Criterion.Types (Config (..))

--- Globals ---

--- Data

f :: S.Vector 1 Double -> S.Vector 1 Double
f x = exp . sin $ 2 * x

fp :: Source # Normal
fp = fromTuple (0, 0.1)

mnx, mxx :: S.Vector 1 Double
mnx = -3
mxx = 3

smlxs, mdmxs, lrgxs :: [S.Vector 1 Double]
smlxs = [0]
mdmxs = range mnx mxx 10
lrgxs = range mnx mxx 100

--- Neural Networks

type SmallNeuralNetwork =
    NeuralNetwork L.Full '[ '(L.Full, R 20 Bernoulli)] (StandardNormal 1) (Euclidean 1)

type MediumNeuralNetwork =
    NeuralNetwork L.Full '[ '(L.Full, R 100 Bernoulli), '(L.Full, R 100 Bernoulli)] (StandardNormal 1) (Euclidean 1)

type LargeNeuralNetwork =
    NeuralNetwork L.Full '[ '(L.Full, R 1000 Bernoulli), '(L.Full, R 1000 Bernoulli)] (StandardNormal 1) (Euclidean 1)

--- Benchmarks

smallName, mediumName, largeName :: String
smallName = "NN: 1x20x1"
mediumName = "NN: 1x100x100x1"
largeName = "NN: 1x1000x1000x1"

-- Record over different test operators
data NeuralNetworks = NeuralNetworks
    { small :: Natural # SmallNeuralNetwork
    , medium :: Natural # MediumNeuralNetwork
    , large :: Natural # LargeNeuralNetwork
    }

forwardPass :: NeuralNetworks -> [S.Vector 1 Double] -> [C.Benchmark]
forwardPass nns xs =
    [ C.bench smallName $ C.nf forwardFun (small nns)
    , C.bench mediumName $ C.nf forwardFun (medium nns)
    , C.bench largeName $ C.nf forwardFun (large nns)
    ]
  where
    forwardFun nn = nn >$>* xs

backprop :: NeuralNetworks -> [S.Vector 1 Double] -> [S.Vector 1 Double] -> [C.Benchmark]
backprop nns xs ys =
    [ C.bench smallName $ C.nf (conditionalLogLikelihoodDifferential (zip ys xs)) (small nns)
    , C.bench mediumName $ C.nf (conditionalLogLikelihoodDifferential (zip ys xs)) (medium nns)
    , C.bench largeName $ C.nf (conditionalLogLikelihoodDifferential (zip ys xs)) (large nns)
    ]

--- Main ---

main :: IO ()
main = do
    smlys <- realize $ mapM (noisyFunction fp f) smlxs
    mdmys <- realize $ mapM (noisyFunction fp f) mdmxs
    lrgys <- realize $ mapM (noisyFunction fp f) lrgxs

    smlmlp :: Natural # SmallNeuralNetwork <- realize $ uniformInitialize (-1, 1)
    mdmmlp :: Natural # MediumNeuralNetwork <- realize $ uniformInitialize (-1, 1)
    lrgmlp :: Natural # LargeNeuralNetwork <- realize $ uniformInitialize (-1, 1)

    let nns = NeuralNetworks smlmlp mdmmlp lrgmlp

    bnchfl <- benchFilePath "backprop.html"
    C.defaultMainWith
        C.defaultConfig
            { reportFile = Just bnchfl
            }
        [ C.bgroup "Forward, Batch 1" (forwardPass nns smlxs)
        , C.bgroup "Forward, Batch 10" (forwardPass nns mdmxs)
        , C.bgroup "Forward, Batch 100" (forwardPass nns lrgxs)
        , C.bgroup "Backward, Batch 1" (backprop nns smlxs smlys)
        , C.bgroup "Backward, Batch 10" (backprop nns mdmxs mdmys)
        , C.bgroup "Backward, Batch 100" (backprop nns lrgxs lrgys)
        ]

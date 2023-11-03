{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Globals ---

--- Data

f :: S.Vector 1 Double -> S.Vector 1 Double
f x = exp . sin $ 2 * x

fp :: Source # Normal
fp = fromTuple (0, 0.1)

mnx, mxx :: S.Vector 1 Double
mnx = -3
mxx = 3

xs :: [S.Vector 1 Double]
xs = concat . replicate 3 $ range mnx mxx 16

--- Neural Network

type LocalNeuralNetwork =
    NeuralNetwork L.Full '[ '(L.Full, R 20 Bernoulli)] (StandardNormal 1) (StandardNormal 1)

--- Training

nepchs :: Int
nepchs = 1000

eps :: Double
eps = 0.03

--- Momentum
mxmu :: Double
mxmu = 0.999

--- Plot

pltrng :: [S.Vector 1 Double]
pltrng = range mnx mxx 1000

finalLineFun :: Natural # LocalNeuralNetwork -> [Double]
finalLineFun mlp = S.head . coordinates <$> mlp >$>* pltrng

--- Main ---

main :: IO ()
main = do
    ys <- realize $ mapM (noisyFunction fp f) xs

    let xys = zip ys xs

    let cost :: Natural # LocalNeuralNetwork -> Double
        cost = conditionalLogLikelihood xys

    let backprop :: Natural # LocalNeuralNetwork -> Natural #* LocalNeuralNetwork
        backprop = conditionalLogLikelihoodDifferential xys

    mlp0 <- realize $ uniformInitialize (-1, 1)

    let gdmlps =
            take nepchs
                $ vanillaGradientSequence backprop eps Classic mlp0
        mtmmlps =
            take nepchs
                $ vanillaGradientSequence backprop eps (defaultMomentumPursuit mxmu) mlp0
        admmlps =
            take nepchs
                $ vanillaGradientSequence backprop eps defaultAdamPursuit mlp0

    let gdln = finalLineFun $ last gdmlps
        mtmln = finalLineFun $ last mtmmlps
        admln = finalLineFun $ last admmlps

    let fxs = S.head . f <$> pltrng
    let gdcsts = cost <$> gdmlps
        mtmcsts = cost <$> mtmmlps
        admcsts = cost <$> admmlps

    let jsonData =
            toJSON
                [ "xys" .= zip (S.head <$> xs) (S.head <$> ys)
                , "inputs" .= (S.head <$> pltrng)
                , "true-outputs" .= fxs
                , "gradient-descent-outputs" .= gdln
                , "gradient-descent-training" .= gdcsts
                , "momentum-outputs" .= mtmln
                , "momentum-training" .= mtmcsts
                , "adam-outputs" .= admln
                , "adam-training" .= admcsts
                ]

    rsltfl <- resultsFilePath "neural-network.json"
    exportJSON rsltfl jsonData

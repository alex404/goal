{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified Criterion.Main as C


--- Globals ---


-- Data --

f :: Double -> Double
f x = exp . sin $ 2 * x

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: [Double]
xs = concat . replicate 5 $ range mnx mxx 8

fp :: Source # Normal
fp = Point $ S.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ S.doubleton 0 0.1

type NeuralNetwork' =
    NeuralNetwork '[ '(Tensor, R 50 Bernoulli)]
    Tensor NormalMean NormalMean

-- Training --

nepchs :: Int
nepchs = 1000

eps :: Double
eps = 0.01

-- Momentum
mxmu :: Double
mxmu = 0.999


--- Main ---


main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let xys = zip ys xs

    let cost :: Natural # NeuralNetwork' -> Double
        cost = conditionalLogLikelihood xys

    let backprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        backprop = conditionalLogLikelihoodDifferential xys

    let sortedBackprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        sortedBackprop = mapConditionalLogLikelihoodDifferential $ conditionalDataMap xys

        sgdmlps0 mlp = take nepchs $ mlp0 : vanillaGradientSequence backprop eps Classic mlp
        mtmmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps (defaultMomentumPursuit mxmu) mlp
        admmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps defaultAdamPursuit mlp
        sadmmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence sortedBackprop eps defaultAdamPursuit mlp

    C.defaultMain
       [ C.bench "sgd" $ C.nf sgdmlps0 mlp0
       , C.bench "momentum" $ C.nf mtmmlps0 mlp0
       , C.bench "adam" $ C.nf admmlps0 mlp0
       , C.bench "sorted-adam" $ C.nf sadmmlps0 mlp0 ]

    let sgdmlps = sgdmlps0 mlp0
        mtmmlps = mtmmlps0 mlp0
        admmlps = admmlps0 mlp0
        sadmmlps = sadmmlps0 mlp0

    let sgdcst = cost $ last sgdmlps
        mtmcst = cost $ last mtmmlps
        admcst = cost $ last admmlps
        sadmcst = cost $ last sadmmlps

    putStrLn "SGD LL:"
    print sgdcst
    putStrLn "Momentum LL:"
    print mtmcst
    putStrLn "Adam LL:"
    print admcst
    putStrLn "Sorted Adam LL:"
    print sadmcst

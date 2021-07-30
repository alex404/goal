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
xs = range mnx mxx 200

fp :: Source # Normal
fp = Point $ S.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ S.doubleton 0 0.0001

type NeuralNetwork' =
    NeuralNetwork '[ '(Tensor, R 1000 Bernoulli), '(Tensor, R 1000 Bernoulli)]
    Tensor NormalMean NormalMean



-- Training --

nepchs :: Int
nepchs = 1

eps :: Double
eps = 0.0001

-- Layout --

main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let xys = zip ys xs

    let cost :: Natural # NeuralNetwork' -> Double
        cost = conditionalLogLikelihood xys

    let backprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        backprop = conditionalLogLikelihoodDifferential xys

        admmlps0 mlp = take nepchs $ vanillaGradientSequence backprop eps defaultAdamPursuit mlp

    let mlp = last $!! admmlps0 mlp0

    C.defaultMain
       [ C.bench "application" $ C.nf cost mlp
       , C.bench "backpropagation" $ C.nf backprop mlp ]

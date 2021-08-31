#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE DeriveGeneric,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import qualified Data.List as L


--- Globals ---


-- Data --

f :: Double -> Double
f x = exp . sin $ 2 * x

fp :: Source # Normal
fp = fromTuple (0,0.1)

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: [Double]
xs = concat . replicate 5 $ range mnx mxx 8

-- Neural Network --

cp :: Source # Normal
cp = fromTuple (0,0.1)

type NeuralNetwork' =
    NeuralNetwork '[ '(Tensor, R 50 Bernoulli)]
    Tensor NormalMean NormalMean

-- Training --

nepchs :: Int
nepchs = 1000

eps :: Double
eps = 0.05

-- Momentum
mxmu :: Double
mxmu = 0.999

-- Plot --

pltrng :: [Double]
pltrng = range mnx mxx 1000

finalLineFun :: Natural # NeuralNetwork' -> [Double]
finalLineFun mlp = S.head . coordinates <$> mlp >$>* pltrng


-- CSV --



--- Main ---


main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    let xys = zip ys xs

    let cost :: Natural # NeuralNetwork' -> Double
        cost = conditionalLogLikelihood xys

    let backprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        backprop = conditionalLogLikelihoodDifferential xys

    mlp0 <- realize $ initialize cp

    let sgdmlps = take nepchs $ vanillaGradientSequence
            backprop eps Classic mlp0
        mtmmlps = take nepchs $ vanillaGradientSequence
            backprop eps (defaultMomentumPursuit mxmu) mlp0
        admmlps = take nepchs $ vanillaGradientSequence
            backprop eps defaultAdamPursuit mlp0

    let sgdln = finalLineFun $ last sgdmlps
        mtmln = finalLineFun $ last mtmmlps
        admln = finalLineFun $ last admmlps

    let smpcsv = zip xs ys

    let fxs = f <$> pltrng
    let rgcsv = L.zip5 pltrng fxs sgdln mtmln admln

    let sgdcst = cost <$> sgdmlps
        mtmcst = cost <$> mtmmlps
        admcst = cost <$> admmlps

    let cstcsv = L.zip3 sgdcst mtmcst admcst

    let ldpth = "."

    goalExport ldpth "samples" smpcsv
    goalExport ldpth "regression" rgcsv
    goalExport ldpth "gradient-ascent" cstcsv

    runGnuplot ldpth "regression-lines"
    runGnuplot ldpth "log-likelihood-ascent"

#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE DeriveGeneric,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Data.List


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

--type Layer1 = MeanNormal (1/1)
--type Layer2 = R 50 Bernoulli
--type Layer3 = MeanNormal (1/1)
--type NeuralNetwork' = NeuralNetworkLayer (Affine Tensor) (Affine Tensor) Layer2 Layer3 Layer1

--type NeuralNetwork' =
--    HiddenNeuralNetwork
--    [Affine Tensor, Affine Tensor]
--    '[R 1000 Bernoulli]
--    (MeanNormal (1/1)) (MeanNormal (1/1))

type NeuralNetwork' = NeuralNetwork
        '[ '((<*),R 50 Bernoulli)] (<*) (MeanNormal (1/1)) (MeanNormal (1/1))

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

data LogLikelihoodAscent = LogLikelihoodAscent
    { sga :: Double
    , momentum :: Double
    , adam :: Double
    , mapAdam :: Double }
    deriving (Generic, Show)

instance ToNamedRecord LogLikelihoodAscent where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered LogLikelihoodAscent where
    headerOrder = goalCSVOrder

data RegressionSamples = RegressionSamples
    { xSample :: Double
    , ySample :: Double }
    deriving (Generic, Show)

instance ToNamedRecord RegressionSamples
instance DefaultOrdered RegressionSamples

data RegressionLines = RegressionLines
    { input :: Double
    , sgaMean :: Double
    , momentumMean :: Double
    , adamMean :: Double
    , mapAdamMean :: Double }
    deriving (Generic, Show)

instance ToNamedRecord RegressionLines where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered RegressionLines where
    headerOrder = goalCSVOrder



--- Main ---


main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let xys = zip ys xs
        xymp = conditionalDataMap xys

    let cost :: Natural # NeuralNetwork' -> Double
        cost = mapConditionalLogLikelihood xymp

    let backprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        backprop = conditionalLogLikelihoodDifferential xys

    let mapBackprop :: Natural # NeuralNetwork' -> Natural #* NeuralNetwork'
        mapBackprop = mapConditionalLogLikelihoodDifferential xymp

        sgdmlps0 mlp = take nepchs $ mlp0 : vanillaGradientSequence backprop eps Classic mlp
        mtmmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps (defaultMomentumPursuit mxmu) mlp
        admmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps defaultAdamPursuit mlp
        sadmmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence mapBackprop eps defaultAdamPursuit mlp

    let sgdmlps = sgdmlps0 mlp0
        mtmmlps = mtmmlps0 mlp0
        admmlps = admmlps0 mlp0
        sadmmlps = sadmmlps0 mlp0

    let sgdln = finalLineFun $ last sgdmlps
        mtmln = finalLineFun $ last mtmmlps
        admln = finalLineFun $ last admmlps
        sadmln = finalLineFun $ last sadmmlps

    let smpcsv = zipWith RegressionSamples xs ys

    let rgcsv = zipWith5 RegressionLines pltrng sgdln mtmln admln sadmln

    let sgdcst = cost <$> sgdmlps
        mtmcst = cost <$> mtmmlps
        admcst = cost <$> admmlps
        sadmcst = cost <$> admmlps

    let cstcsv = zipWith4 LogLikelihoodAscent sgdcst mtmcst admcst sadmcst

    let ldpth = "."

    goalExportNamed ldpth "samples" smpcsv
    goalExportNamed ldpth "regression" rgcsv
    goalExportNamed ldpth "gradient-ascent" cstcsv

    runGnuplot ldpth "regression-lines"
    runGnuplot ldpth "log-likelihood-ascent"

{-# LANGUAGE DeriveGeneric,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Data.List

-- Qualified --

import qualified Criterion.Main as C


--- Globals ---


expnm :: String
expnm = "regression"

-- Data --

f :: Double -> Double
f x = exp . sin $ 2 * x

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: [Double]
xs = range mnx mxx 20

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
    [Affine Tensor, Affine Tensor] [MeanNormal (1/1), R 50 Bernoulli, MeanNormal (1/1)]

-- Training --

nepchs :: Int
nepchs = 1000

eps :: Double
eps = -0.05

-- Momentum
mxmu :: Double
mxmu = 0.999

-- Plot --

pltrng :: [Double]
pltrng = range mnx mxx 1000

finalLineFun :: Mean #> Natural # NeuralNetwork' -> [Double]
finalLineFun mlp = S.head . coordinates <$> mlp >$>* (xs ++ pltrng)


-- CSV --

data CrossEntropyDescent = CrossEntropyDescent
    { sgdConditionalCrossEntropy :: Double
    , momentumConditionalCrossEntropy :: Double
    , adamConditionalCrossEntropy :: Double }
    deriving (Generic, Show)

instance FromNamedRecord CrossEntropyDescent
instance ToNamedRecord CrossEntropyDescent
instance DefaultOrdered CrossEntropyDescent
instance NFData CrossEntropyDescent

data RegressionLines = RegressionLines
    { intput :: Double
    , output :: Maybe Double
    , sgdMeanOutput :: Double
    , momentumMeanOutput :: Double
    , adamMeanOutput :: Double }
    deriving (Generic, Show)

instance FromNamedRecord RegressionLines
instance ToNamedRecord RegressionLines
instance DefaultOrdered RegressionLines
instance NFData RegressionLines


--- Main ---


main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let cost = stochasticConditionalCrossEntropy xs ys

    let !mxs = sufficientStatistic <$> xs
        !mys = sufficientStatistic <$> ys

    let backprop :: Mean #> Natural # NeuralNetwork'
                 -> CotangentPair (Mean #> Natural) NeuralNetwork'
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential0 mxs mys p

        sgdmlps0 mlp = take nepchs $ vanillaGradientSequence backprop eps Classic mlp
        mtmmlps0 mlp = take nepchs
            $ vanillaGradientSequence backprop eps (defaultMomentumPursuit mxmu) mlp
        admmlps0 mlp = take nepchs
            $ vanillaGradientSequence backprop eps defaultAdamPursuit mlp

    goalCriterionMain expnm
       [ C.bench "sgd" $ C.nf sgdmlps0 mlp0
       , C.bench "momentum" $ C.nf mtmmlps0 mlp0
       , C.bench "adam" $ C.nf admmlps0 mlp0 ]

    let sgdmlps = sgdmlps0 mlp0
        mtmmlps = mtmmlps0 mlp0
        admmlps = admmlps0 mlp0

    let sgdln = finalLineFun $ last sgdmlps
        mtmln = finalLineFun $ last mtmmlps
        admln = finalLineFun $ last admmlps
        smps = (Just <$> ys) ++ repeat Nothing

    let rgcsv = zipWith5 RegressionLines (xs ++ pltrng) smps sgdln mtmln admln

    goalWriteAnalysis "benchmarks" "regression" "regression-lines" Nothing rgcsv

    let sgdcst = cost <$> sgdmlps
        mtmcst = cost <$> mtmmlps
        admcst = cost <$> admmlps

    let cstcsv = zipWith3 CrossEntropyDescent sgdcst mtmcst admcst

    goalWriteAnalysis "benchmarks" "regression" "cross-entropy-descent" Nothing cstcsv


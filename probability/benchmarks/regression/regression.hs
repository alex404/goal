{-# LANGUAGE DeriveGeneric,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Data.List
import Paths_goal_probability

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
        '[ '(Affine Tensor,R 50 Bernoulli)] (Affine Tensor) (MeanNormal (1/1)) (MeanNormal (1/1))

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

finalLineFun :: Mean #> Natural # NeuralNetwork' -> [Double]
finalLineFun mlp = S.head . coordinates <$> mlp >$>* pltrng


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

data RegressionSamples = RegressionSamples
    { xSample :: Double
    , ySample :: Double }
    deriving (Generic, Show)

instance FromNamedRecord RegressionSamples
instance ToNamedRecord RegressionSamples
instance DefaultOrdered RegressionSamples
instance NFData RegressionSamples

data RegressionLines = RegressionLines
    { input :: Double
    , sgdMeanOutput :: Double
    , momentumMeanOutput :: Double
    , adamMeanOutput :: Double }
    deriving (Generic, Show)

instance FromNamedRecord RegressionLines
instance ToNamedRecord RegressionLines
instance DefaultOrdered RegressionLines
instance NFData RegressionLines

expmnt :: Experiment
expmnt = Experiment "benchmarks" "regression"

--- Main ---


main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let xys = zip ys xs

    let cost :: Mean #> Natural # NeuralNetwork' -> Double
        cost = conditionalLogLikelihood xys

    let backprop :: Mean #> Natural # NeuralNetwork' -> Mean #> Natural #* NeuralNetwork'
        backprop = conditionalLogLikelihoodDifferential xys

        sgdmlps0 mlp = take nepchs $ mlp0 : vanillaGradientSequence backprop eps Classic mlp
        mtmmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps (defaultMomentumPursuit mxmu) mlp
        admmlps0 mlp = take nepchs
            $ mlp0 : vanillaGradientSequence backprop eps defaultAdamPursuit mlp

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

    let smpcsv = zipWith RegressionSamples xs ys

    let rgcsv = zipWith4 RegressionLines pltrng sgdln mtmln admln

    let sgdcst = cost <$> sgdmlps
        mtmcst = cost <$> mtmmlps
        admcst = cost <$> admmlps

    let cstcsv = zipWith3 CrossEntropyDescent sgdcst mtmcst admcst

    goalExportNamed True expmnt Nothing smpcsv
    goalExportNamed False expmnt Nothing rgcsv
    goalExportNamed False expmnt Nothing cstcsv

    rglngpi <- getDataFileName "benchmarks/regression/regression-lines.gpi"
    cedsgpi <- getDataFileName "benchmarks/regression/cross-entropy-descent.gpi"

    runGnuplot expmnt Nothing defaultGnuplotOptions rglngpi
    runGnuplot expmnt Nothing defaultGnuplotOptions cedsgpi

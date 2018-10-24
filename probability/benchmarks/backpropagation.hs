{-# LANGUAGE BangPatterns,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified Criterion.Main as C

--- Globals ---


expnm :: String
expnm = "backpropagation"

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

--type Layer1 = MeanNormal (1/1)
--type Layer2 = R 1000 Bernoulli
--type Layer3 = R 1000 Bernoulli
--type Layer4 = MeanNormal (1/1)
--
--type NeuralNetwork' =
--    NeuralNetworkLayer
--        (Affine Tensor)
--        (NeuralNetworkLayer (Affine Tensor) (Affine Tensor) Layer2)
--        Layer3 Layer4 Layer1

--type NeuralNetwork' = NeuralNetworkLayer (Affine Tensor) (Affine Tensor) Layer2 Layer4 Layer1

type NeuralNetwork' = NeuralNetwork
        [Affine Tensor, Affine Tensor, Affine Tensor]
        [MeanNormal (1/1), R 1000 Bernoulli, R 1000 Bernoulli, (MeanNormal (1/1))]


-- Training --

nepchs :: Int
nepchs = 1

eps :: Double
eps = -0.0001

-- Layout --

main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let !mxs = sufficientStatistic <$> xs
        !mys = sufficientStatistic <$> ys

    let cost :: Mean #> Natural # NeuralNetwork' -> Double
        cost = stochasticConditionalCrossEntropy xs ys

    let backprop :: Point (Mean #> Natural) NeuralNetwork' -> CotangentPair (Mean #> Natural) NeuralNetwork'
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential0 mxs mys p

        admmlps0 mlp = take nepchs $ vanillaGradientSequence backprop eps defaultAdamPursuit mlp

    let mlp = last $!! admmlps0 mlp0

    goalCriterionMain expnm
       [ C.bench "application" $ C.nf cost mlp
       , C.bench "backpropagation" $ C.nf backprop mlp ]

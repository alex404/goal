{-# LANGUAGE BangPatterns,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

-- Qualified --

import qualified Criterion.Main as C

--- Globals ---


-- Data --

f :: Double -> Double
f x = exp . sin $ 2 * x

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: B.Vector 200 Double
xs = B.range mnx mxx

fp :: Source # Normal
fp = Point $ S.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ S.doubleton 0 0.0001

type Layer1 = MeanNormal (1/1)
type Layer2 = R 1000 Bernoulli
type Layer3 = R 1000 Bernoulli
type Layer4 = MeanNormal (1/1)

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

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8

-- Layout --

main :: IO ()
main = do

    ys <- realize $ B.mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let !mxs = joinBoxedReplicated $ sufficientStatistic <$> xs
        !mys = joinBoxedReplicated $ sufficientStatistic <$> ys

    let cost :: Mean ~> Natural # NeuralNetwork' -> Double
        cost = stochasticConditionalCrossEntropy xs ys

    let backprop :: Point (Mean ~> Natural) NeuralNetwork' -> CotangentPair (Mean ~> Natural) NeuralNetwork'
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential0 mxs mys p

        admmlps0 mlp = take nepchs $ vanillaAdamSequence eps b1 b2 rg backprop mlp

    let mlp = last $!! admmlps0 mlp0

--    print $ cost mlp
--    print . S.sum . coordinates $ backprop mlp
    C.defaultMain
       [ C.bench "application" $ C.nf cost mlp
       , C.bench "backpropagation" $ C.nf backprop mlp ]

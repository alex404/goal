{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

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

xs :: B.Vector 20 Double
xs = B.range mnx mxx

fp :: Source # Normal
fp = Point $ B.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ B.doubleton 0 0.1

type NN = MeanNormal (1/1) <*< R 100 Bernoulli <* MeanNormal (1/1)

-- Training --

nepchs :: Int
nepchs = 10

eps :: Double
eps = -0.05

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8

-- Layout --

main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let cost :: RealFloat x => Point (Mean ~> Natural) NN x -> x
        cost = stochasticConditionalCrossEntropy xs ys

    let cost2 :: Mean ~> Natural # NN -> Double
        cost2 foo = average . B.zipWith stochasticCrossEntropy (B.singleton <$> ys) $ foo >>$>* xs

    let backprop :: RealFloat x => Point (Mean ~> Natural) NN x -> CotangentVector (Mean ~> Natural) NN x
        backprop = differential (stochasticConditionalCrossEntropy xs ys)

    let backprop2 :: Point (Mean ~> Natural) NN Double -> CotangentVector (Mean ~> Natural) NN Double
        backprop2 = backpropagation xs ys

        admmlps0 mlp = take nepchs $ vanillaAdamSequence eps b1 b2 rg cost mlp

    let mlp = last $!! admmlps0 mlp0

    C.defaultMain
       [ C.bench "application" $ C.nf cost mlp
       , C.bench "application2" $ C.nf cost2 mlp
       , C.bench "backpropagation" $ C.nf backprop mlp
       , C.bench "backpropagation2" $ C.nf backprop2 mlp ]


    putStrLn "Mean squared error between backprop gradients:"
    print . average $ (^(2 :: Int)) <$> backprop mlp <-> backprop2 mlp


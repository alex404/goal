{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds,Arrows #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import Goal.Datasets.MNIST


-- Globals ---


-- Network --

ip :: Source # Normal
ip = Point $ doubleton 0 0.001


type MLP = Categorical Int 10 <*< Replicated 100 Bernoulli <* Replicated MNISTSize Poisson

-- Data --


-- Training --

type NBatch = 20


nepchs,tbtch :: Int
nepchs = 2
tbtch = 100

eps :: Double
eps = -0.0001

-- Momentum
mxmu :: Double
mxmu = 0.999

mu :: Int -> Double
mu = defaultMomentumSchedule mxmu

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8


-- Functions --

classifications :: KnownNat n => Vector n (Vector MNISTSize Int) -> Mean ~> Natural # MLP -> Vector n Int
classifications xs mlp =
    maxIndexV <$> classifications0 mlp xs

classifications0 :: KnownNat n => Mean ~> Natural # MLP -> Vector n (Vector MNISTSize Int) -> Vector n (Vector 10 Double)
classifications0 mlp xs =
    zipWithV fmap (density <$> mlp >$>* xs) (replicateV $ generateV id)

classification
    :: KnownNat n
    => Vector n (Vector MNISTSize Int,Int)
    -> Mean ~> Natural # MLP
    -> Double
classification vxys mlp =
    let (xs,ys) = unzipV vxys
        classy i j = if i == j then 1 else 0
     in (/ fromIntegral (length vxys)) . sum . zipWithV classy ys $ classifications xs mlp

-- Main --


main :: IO ()
main = do

    txys <- mnistTrainingData

    mlp0 <- realize $ initialize ip

    let tstrm = breakStream txys

    let trncrc :: Circuit (Vector NBatch (Vector MNISTSize Int,Int)) (Mean ~> Natural # MLP)
        trncrc = accumulateCircuit0 mlp0 $ proc (xys,mlp) -> do
            let dmlp = differential (stochasticConditionalCrossEntropy xys) mlp
                dmlp' = joinTangentPair mlp (breakChart dmlp)
            adamAscent eps b1 b2 rg -< dmlp'
            --gradientAscent eps -< dmlp'

    vxys0 <- mnistTestData

    let vxys' :: Vector 10000 (Vector MNISTSize Int,Int)
        vxys' = strongVectorFromList vxys0

    let vxys :: Vector 10 (Vector MNISTSize Int,Int)
        vxys = fst $ splitV vxys'

    let (vxs,vys) = unzipV vxys

    let ces = take nepchs . takeEvery tbtch . stream tstrm $ trncrc >>> arr (classifications vxs)

    --let ces = stochasticConditionalCrossEntropy vxys <$> mlps
    sequence_ $ print . zipV vys <$> ces


{-
    let celyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..] ces ]

    goalRenderableToSVG "mnist" "cross-entropy-descent" 1000 500 $ toRenderable celyt
    -}
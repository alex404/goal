{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds,Arrows #-}

--- Imports ---


-- Goal --

import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

-- Globals ---


-- Network --

ip :: Source # Normal
ip = Point $ S.doubleton 0 0.001

--data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) (ih :: Nat) (oh :: Nat) om im
--type NKernels = 10
--type ConvolutionalLayer = Affine (Convolutional 3 Height Width)

type MLP = NeuralNetwork
    [Affine Tensor, Affine Tensor, Affine Tensor]
    [Categorical Int 10, Replicated 128 Bernoulli, Replicated 256 Bernoulli, Replicated Length (MeanNormal (1/1)) ]
--type MLP = NeuralNetwork Categorical Int 10 <*< Replicated 100 Bernoulli <*< ConvolutionalLayer


-- Data --


-- Training --

nbtch :: Int
nbtch = 10

nepchs,epchn :: Int
nepchs = 10
epchn = 100

eps :: Double
eps = -0.0005


-- Functions --

classifications :: [B.Vector Length Double] -> Mean #> Natural # MLP -> [Int]
classifications xs mlp = fromIntegral . B.maxIndex <$> classifications0 mlp xs

classifications0 :: Mean #> Natural # MLP -> [B.Vector Length Double] -> [B.Vector 10 Double]
classifications0 mlp zs = [ density nx <$> B.generate finiteInt | nx <- mlp >$>* zs ]

--l2norm :: Point (Mean #> Natural) MLP -> Dou
--l2norm mlp = sqrt . sum $ (^(2:: Int)) <$> mlp

accuracy
    :: [(B.Vector Length Double,Int)]
    -> Mean #> Natural # MLP
    -> (Double,Double)
accuracy vxys mlp =
    let (xs,ys) = unzip vxys
        classy i j = if i == j then 1 else 0
     in (average . zipWith classy ys $ classifications xs mlp, S.maximum . S.map abs $ coordinates mlp)

forcePoint :: Point c m -> Point c m
forcePoint (Point xs) = Point $ S.force xs

-- Main --


main :: IO ()
main = do

    txys <- mnistTrainingData

    mlp0 <- realize $ initialize ip

    vxys0 <- mnistTestData

    let trngn = generator (breakEvery nbtch $ cycle txys)


    let vxys :: [(B.Vector Length Double,Int)]
        vxys = take 1000 vxys0

    let trncrc :: [(B.Vector Length Double,Int)] >>> (Double,Double)
        trncrc = accumulateCircuit mlp0 $ proc (xys,mlp) -> do
            let (xs,ys) = unzip xys
                dmlp1 = stochasticConditionalCrossEntropyDifferential xs ys $ forcePoint mlp
                --dmlp1 = differential (stochasticConditionalCrossEntropy xs ys) mlp
                --dmlp2 = differential l2norm mlp
                --dmlp = convexCombination 0.99 dmlp1 dmlp2
                --dmlp' = joinTangentPair mlp (breakChart dmlp)
            mlp' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair mlp $ breakPoint (forcePoint dmlp1)
            returnA -< (accuracy vxys mlp',mlp')
            --momentumAscent eps mu -< dmlp'
            --gradientAscent eps -< dmlp'

    let ces = take nepchs . takeEvery epchn . streamChain $ trncrc <<< trngn

    --let ces = stochasticConditionalCrossEntropy vxys <$> mlps
    sequence_ $ print <$> ces
    --print $ mlp0 >.>* (fst $ B.head vxys)


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

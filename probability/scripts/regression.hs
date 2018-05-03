{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

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

xs :: B.Vector 20 Double
xs = B.range mnx mxx

fp :: Source # Normal
fp = Point $ S.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ S.doubleton 0 0.1

type Layer1 = MeanNormal (1/1)
type Layer2 = R 50 Bernoulli
type Layer3 = MeanNormal (1/1)
type NeuralNetwork = (Layer3 <*< Layer2) (Affine Product) Layer1

type (m <*< n) g o = NeuralNetworkLayer (Affine Product) g n m o
infixr 3 <*<

-- Training --

nepchs :: Int
nepchs = 1000

eps :: Double
eps = -0.05

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

-- Plot --

pltrng :: B.Vector 1000 Double
pltrng = B.range mnx mxx

-- Layout --

main :: IO ()
main = do

    ys <- realize $ B.mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let cost = stochasticConditionalCrossEntropy xs ys

    let !mxs = joinBoxedReplicated $ sufficientStatistic <$> xs
        !mys = joinBoxedReplicated $ sufficientStatistic <$> ys

    let backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential0 mxs mys p

        sgdmlps0 mlp = take nepchs $ vanillaGradientSequence eps backprop mlp
        mtmmlps0 mlp = take nepchs $ vanillaMomentumSequence eps mu backprop mlp
        admmlps0 mlp = take nepchs $ vanillaAdamSequence eps b1 b2 rg backprop mlp

    C.defaultMain
       [ C.bench "sgd" $ C.nf sgdmlps0 mlp0
       , C.bench "momentum" $ C.nf mtmmlps0 mlp0
       , C.bench "adam" $ C.nf admmlps0 mlp0 ]

    let sgdmlps = sgdmlps0 mlp0
        mtmmlps = mtmmlps0 mlp0
        admmlps = admmlps0 mlp0

    let finalLineFun :: Mean ~> Natural # NeuralNetwork -> [(Double,Double)]
        finalLineFun mlp =
            let ys' = S.map (S.head . coordinates) . splitReplicated $ mlp >$>* pltrng
             in zip (B.toList pltrng) (S.toList ys')

        lyt1 = execEC $ do

            goalLayout

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip (B.toList pltrng) (B.toList $ f <$> pltrng)]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [finalLineFun $ last sgdmlps]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [finalLineFun $ last mtmmlps]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque green)
                plot_lines_values .= [finalLineFun $ last admmlps]

            plot . liftEC $ do

                plot_points_style .=  filledCircles 5 (opaque black)
                plot_points_values .= toList (zip (B.toList xs) (B.toList ys))

        lyt2 = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..] $ cost <$> sgdmlps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip [0..] $ cost <$> mtmmlps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque green)
                plot_lines_values .= [ zip [0..] $ cost <$> admmlps ]

    goalRenderableToSVG "probability/backpropagation" "regression" 500 200 $ toRenderable lyt1
    goalRenderableToSVG "probability/backpropagation" "descent" 500 200 $ toRenderable lyt2

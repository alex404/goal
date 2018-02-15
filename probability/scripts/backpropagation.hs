{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

-- Qualified --

import qualified Criterion.Main as C

--- Globals ---


-- Data --

f :: Double -> Double
f x = exp . sin $ 2 * x

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: Vector 20 Double
xs = rangeV mnx mxx

fp :: Source # Normal
fp = Point $ doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ doubleton 0 0.1

type NeuralNetwork = MeanNormal (1/1) <*< R 50 Bernoulli <* MeanNormal (1/1)

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

pltrng :: Vector 1000 Double
pltrng = rangeV mnx mxx

-- Layout --

main :: IO ()
main = do

    ys <- realize $ mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let cost :: RealFloat x => Point (Mean ~> Natural) NeuralNetwork x -> x
        cost = stochasticConditionalCrossEntropy (zipV xs ys)

        sgdmlps0 mlp = take nepchs $ vanillaGradientSequence eps cost mlp
        mtmmlps0 mlp = take nepchs $ vanillaMomentumSequence eps mu cost mlp
        admmlps0 mlp = take nepchs $ vanillaAdamSequence eps b1 b2 rg cost mlp

    C.defaultMain
       [ C.bench "sgd" $ C.nf sgdmlps0 mlp0
       , C.bench "momentum" $ C.nf mtmmlps0 mlp0
       , C.bench "adam" $ C.nf admmlps0 mlp0 ]

    let sgdmlps = sgdmlps0 mlp0
        mtmmlps = mtmmlps0 mlp0
        admmlps = admmlps0 mlp0

    let finalLineFun :: Mean ~> Natural # NeuralNetwork -> [(Double,Double)]
        finalLineFun mlp = toList . zipV pltrng $ headV . coordinates <$> mlp >$>* pltrng

        lyt1 = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "x"
            layout_y_axis . laxis_title .= "y"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [toList $ zipV pltrng (f <$> pltrng)]

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
                plot_points_values .= toList (zipV xs ys)

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

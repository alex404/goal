{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified Criterion.Main as C

--- Globals ---


-- Data --

f :: S.Vector 1 Double -> S.Vector 1 Double
f x = exp . sin $ 2 * x

mnx,mxx :: Double
mnx = -3
mxx = 3

xs :: S.Vector 20 (S.Vector 1 Double)
xs = S.map S.singleton $ S.range mnx mxx

fp :: Source # Normal
fp = Point $ S.doubleton 0 0.1

-- Neural Network --

cp :: Source # Normal
cp = Point $ S.doubleton 0 0.1

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

pltrng :: S.Vector 1000 (S.Vector 1 Double)
pltrng = S.map S.singleton $ S.range mnx mxx

toLine :: S.Vector k (S.Vector 1 Double) -> [Double]
toLine xs' = S.head <$> S.toList xs'

-- Layout --

main :: IO ()
main = do

    ys <- realize $ S.mapM (noisyFunction fp f) xs

    mlp0 <- realize $ initialize cp

    let cost = stochasticConditionalCrossEntropy xs ys

    let backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs ys p

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
            let ys' = S.map coordinates . splitReplicated $ mlp >$>* pltrng
             in zip (toLine pltrng) (toLine ys')

        lyt1 = execEC $ do

            goalLayout

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip (toLine pltrng) (toLine $ S.map f pltrng)]

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
                plot_points_values .= toList (zip (toLine xs) (toLine ys))

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

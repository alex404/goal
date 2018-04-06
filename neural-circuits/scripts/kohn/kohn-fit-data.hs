{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.NeuralCircuits

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import Data.List

--- Globals ---

dr,flnm,sbdr,preflnm,pstflnm :: String
dr = "adaptation/small40"
flnm = "112r36/"
sbdr = "neural-circuits/kohn-data/" ++ flnm
preflnm = "prestms"
pstflnm = "pststms"

type NStimuli = 8
type NNeurons = 13

pltsmps :: B.Vector 100 Double
pltsmps = B.range 0 (2*pi)

unsafeFromListS :: (KnownNat k, Storable x) => [x] -> S.Vector k x
unsafeFromListS = fromJust . S.fromList

unsafeFromListB :: KnownNat k => [x] -> B.Vector k x
unsafeFromListB = fromJust . B.fromList

sps :: S.Vector NNeurons (Source # VonMises)
sps = unsafeFromListS [ Point $ S.doubleton mu 1 | mu <- tail $ range 0 (2*pi) (1 + natValInt (Proxy :: Proxy NNeurons)) ]

ppc0 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
ppc0 = vonMisesPopulationEncoder sps 1

sinusoid :: Double -> S.Vector 3 Double
sinusoid x =
     fromJust $ S.fromList [cos x,sin x,1]

sumOfTuningCurves
    :: (KnownNat k, KnownNat j)
    => B.Vector k (B.Vector j (Double, Double))
    -> B.Vector j (Double, Double)
sumOfTuningCurves tcs =
    let tcs' = B.unzip <$> tcs
        xs' = fst . head $ B.toList tcs'
        cs' = sum $ snd <$> tcs'
     in B.zip xs' cs'

-- Training --

nepchs :: Int
nepchs = 10000

eps :: Double
eps = -0.005

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8


--- Functions ---


mapToSample
    :: M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> (B.Vector NStimuli (Sample VonMises), B.Vector NStimuli (Sample (R NNeurons Poisson)))
mapToSample stmmp =
    let xs = unsafeFromListB $ M.keys stmmp
        ys = unsafeFromListB
            [ unsafeFromListB . map snd . sort . M.toList $ round . (/(100 :: Double)) . genericLength <$> nmp | nmp <- M.elems stmmp ]
     in (xs, ys)

tclyt :: (Mean ~> Natural # R NNeurons Poisson <* VonMises)
      -> Double
      -> LayoutLR Double Double Double
tclyt lkl adpt = execEC $ do

    let tcs = tuningCurves pltsmps lkl
        stcs = sumOfTuningCurves tcs
        cs = snd <$> stcs
        alphs = linearLeastSquares (S.map sinusoid $ G.convert pltsmps) $ G.convert cs
        sinsmps = S.dotProduct alphs . sinusoid <$> pltsmps
        ssres = sum $ (^(2::Int)) <$> (cs - sinsmps)
        sstot = sum $ (^(2::Int)) . subtract (average cs) <$> cs
        rmse = sqrt $ ssres / fromIntegral (length cs)
        r2 = 1 - (ssres/sstot)
        amp =
            let x1:x2:_ = S.toList alphs
             in sqrt $ x1^2 + x2^2

    layoutlr_title .= ("Amplitude: " ++ show amp ++ ", RMSE: " ++ show rmse ++ ", r^2: " ++ show r2)
    goalLayoutLR
    radiansAbscissaLR

    layoutlr_x_axis . laxis_title .= "Stimulus"
    layoutlr_left_axis . laxis_title .= "Firing Rate"
    layoutlr_right_axis . laxis_title .= "Total Rate"
    layoutlr_left_axis . laxis_generate .= scaledAxis def (0,50)
    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,200)

    let adptpi0 = 2*pi*adpt/360
        adptpi =
            2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)
    plotRight . return $ vlinePlot "adaptor" (solidLine 4 $ opaque blue) adptpi

    plotLeft . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tcs)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList stcs ]

    plotRight . liftEC $ do
        plot_lines_style .= dashedLine 3 [4,4] (opaque black)
        plot_lines_values .= [ B.toList $ B.zip pltsmps sinsmps ]

--- Main ---


main :: IO ()
main = do

    adpt <- getAdaptor dr flnm

    (prestms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr preflnm
    (pststms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr pstflnm

    putStr "Number of Neurons: "
    print . length. M.keys . head $ toList prestms

    let (xs,prens) = mapToSample prestms
        pstns = snd $ mapToSample pststms

    let cost = stochasticConditionalCrossEntropy xs

    let backprop ns p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs ns p

        preppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop prens) ppc0
        pstppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop pstns) ppc0
        preppc = last preppcs
        pstppc = last pstppcs

        nlllyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_title .= "pre"
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..] $ cost prens <$> preppcs ]

            plot . liftEC $ do

                plot_lines_title .= "post"
                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip [0..] $ cost pstns <$> pstppcs ]

    goalRenderableToSVG (sbdr ++ "/fit") "initial-population" 1200 800 . toRenderable $ tclyt ppc0 adpt
    goalRenderableToSVG (sbdr ++ "/fit") "pre-population" 1200 800 . toRenderable $ tclyt preppc adpt
    goalRenderableToSVG (sbdr ++ "/fit") "post-population" 1200 800 . toRenderable $ tclyt pstppc adpt
    goalRenderableToSVG (sbdr ++ "/fit") "negative-log-likelihood" 1200 800 . toRenderable $ nlllyt

-- Data --

--f :: Double -> Double
--f x = exp . sin $ 2 * x
--
--mnx,mxx :: Double
--mnx = -3
--mxx = 3
--
--xs :: B.Vector 20 Double
--xs = B.range mnx mxx
--
--fp :: Source # Normal
--fp = Point $ S.doubleton 0 0.1
--
---- Neural Network --
--
--cp :: Source # Normal
--cp = Point $ S.doubleton 0 0.1
--
--type NeuralNetwork = MeanNormal (1/1) <*< R 50 Bernoulli <* MeanNormal (1/1)
--
---- Training --
--
--nepchs :: Int
--nepchs = 1000
--
--eps :: Double
--eps = -0.05
--
---- Momentum
--mxmu :: Double
--mxmu = 0.999
--
--mu :: Int -> Double
--mu = defaultMomentumSchedule mxmu
--
---- Adam
--b1,b2,rg :: Double
--b1 = 0.9
--b2 = 0.999
--rg = 1e-8
--
---- Plot --
--
--pltrng :: B.Vector 1000 Double
--pltrng = B.range mnx mxx
--
---- Layout --
--
--main :: IO ()
--main = do
--
--    ys <- realize $ mapM (noisyFunction fp f) xs
--
--    mlp0 <- realize $ initialize cp
--
--    let cost = stochasticConditionalCrossEntropy xs ys
--
--    let backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs ys p
--
--        sgdmlps0 mlp = take nepchs $ vanillaGradientSequence eps backprop mlp
--        mtmmlps0 mlp = take nepchs $ vanillaMomentumSequence eps mu backprop mlp
--        admmlps0 mlp = take nepchs $ vanillaAdamSequence eps b1 b2 rg backprop mlp
--
--    let sgdmlps = sgdmlps0 mlp0
--        mtmmlps = mtmmlps0 mlp0
--        admmlps = admmlps0 mlp0
--
--    let finalLineFun :: Mean ~> Natural # NeuralNetwork -> [(Double,Double)]
--        finalLineFun mlp = toList . B.zip pltrng $ S.head . coordinates <$> G.convert (mlp >$>* pltrng)
--
--        lyt1 = execEC $ do
--
--            goalLayout
--            layout_x_axis . laxis_title .= "x"
--            layout_y_axis . laxis_title .= "y"
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque black)
--                plot_lines_values .= [toList $ B.zip pltrng (f <$> pltrng)]
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque red)
--                plot_lines_values .= [finalLineFun $ last sgdmlps]
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque blue)
--                plot_lines_values .= [finalLineFun $ last mtmmlps]
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque green)
--                plot_lines_values .= [finalLineFun $ last admmlps]
--
--            plot . liftEC $ do
--
--                plot_points_style .=  filledCircles 5 (opaque black)
--                plot_points_values .= toList (B.zip xs ys)
--
--        lyt2 = execEC $ do
--
--            goalLayout
--            layout_x_axis . laxis_title .= "Epochs"
--            layout_y_axis . laxis_title .= "Negative Log-Likelihood"
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque red)
--                plot_lines_values .= [ zip [(0 :: Int)..] $ cost <$> sgdmlps ]
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque blue)
--                plot_lines_values .= [ zip [0..] $ cost <$> mtmmlps ]
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque green)
--                plot_lines_values .= [ zip [0..] $ cost <$> admmlps ]
--
--    goalRenderableToSVG "probability/backpropagation" "regression" 500 200 $ toRenderable lyt1
--    goalRenderableToSVG "probability/backpropagation" "descent" 500 200 $ toRenderable lyt2

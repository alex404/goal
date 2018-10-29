{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies,FlexibleContexts #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Data.Char


--- Globals ---

nsmps :: Int
nsmps = 20

lineFun1 :: AbsolutelyContinuous Source m
         => (SamplePoint m -> Double) -> Sample m -> Source # m -> [(Double,Double)]
lineFun1 toDouble rng p = zip (toDouble <$> rng) $ density p <$> rng

lineFun2 :: AbsolutelyContinuous Natural m
         => (SamplePoint m -> Double) -> Sample m -> Natural # m -> [(Double,Double)]
lineFun2 toDouble rng p = zip (toDouble <$> rng) $ density p <$> rng

-- Bernoulli --

ttlB :: String
ttlB = "Bernoulli"

mnB,mxB :: Double
(mnB,mxB) = (0,1)

bnsB :: Int
bnsB = 2

truB :: Source # Bernoulli
truB = Point $ S.singleton 0.7

toDoubleB :: Bool -> Double
toDoubleB b = if b then 1 else 0

rngB :: [Bool]
rngB = pointSampleSpace truB

-- Categorical --

ttlC :: String
ttlC = "Categorical"

mnC,mxC :: Double
(mnC,mxC) = (0,4)

bnsC :: Int
bnsC = 5

truC :: Source # Categorical Int 5
truC = Point . fromJust $ S.fromList [0.1,0.4,0.1,0.2]

toDoubleC :: Int -> Double
toDoubleC = fromIntegral

rngC :: [Int]
rngC = pointSampleSpace truC

-- Poisson --

ttlP :: String
ttlP = "Poisson"

mnP,mxP :: Double
(mnP,mxP) = (0,20)

bnsP :: Int
bnsP = 20

truP :: Source # Poisson
truP = Point $ S.singleton 5

toDoubleP :: Int -> Double
toDoubleP = fromIntegral

rngP :: [Int]
rngP = [0..20]

-- Normal --

ttlN :: String
ttlN = "Normal"

mnN,mxN :: Double
(mnN,mxN) = (-3,7)

bnsN :: Int
bnsN = 20

truN :: Source # Normal
truN = Point $ S.doubleton 2 0.7

toDoubleN :: Double -> Double
toDoubleN = id

rngN :: [Double]
rngN = range (-3) 7 100

-- Layout --

generateLayout
    :: ( Transition Source Mean m, Transition Source Natural m, ClosedFormExponentialFamily m
       , MaximumLikelihood Source m, AbsolutelyContinuous Source m, Generative Source m
       , AbsolutelyContinuous Natural m )
       => String
       -> Int
       -> Double
       -> Double
       -> (SamplePoint m -> Double)
       -> Sample m
       -> Source # m
       -> IO (LayoutLR Double Int Double)
generateLayout ttl nb mn mx toDouble rng p = do

    smps <- realize $ sample nsmps p

    let mle1 = mle smps

    return . execEC $ do

        goalLayoutLR

        histogramLayoutLR nb mn mx

        layoutlr_title .= (ttl ++ "; KLD: " ++ take 5 (showFFloat (Just 3) (transition2 relativeEntropy mle1 p) ""))
        layoutlr_left_axis . laxis_title .= "Sample Count"
        layoutlr_right_axis . laxis_title .= "Probability Mass"
        layoutlr_x_axis . laxis_title .= "Value"

        plotLeft . fmap plotBars . liftEC $ do
            void $ histogramPlot nb mn mx [toDouble <$> smps]
            plot_bars_titles .= ["Samples"]
            plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

        plotRight . liftEC $ do
            plot_lines_style .= solidLine 3 (opaque black)
            plot_lines_title .= "True"
            plot_lines_values .= [lineFun1 toDouble rng p]

        plotRight . liftEC $ do
            plot_lines_style .= solidLine 3 (opaque red)
            plot_lines_title .= "MLE"
            plot_lines_values .= [ lineFun1 toDouble rng mle1 ]

        plotRight . liftEC $ do
            plot_lines_style .= dashedLine 3 [8,8] (opaque purple)
            plot_lines_title .= "MLE"
            plot_lines_values .= [ lineFun2 toDouble rng $ transition mle1 ]

main :: IO ()
main = do

    rnblB <- toRenderable <$> generateLayout ttlB bnsB mnB mxB toDoubleB rngB truB
    rnblC <- toRenderable <$> generateLayout ttlC bnsC mnC mxC toDoubleC rngC truC
    rnblP <- toRenderable <$> generateLayout ttlP bnsP mnP mxP toDoubleP rngP truP
    rnblN <- toRenderable <$> generateLayout ttlN bnsN mnN mxN toDoubleN rngN truN

    goalRenderableToSVG "probability/univariate" (toLower <$> ttlB) 600 300 rnblB
    goalRenderableToSVG "probability/univariate" (toLower <$> ttlC) 800 400 rnblC
    goalRenderableToSVG "probability/univariate" (toLower <$> ttlP) 800 400 rnblP
    goalRenderableToSVG "probability/univariate" (toLower <$> ttlN) 800 400 rnblN



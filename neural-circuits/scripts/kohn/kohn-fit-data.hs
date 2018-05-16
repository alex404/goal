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
dr = "adaptation/big40"
flnm = "112r32/"
sbdr = "neural-circuits/kohn-data/" ++ flnm
preflnm = "prestms"
pstflnm = "pststms"
presflnm = "prestrm"
pstsflnm = "pststrm"

type NStimuli = 8
type NNeurons = 126
type PreTrials = 400
type PostTrials = 320

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
    :: Double
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> (B.Vector NStimuli (SamplePoint VonMises), B.Vector NStimuli (SamplePoint (R NNeurons Poisson)))
mapToSample nrm stmmp =
    let xs = unsafeFromListB $ M.keys stmmp
        ys = unsafeFromListB
            [ unsafeFromListB . map snd . sort . M.toList $ round . (/nrm) . genericLength <$> nmp | nmp <- M.elems stmmp ]
     in (xs, ys)

streamToRawSum
    :: KnownNat k
    => [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> (B.Vector k Stimulus, B.Vector k Double)
streamToRawSum sstrm =
    let (xs,ns0) = B.unzip $ unsafeFromListB sstrm
     in ( xs, genericLength . M.foldr (++) [] <$> ns0)


rwlyt :: (KnownNat k, 1 <= k) => (B.Vector k Stimulus, B.Vector k Double) -> Layout Double Double
rwlyt (xs,sms) = execEC $ do

    let alphs = linearLeastSquares (S.map sinusoid $ G.convert xs) $ G.convert sms
        --ssres = sum $ (^(2::Int)) <$> (sms - sinsmps)
        --sstot = sum $ (^(2::Int)) . subtract (average sms) <$> sms
        --rmse = sqrt $ ssres / fromIntegral (length sms)
        --r2 = 1 - (ssres/sstot)
        amp =
            let x1:x2:_ = S.toList alphs
             in sqrt $ x1^(2::Int) + x2^(2::Int)
        sinsmps = S.dotProduct alphs . sinusoid <$> pltsmps

    --layout_title .= ("Amplitude: " ++ take 4 (show amp) ++ ", RMSE: " ++ take 4 (show rmse) ++ ", r^2: " ++ take 4 (show r2))
    goalLayout
    radiansAbscissa

    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Total Rate"

    plot . liftEC $ do
        plot_points_style .= hollowCircles 4 2 (opaque black)
        plot_points_values .= toList (B.zip xs sms)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [ B.toList $ B.zip pltsmps sinsmps ]

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
             in sqrt $ x1^(2::Int) + x2^(2::Int)

    layoutlr_title .= ("Amplitude: " ++ take 4 (show amp) ++ ", RMSE: " ++ take 4 (show rmse) ++ ", r^2: " ++ take 4 (show r2))
    goalLayoutLR
    radiansAbscissaLR

    layoutlr_x_axis . laxis_title .= "Stimulus"
    layoutlr_left_axis . laxis_title .= "Firing Rate"
    layoutlr_right_axis . laxis_title .= "Total Rate"
    layoutlr_left_axis . laxis_generate .= scaledAxis def (0,300)
    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,5000)

    let adptpi0 = 2*pi*adpt/360
        adptpi =
            2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)
    plotRight . return $ vlinePlot "" (solidLine 4 $ opaque black) adptpi

    plotLeft . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tcs)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList stcs ]

    plotRight . liftEC $ do
        plot_lines_style .= dashedLine 3 [8,10] (opaque blue)
        plot_lines_values .= [ B.toList $ B.zip pltsmps sinsmps ]


--- Main ---


main :: IO ()
main = do

    adpt <- getAdaptor dr flnm

    (prestms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr preflnm
    (pststms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile sbdr pstflnm
    (prestrm :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile sbdr presflnm
    (pststrm :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile sbdr pstsflnm

    putStr "Number of Neurons: "
    print . length. M.keys . head $ toList prestms
    putStr "Stimuli: "
    print $ M.keys prestms
    putStr "Number of Pre-Trials: "
    print $ length prestrm
    putStr "Number of Post-Trials: "
    print $ length pststrm

    let (xs,prens) = mapToSample 25 prestms
        pstns = snd $ mapToSample 20 pststms

        prerwsm :: (B.Vector PreTrials Stimulus, B.Vector PreTrials Double)
        prerwsm = streamToRawSum prestrm
        pstrwsm :: (B.Vector PostTrials Stimulus, B.Vector PostTrials Double)
        pstrwsm = streamToRawSum pststrm

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

    --goalRenderableToSVG (sbdr ++ "/fit") "initial-population" 400 300 . toRenderable $ tclyt ppc0 adpt
    goalRenderableToSVG (sbdr ++ "/fit") ("pre-population-" ++ (reverse . tail $ reverse flnm)) 1200 600 . toRenderable $ tclyt preppc adpt
    goalRenderableToSVG (sbdr ++ "/fit") ("post-population-" ++ (reverse . tail $ reverse flnm)) 1200 600 . toRenderable $ tclyt pstppc adpt
    goalRenderableToSVG (sbdr ++ "/fit") ("raw-pre-population-" ++ (reverse . tail $ reverse flnm)) 1200 600 . toRenderable $ rwlyt prerwsm
    goalRenderableToSVG (sbdr ++ "/fit") ("raw-post-population-" ++ (reverse . tail $ reverse flnm)) 1200 600 . toRenderable $ rwlyt pstrwsm
    --goalRenderableToSVG (sbdr ++ "/fit") "negative-log-likelihood" 400 300 . toRenderable $ nlllyt
    --goalRenderableToSVG (sbdr ++ "/fit") "negative-log-likelihood" 400 300 . toRenderable $ nlllyt

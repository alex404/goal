{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import KohnLab

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M

--- Globals ---


-- Training --

nepchs :: Int
nepchs = 2000

eps :: Double
eps = -0.005

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8

ftsmps :: S.Vector 100 Double
ftsmps = S.init $ S.range 0 (2*pi)

bftsmps :: B.Vector 100 Double
bftsmps = G.convert ftsmps

--- Functions ---


unsafeFromListB :: KnownNat k => [x] -> B.Vector k x
unsafeFromListB = fromJust . B.fromList

vonMisesFits :: KnownNat nn
      => (Mean ~> Natural # R nn Poisson <* VonMises)
      -> LayoutLR Double Double Double
vonMisesFits lkl = execEC $ do

    let (rho0,rprms) = populationCodeRectificationParameters lkl bftsmps
        stcs = sumOfTuningCurves lkl bftsmps

--        r2 = rSquared ftsmps stchts
--        mse = meanSquaredError ftsmps stchts
--        amp = let x1:x2:_ = S.toList alphs
--               in sqrt $ square x1 + square x2
--
--    layoutlr_title .= ("Amplitude: " ++ take 4 (show amp) ++ ", r^2: " ++ show (roundSD 3 r2))

    goalLayoutLR
    radiansAbscissaLR

    layoutlr_x_axis . laxis_title .= "Stimulus"
    layoutlr_left_axis . laxis_title .= "Firing Rate"
    layoutlr_right_axis . laxis_title .= "Total Rate"
--    layoutlr_left_axis . laxis_generate .= scaledAxis def (0,300)
--    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,5000)

    plotRight . return $ vlinePlot "" (solidLine 2 $ opaque black) pi
    plotRight . return $ hlinePlot "" (solidLine 2 $ opaque white) 0

    plotLeft . liftEC $ do
        plot_lines_title .= "tuning curves"
        plot_lines_style .= solidLine 2 (blue `withOpacity` 0.3)
        plot_lines_values .= toList (toList <$> tuningCurves bftsmps lkl)

    plotRight . liftEC $ do
        plot_lines_title .= "tuning curve sum"
        plot_lines_style .= solidLine 2 (opaque black)
        plot_lines_values .= [ zip (S.toList ftsmps) (S.toList stcs) ]

    plotRight . liftEC $ do
        plot_lines_title .= "sum fit"
        plot_lines_style .= dashedLine 3 [20,20] (opaque red)
        plot_lines_values .= [ toList . B.zip bftsmps $ rectificationCurve rho0 rprms bftsmps ]

streamToTrainingSample
    :: (KnownNat nn, KnownNat t)
    => Proxy t
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> (B.Vector t (SamplePoint VonMises), B.Vector t (SamplePoint (R nn Poisson)))
streamToTrainingSample _ stmstrm =
    let xs = unsafeFromListB $ fst <$> stmstrm
        ys = unsafeFromListB $ unsafeFromListB . toList . fmap length . snd <$> stmstrm
     in (xs, ys)


--- Main ---


fitData
    :: forall nn t1 t2
    . (KnownNat nn, KnownNat t1, KnownNat t2, 1 <= nn, 1 <= t1, 1 <= t2)
    => KohnExperiment nn t1 t2
    -> IO ()
fitData kxp = do

    let kpdr = kohnProjectPath kxp

    --(stmttls0 :: M.Map Stimulus (Int,M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kdr "stmttls0"
    --(stmttls1 :: M.Map Stimulus (Int,M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kdr "stmttls1"
    (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm0"
    (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm1"

    let sps :: S.Vector nn (Source # VonMises)
        sps = fromJust $ S.fromList
            [ Point $ S.doubleton mu 1 | mu <- tail $ range 0 (2*pi) (1 + natValInt (Proxy :: Proxy nn)) ]

    let ppc0 :: Mean ~> Natural # R nn Poisson <* VonMises
        ppc0 = vonMisesPopulationEncoder False (Left 1) sps

    let (xs0,ys0) = streamToTrainingSample (Proxy :: Proxy t1) stmstrm0
        (xs1,ys1) = streamToTrainingSample (Proxy :: Proxy t2) stmstrm1

    let backprop xs ys p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs ys p

        preppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop xs0 ys0) ppc0
        pstppcs = take nepchs $ vanillaAdamSequence eps b1 b2 rg (backprop xs1 ys1) ppc0
        preppc = last preppcs
        pstppc = last pstppcs

        nlllyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_title .= "pre"
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..]
                    $ stochasticConditionalCrossEntropy xs0 ys0 <$> preppcs ]

            plot . liftEC $ do

                plot_lines_title .= "post"
                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip [0..]
                    $ stochasticConditionalCrossEntropy xs1 ys1 <$> pstppcs ]

    let prettl = "Pre-Adapt; Dataset: " ++ experiment kxp
        pstttl = "Post-Adapt; Dataset: " ++ experiment kxp

    goalRenderableToPDF kpdr ("fit-pre-population" ++ experiment kxp) 400 200 . toRenderable
        . (layoutlr_title .~ prettl) $ vonMisesFits preppc

    goalRenderableToPDF kpdr ("fit-post-population" ++ experiment kxp) 400 200 . toRenderable
        . (layoutlr_title .~ pstttl) $ vonMisesFits pstppc
    goalRenderableToPDF kpdr "fit-negative-log-likelihood" 400 200 . toRenderable $ nlllyt


main :: IO ()
main = do
    fitData experiment112l44
    fitData experiment112l45
    fitData experiment112r35
    fitData experiment112r36
    fitData experiment105r62
    fitData experiment107l114
    fitData experiment112l16
    fitData experiment112r32
    --fitData big40Pooled
    fitData small40Pooled

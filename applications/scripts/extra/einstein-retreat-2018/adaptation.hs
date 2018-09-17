{-# LANGUAGE TypeOperators,DataKinds,FlexibleContexts,TypeFamilies #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import Goal.Geometry
import Goal.Probability

--- Program ---


-- Globals --

-- Stimulus
x0 :: Double
x0 = pi

sx0 :: Source # VonMises
sx0 = Point $ S.doubleton pi 0

nx0 :: Natural # VonMises
nx0 = transition sx0

-- Tuning Curves
type NNeurons = 6

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

mus :: B.Vector NNeurons Double
mus = B.init $ B.range mnx mxx

prc :: Double
prc = 2

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map (Point . flip S.doubleton prc) (G.convert mus)

lkl0 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl0 = vonMisesPopulationEncoder sps 1

dec :: Mean ~> Natural # VonMises <* R NNeurons Poisson
dec = joinAffine (transition sx0) . transpose . snd $ splitAffine lkl0

rho0 :: Double
rprms0 :: Natural # VonMises
(rho0,rprms0) = populationCodeRectificationParameters lkl0 mus

rprms1 :: Natural # VonMises
rprms1 = toNatural (Point (S.doubleton 0 6) :: Source # VonMises)

lkl1 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl1 = rectifyPopulationCode rho0 rprms1 mus lkl0

-- Plot
pltsmps :: B.Vector 500 Double
pltsmps = B.range mnx mxx


-- Functions --


tclyt
    :: (Mean ~> Natural # R NNeurons Poisson <* VonMises)
    -> Natural # VonMises
    -> Maybe (Natural # VonMises)
    -> Layout Double Double
tclyt lkl rprms mrprms = execEC $ do

    goalLayout
    radiansAbscissa

    let mxplty = 20
        tcs = toList (toList <$> tuningCurves pltsmps lkl)
        tc = tcs !! 3
        tcs' = take 3 tcs ++ drop 4 tcs

    layout_y_axis . laxis_generate .= scaledAxis def (0,mxplty)
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Firing Rate"

    plot . liftEC $ do
        plot_lines_title .= "tuning curves"
        plot_lines_style .= solidLine 4 (opaque blue)
        plot_lines_values .= [tc]

    plot . liftEC $ do
        plot_lines_style .= solidLine 2 (blue `withOpacity` 0.5)
        plot_lines_values .= tcs'

    plot . liftEC $ do
        plot_lines_title .= "t.c. sum"
        plot_lines_style .= solidLine 2 (opaque black)
        plot_lines_values .= [ zip (B.toList pltsmps) . S.toList $ sumOfTuningCurves lkl pltsmps ]

    case mrprms of
      (Just rprms') ->
        plot . liftEC $ do
            plot_lines_title .= "constant fit"
            plot_lines_style .= dashedLine 2 [20,20] (opaque purple)
            plot_lines_values .=
                [ toList . B.zip pltsmps $ rectificationCurve rho0 rprms' pltsmps ]
      Nothing -> return ()

    plot . liftEC $ do
        plot_lines_title .= "sum fit"
        plot_lines_style .= dashedLine 3 [20,20] (opaque red)
        plot_lines_values .=
            [ toList . B.zip pltsmps $ rectificationCurve rho0 rprms pltsmps ]

byslyt
    :: Mean ~> Natural # R NNeurons Poisson <* VonMises
    -> Natural # VonMises
    -> Maybe (Natural # VonMises)
    -> B.Vector NNeurons Int
    -> LayoutLR Double Int Double
byslyt lkl rprms mrprms z = execEC $ do

    goalLayoutLR
    radiansAbscissaLR

    let pstr = fromOneHarmonium . rectifiedBayesRule rprms lkl z $ toOneHarmonium nx0
        uposterior = unnormalizedDensity pstr
        nrm = integrate 1e-2000 uposterior 0 (2*pi)
        posterior x = uposterior x / nrm

    let mposterior = do
            rprms' <- mrprms
            let mpstr = fromOneHarmonium . rectifiedBayesRule rprms' lkl z $ toOneHarmonium nx0
                muposterior = unnormalizedDensity mpstr
                mnrm = integrate 1e-2000 muposterior 0 (2*pi)
            return (\x -> muposterior x / mnrm)

    let mxlft = 10
        mxrght = 2.5

    layoutlr_x_axis . laxis_title .= "Stimulus/Preferred Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Posterior Density"

    layoutlr_left_axis . laxis_title .= "Spike Count"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    layoutlr_plots
        .= [ Left $ vlinePlot "stimulus" (solidLine 2 $ opaque black) x0 ]

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 2 (opaque red)
        plot_lines_values .= [toList . B.zip pltsmps $ posterior <$> pltsmps]

    case mposterior of
      (Just posterior') ->
          plotRight . liftEC $ do
              plot_lines_title .= "false posterior"
              plot_lines_style .= solidLine 2 (opaque purple)
              plot_lines_values .= [toList . B.zip pltsmps $ posterior' <$> pltsmps]
      Nothing -> return ()

    plotRight . liftEC $ do
        plot_lines_title .= "posterior"
        plot_lines_style .= solidLine 2 (opaque red)
        plot_lines_values .= [toList . B.zip pltsmps $ posterior <$> pltsmps]

    plotLeft . liftEC $ do
        plot_points_title .= "spikes"
        plot_points_style .= filledCircles 3 (opaque black)
        plot_points_values .= toList (B.zip mus z)

-- Main --

main :: IO ()
main = do

    void . goalRenderableToPDF "einstein-retreat-2018" "adaptive-tuning-curves0"
        400 200 . toRenderable .  (layout_title .~ "Homogeneous Population")
        $ tclyt lkl0 rprms0 Nothing
    void . goalRenderableToPDF "einstein-retreat-2018" "adaptive-tuning-curves1"
        400 200 . toRenderable . (layout_title .~ "Heterogeneous Population")
        $ tclyt lkl1 rprms1 (Just rprms0)

    z0 <- realize . samplePoint $ lkl0 >.>* x0
    z1 <- realize . samplePoint $ lkl1 >.>* x0

    void . goalRenderableToPDF "einstein-retreat-2018" "ppc-inference0"
        400 200 . toRenderable . (layoutlr_title .~ "Homogeneous Inference")
        $ byslyt lkl0 rprms0 Nothing z0
    void . goalRenderableToPDF "einstein-retreat-2018" "ppc-inference1"
        400 200 . toRenderable . (layoutlr_title .~ "Heterogeneous Inference")
        $ byslyt lkl0 rprms1 (Just rprms0) z1
--
--    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference1" 300 200 . toRenderable . byslyt $ take 1 rprmzs
--    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference2" 300 200 . toRenderable . byslyt $ take 2 rprmzs
--    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference3" 300 200 . toRenderable $ byslyt rprmzs

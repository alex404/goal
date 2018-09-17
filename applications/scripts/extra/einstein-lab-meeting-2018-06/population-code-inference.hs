{-# LANGUAGE DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import Goal.Geometry
import Goal.Probability


--- Program ---

prjdr :: String
prjdr = "einstein-lab-meeting-2018-06"

-- Globals --

x0 :: Double
x0 = pi + 1

prr :: Natural # VonMises
prr = zero

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

pltsmps :: B.Vector 200 Double
pltsmps = B.range mnx mxx

type NNeurons = 6

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mnx mxx

xsmps :: B.Vector 50 Double
xsmps = B.init $ B.range mnx mxx

kp :: Double
kp = 2

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map (Point . flip S.doubleton kp) mus

gn0 :: Double
gn0 = 0.5

lkl0 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl0 = vonMisesPopulationEncoder sps gn0

rho0 :: Double
rprms0 :: Natural # VonMises
(rho0,rprms0) = populationCodeRectificationParameters lkl0 xsmps

rprms1,rprms2 :: Natural # VonMises
rprms1 = Point $ S.doubleton 1 0
rprms2 = Point $ S.doubleton 2 0

lkl1,lkl2 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl1 = rectifyPopulationCode rho0 rprms1 xsmps lkl0
lkl2 = rectifyPopulationCode rho0 rprms2 xsmps lkl0

-- Functions --

tclyt
    :: Mean ~> Natural # R NNeurons Poisson <* VonMises
    -> Natural # VonMises
    -> B.Vector NNeurons Int
    -> LayoutLR Double Int Double
tclyt lkl rprms rs = execEC $ do

    goalLayoutLR

    let mxlft = 6
        mxrght = 10

    radiansAbscissaLR
    layoutlr_x_axis . laxis_title .= "(Preferred) Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Rate"
    layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Response"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    plotRight . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves pltsmps lkl)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ zip (B.toList pltsmps) . S.toList $ sumOfTuningCurves lkl pltsmps ]

    plotRight . liftEC $ do
        plot_lines_style .= dashedLine 4 [20,20] (opaque red)
        plot_lines_values .=
            [ toList . B.zip pltsmps $ rectificationCurve rho0 rprms pltsmps ]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 5 (opaque black)
        plot_points_values .= zip (S.toList mus) (toList rs)

blflyt :: B.Vector NNeurons Int -> B.Vector NNeurons Int -> B.Vector NNeurons Int -> Layout Double Double
blflyt z0 z1 z2 = execEC $ do

    goalLayout

    let pst1 = fromOneHarmonium . rectifiedBayesRule rprms0 lkl0 z0 $ toOneHarmonium prr
        pst2 = fromOneHarmonium . rectifiedBayesRule rprms1 lkl1 z1 $ toOneHarmonium pst1
        pst3 = fromOneHarmonium . rectifiedBayesRule rprms2 lkl2 z2 $ toOneHarmonium pst2

    let [clr0,clr1,clr2,clr3] = rgbaGradient (1,0,0,0.4) (1,0,0,1) 4

    radiansAbscissa
    layout_x_axis . laxis_title .= "Stimulus"

    layout_y_axis . laxis_title .= "Beliefs"

    plot . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr0
        plot_lines_values .= [toList (B.zip pltsmps $ density prr <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr1
        plot_lines_values .= [toList (B.zip pltsmps $ density pst1 <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr2
        plot_lines_values .= [toList (B.zip pltsmps $ density pst2 <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr3
        plot_lines_values .= [toList (B.zip pltsmps $ density pst3 <$> pltsmps)]


-- Main --


main :: IO ()
main = do

    z0 <- realize . samplePoint $ lkl0 >.>* x0
    z1 <- realize . samplePoint $ lkl1 >.>* x0
    z2 <- realize . samplePoint $ lkl2 >.>* x0

--    void $ goalRenderableToPDF prjdr "simple-population-response" 500 200
--        . toRenderable $ tclyt lkl0 rprms0 z0

    void $ goalRenderableToPDF prjdr "population-response0" 300 150
        . toRenderable $ tclyt lkl0 rprms0 z0
    void $ goalRenderableToPDF prjdr "population-response1" 300 150
        . toRenderable $ tclyt lkl1 rprms1 z1
    void $ goalRenderableToPDF prjdr "population-response2" 300 150
        . toRenderable $ tclyt lkl2 rprms2 z2

    void $ goalRenderableToPDF prjdr "dynamic-beliefs-small" 300 150
        . toRenderable $ blflyt z0 z1 z2

--    void $ goalRenderableToPDF prjdr "dynamic-beliefs" 500 200
--        . toRenderable $ blflyt z0 z1 z2

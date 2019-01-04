{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators #-}

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

prjdr :: String
prjdr = "tuning-curve-sums/tradeoff"

type SampleSize = 1000

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

pltsmps :: B.Vector 200 Double
pltsmps = B.range mnx mxx

type NNeurons = 8

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mnx mxx

xsmps :: B.Vector 50 Double
xsmps = B.init $ B.range mnx mxx

kp :: Double
kp = 1

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map (Point . flip S.doubleton kp) mus

gn0 :: Double
gn0 = 1

lkl0 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl0 = vonMisesPopulationEncoder sps gn0

rho0 :: Double
rprms0 :: Natural # VonMises
(rho0,rprms0) = populationCodeRectificationParameters lkl0 xsmps

rprms1,rprms2 :: Natural # VonMises
rprms1 = Point $ S.doubleton 0 0
rprms2 = Point $ S.doubleton 0 4

lkl1,lkl2 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl1 = rectifyPopulationCode rho0 rprms1 xsmps lkl0
lkl2 = rectifyPopulationCode rho0 rprms2 xsmps lkl0

prr1, prr2 :: Natural # VonMises
prr1 = Point $ S.doubleton 0 4
prr2 = Point $ S.doubleton 0 0

hrm1, hrm2 :: Natural # Harmonium Tensor (R NNeurons Poisson) VonMises
hrm1 = joinBottomHarmonium lkl1 $ toOneHarmonium prr1
hrm2 = joinBottomHarmonium lkl2 $ toOneHarmonium prr2

-- Functions --


tclyt
    :: Mean ~> Natural # R NNeurons Poisson <* VonMises
    -> Natural # VonMises
    -> Layout Double Double
tclyt lkl rprms = execEC $ do

    goalLayout

    radiansAbscissa
    layout_x_axis . laxis_title .= "(Preferred) Stimulus"
    layout_y_axis . laxis_title .= "Rate"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves pltsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ zip (B.toList pltsmps) . S.toList $ sumOfTuningCurves lkl pltsmps ]

    plot . liftEC $ do
        plot_lines_style .= dashedLine 4 [20,20] (opaque red)
        plot_lines_values .=
            [ toList . B.zip pltsmps $ rectificationCurve rho0 rprms pltsmps ]

prrlyt :: Layout Double Double
prrlyt = execEC $ do

    goalLayout
    radiansAbscissa

    layout_x_axis . laxis_title .= "Preferred Stimulus"
    layout_y_axis . laxis_title .= "Belief"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_title .= "Non-Flat Prior + Homogeneous"
        plot_lines_values .= [ B.toList . B.zip pltsmps $ density prr1 <$> pltsmps ]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_title .= "Flat Prior + Imhomogeneous"
        plot_lines_values .= [ B.toList . B.zip pltsmps $ density prr2 <$> pltsmps ]

smplyt
    :: Sample SampleSize (R NNeurons Poisson)
    -> Sample SampleSize (R NNeurons Poisson)
    -> Layout Double Int
smplyt smps1 smps2 = execEC $ do

    goalLayout
    radiansAbscissa

    let smp1 = B.foldl1 (+) smps1
        smp2 = B.foldl1 (+) smps2

    layout_x_axis . laxis_title .= "Preferred Stimulus"
    layout_y_axis . laxis_title .= "Total Spikes"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_title .= "Non-Flat Prior + Homogeneous"
        plot_lines_values .= [ loopRadiansPlotData . B.toList $ B.zip (G.convert mus) smp1 ]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_title .= "Flat Prior + Imhomogeneous"
        plot_lines_values .= [ loopRadiansPlotData . B.toList $ B.zip (G.convert mus) smp2 ]


-- Main --


main :: IO ()
main = do

    (smps1 :: Sample SampleSize (R NNeurons Poisson)) <-
        fmap (B.map hHead) . realize $ sampleRectified (toSingletonSum rprms1) hrm1
    (smps2 :: Sample SampleSize (R NNeurons Poisson)) <-
        fmap (B.map hHead) . realize $ sampleRectified (toSingletonSum rprms2) hrm2

    void . goalRenderableToSVG prjdr "ppc1" 600 300 . toRenderable $ tclyt lkl1 rprms1
    void . goalRenderableToSVG prjdr "ppc2" 600 300 . toRenderable $ tclyt lkl2 rprms2

    void . goalRenderableToSVG prjdr "priors" 600 300 $ toRenderable prrlyt

    void $ goalRenderableToSVG prjdr "sample-comparison" 600 300 . toRenderable $ smplyt smps1 smps2

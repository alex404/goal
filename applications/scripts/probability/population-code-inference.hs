{-# LANGUAGE DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import Goal.Geometry
import Goal.Probability


--- Program ---


-- Globals --

x0 :: Double
x0 = pi + 1

mcts :: Source # Categorical Int 3
mcts = Point $ S.doubleton 0.5 0.2

sp0 :: Source # VonMises
sp0 = Point $ S.doubleton pi 10

sp1 :: Source # VonMises
sp1 = Point $ S.doubleton 1 10

sp2 :: Source # VonMises
sp2 = Point $ S.doubleton 5 10

prr :: Natural # Harmonium Tensor VonMises (Categorical Int 3)
prr = buildMixtureModel (S.map toNatural $ S.fromTuple (sp0,sp1,sp2)) (toNatural mcts)

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

pltsmps :: [Double]
pltsmps = range mnx mxx 200

type NNeurons = 6

type Neurons = R NNeurons Poisson

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mnx mxx

xsmps :: [Double]
xsmps = init $ range mnx mxx 50

kp :: Double
kp = 2

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map (Point . flip S.doubleton kp) mus

gn0 :: Double
gn0 = 0.5

lkl0 :: Mean #> Natural # R NNeurons Poisson <* VonMises
lkl0 = vonMisesPopulationEncoder False (Left gn0) sps

rho0 :: Double
rprms0 :: Natural # VonMises
(rho0,rprms0) = populationCodeRectificationParameters lkl0 xsmps

rprms1,rprms2 :: Natural # VonMises
rprms1 = Point $ S.doubleton 1 0
rprms2 = Point $ S.doubleton 2 0

lkl1,lkl2 :: Mean #> Natural # R NNeurons Poisson <* VonMises
lkl1 = rectifyPopulationCode rho0 rprms1 xsmps lkl0
lkl2 = rectifyPopulationCode rho0 rprms2 xsmps lkl0

numericalPosterior :: [B.Vector NNeurons Int] -> Double -> Double
numericalPosterior zs y =
    let lkls x = [ density (lkl >.>* x) z  |  (lkl,z) <- zip [lkl0,lkl1,lkl2] zs]
        uposterior x = product $ mixtureDensity prr x : lkls x
        nrm = integrate 1e-500 uposterior (-pi) pi
     in uposterior y / nrm

-- Functions --

tclyt
    :: Mean #> Natural # R NNeurons Poisson <* VonMises
    -> Natural # VonMises
    -> B.Vector NNeurons Int
    -> LayoutLR Double Int Double
tclyt lkl rprms rs = execEC $ do

    goalLayoutLR

    let mxlft = 10
        mxrght = 10

    radiansAbscissaLR
    layoutlr_x_axis . laxis_title .= "(Preferred) Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Rate"
    --layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Response"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    plotRight . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tuningCurves pltsmps lkl)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ zip pltsmps $ sumOfTuningCurves lkl pltsmps ]

    plotRight . liftEC $ do
        plot_lines_style .= dashedLine 4 [20,20] (opaque red)
        plot_lines_values .=
            [ zip pltsmps $ rectificationCurve rho0 rprms pltsmps ]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 5 (opaque black)
        plot_points_values .= zip (S.toList mus) (toList rs)

blflyt :: SamplePoint Neurons -> SamplePoint Neurons -> SamplePoint Neurons -> Layout Double Double
blflyt z0 z1 z2 = execEC $ do

    goalLayout

    let pst1 = rectifiedBayesRule rprms0 lkl0 z0 prr
        pst2 = rectifiedBayesRule rprms1 lkl1 z1 pst1
        pst3 = rectifiedBayesRule rprms2 lkl2 z2 pst2
        pst0' = numericalPosterior []
        pst1' = numericalPosterior [z0]
        pst2' = numericalPosterior [z0,z1]
        pst3' = numericalPosterior [z0,z1,z2]

    let [clr0,clr1,clr2,clr3] = rgbaGradient (0,0,1,1) (1,0,0,1) 4

    radiansAbscissa
    layout_x_axis . laxis_title .= "Stimulus"

    layout_y_axis . laxis_title .= "Beliefs"

    plot . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plot . liftEC $ do
        plot_lines_style .= solidLine 6 (black `withOpacity` 0.5)
        plot_lines_values .= [toList (zip pltsmps $ pst0' <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 6 (black `withOpacity` 0.5)
        plot_lines_values .= [toList (zip pltsmps $ pst1' <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 6 (black `withOpacity` 0.5)
        plot_lines_values .= [toList (zip pltsmps $ pst2' <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 6 (black `withOpacity` 0.5)
        plot_lines_values .= [toList (zip pltsmps $ pst3' <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr0
        plot_lines_values .= [toList (zip pltsmps $ mixtureDensity prr <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr1
        plot_lines_values .= [toList (zip pltsmps $ mixtureDensity pst1 <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr2
        plot_lines_values .= [toList (zip pltsmps $ mixtureDensity pst2 <$> pltsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 clr3
        plot_lines_values .= [toList (zip pltsmps $ mixtureDensity pst3 <$> pltsmps)]


-- Main --


main :: IO ()
main = do

    z0 <- realize . samplePoint $ lkl0 >.>* x0
    z1 <- realize . samplePoint $ lkl1 >.>* x0
    z2 <- realize . samplePoint $ lkl2 >.>* x0

    void $ goalRenderableToSVG "probability/population-code-inference" "population-response0" 800 600
        . toRenderable $ tclyt lkl0 rprms0 z0
    void $ goalRenderableToSVG "probability/population-code-inference" "population-response1" 800 600
        . toRenderable $ tclyt lkl1 rprms1 z1
    void $ goalRenderableToSVG "probability/population-code-inference" "population-response2" 800 600
        . toRenderable $ tclyt lkl2 rprms2 z2

    void $ goalRenderableToSVG "probability/population-code-inference" "dynamic-beliefs" 800 600
        . toRenderable $ blflyt z0 z1 z2

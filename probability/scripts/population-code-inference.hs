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
x0 = pi

sp0 :: Source # VonMises
sp0 = Point $ S.doubleton pi 4

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

xsmp :: Int
xsmp = 200

xsmps :: B.Vector 200 Double
xsmps = B.range mnx mxx

mus :: S.Vector 8 Double
mus = S.range mnx mxx

kp :: Double
kp = 2

sps :: S.Vector 7 (Source # VonMises)
sps = S.tail $ S.map (Point . flip S.doubleton kp) mus

gn1,gn2 :: Double
gn1 = 1
gn2 = 0.5

lkl1,lkl2 :: Mean ~> Natural # R 7 Poisson <* VonMises
lkl1 = vonMisesPopulationEncoder sps gn1
lkl2 = vonMisesPopulationEncoder sps gn2

decoder :: Mean ~> Natural # Tensor VonMises (R 7 Poisson)
decoder =
    transpose . snd $ splitAffine lkl1

-- Functions --

normalizedVonMises :: Natural # VonMises -> Double -> Double
normalizedVonMises nx x =
    let uposterior = unnormalizedDensity nx
        nrm = integrate 1e-5000 uposterior 0 (2*pi)
     in uposterior x / nrm

rtrns :: [a] -> [a]
rtrns rs = last rs : rs

lyt1 :: Mean ~> Natural # R 7 Poisson <* VonMises -> B.Vector 7 Int -> LayoutLR Double Int Double
lyt1 lkl rs = execEC $ do

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
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 5 (opaque black)
        plot_points_values .= zip (S.toList mus) (rtrns $ toList rs)

lyt2 :: B.Vector 7 Int -> B.Vector 7 Int -> Layout Double Double
lyt2 rs1 rs2 = execEC $ do

    goalLayout

    let mx = 4

        np0 = toNatural sp0
        np1 = (decoder >.>* rs1) <+> np0
        np2 = (decoder >.>* rs2) <+> np1

    radiansAbscissa
    layout_x_axis . laxis_title .= "Stimulus"

    layout_y_axis . laxis_generate .= scaledAxis def (0,mx)
    layout_y_axis . laxis_title .= "Beliefs"

    plot . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (blue `withOpacity` 0.5)
        plot_lines_values .= [toList (B.zip xsmps $ normalizedVonMises np0 <$> xsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (blue `withOpacity` 0.75)
        plot_lines_values .= [toList (B.zip xsmps $ normalizedVonMises np1 <$> xsmps)]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (blue `withOpacity` 1)
        plot_lines_values .= [toList (B.zip xsmps $ normalizedVonMises np2 <$> xsmps)]


-- Main --


main :: IO ()
main = do

    rs1 <- realize . samplePoint $ lkl1 >.>* x0
    rs2 <- realize . samplePoint $ lkl2 >.>* x0

    void $ goalRenderableToSVG "probability/population-code-inference" "population-response1" 800 600 . toRenderable $ lyt1 lkl1 rs1
    void $ goalRenderableToSVG "probability/population-code-inference" "population-response2" 800 600 . toRenderable $ lyt1 lkl2 rs2
    void $ goalRenderableToSVG "probability/population-code-inference" "dynamic-beliefs" 800 600 . toRenderable $ lyt2 rs1 rs2

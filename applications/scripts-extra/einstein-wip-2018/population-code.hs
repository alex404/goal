{-# LANGUAGE DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import Goal.Geometry
import Goal.Probability

import Goal.NeuralCircuits

--- Program ---


-- Globals --

x0 :: Double
x0 = pi

sp0 :: Source # VonMises
sp0 = Point $ S.doubleton 0 0

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

gn :: Double
gn = 1

lkl :: Mean ~> Natural # R 7 Poisson <* VonMises
lkl = vonMisesPopulationEncoder sps gn

decoding :: B.Vector 7 Int -> Double -> Double
decoding r x = density (lkl >.>* x) r



-- Functions --

rtrns :: [a] -> [a]
rtrns rs = last rs : rs

lyt1 :: B.Vector 7 Int -> LayoutLR Double Int Double
lyt1 rs = execEC $ do

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

lyt2 :: B.Vector 7 Int -> LayoutLR Double Int Double
lyt2 rs = execEC $ do

    goalLayoutLR

    let mxlft = 10
        mxrght = 0.001


    radiansAbscissaLR
    layoutlr_x_axis . laxis_title .= "Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "(Unnormalized) Beliefs"
    layoutlr_right_axis_visibility .= AxisVisibility True False False
    --layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Response"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    plotRight . return $ vlinePlot "" (solidLine 3 $ opaque black) x0

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [toList (B.zip xsmps $ decoding rs <$> xsmps)]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 5 (opaque black)
        plot_points_values .= zip (S.toList mus) (rtrns $ toList rs)



-- Main --

main :: IO ()
main = do

    rs <- realize . sample $ lkl >.>* x0

    void $ goalRenderableToPDF "neural-circuits" "population-response" 500 200 . toRenderable $ lyt1 rs
    void $ goalRenderableToPDF "neural-circuits" "population-decoding" 500 200 . toRenderable $ lyt2 rs

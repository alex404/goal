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


--- Globals ---


-- Stimulus --

x0 :: Double
x0 = pi

-- PPC --

mu :: Double
mu = 2.5

sp :: S.Vector 1 (Source # VonMises)
sp = S.singleton . Point $ S.doubleton mu 2

lkl :: Mean ~> Natural # R 1 Poisson <* VonMises
lkl = vonMisesPopulationEncoder sp 1

-- Plotting --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

xsmp :: Int
xsmp = 200

xsmps :: B.Vector 200 Double
xsmps = B.range mnx mxx

-- Functions --

rtrns :: [a] -> [a]
rtrns rs = last rs : rs

lyt :: Int -> LayoutLR Double Int Double
lyt r = execEC $ do

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
        plot_points_values .= [(mu,r)]


-- Main --

main :: IO ()
main = do

    rs <- realize . sample $ lkl >.>* x0

    void $ goalRenderableToPDF "neural-circuits" "tuning-curve" 500 200 . toRenderable . lyt $ B.head rs

{-# LANGUAGE DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability

--- Program ---


-- Globals --

x0 :: Double
x0 = 0

sp0 :: Source # VonMises
sp0 = Point $ doubleton 0 0

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = -pi
mxx = pi

xsmp :: Int
xsmp = 200

xsmps :: Vector 200 Double
xsmps = rangeV mnx mxx

mus :: Vector 8 Double
mus = rangeV mnx mxx

kp :: Double
kp = 2

sps :: Vector 7 (Source # VonMises)
sps = tailV $ Point . flip doubleton kp <$> mus

gn :: Double
gn = 1

lkl :: Mean ~> Natural # R 7 Poisson <+ VonMises
lkl = vonMisesPopulationEncoder sps gn

-- Functions --

rtrns :: [a] -> [a]
rtrns rs = last rs : rs

lyt :: Vector 7 Int -> LayoutLR Double Int Double
lyt rs = execEC $ do

    goalLayoutLR

    let axs = [-pi,-pi/2,0,pi/2,pi]
        mxlft = 10
        mxrght = 10

    layoutlr_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(-pi,"-π"),(-pi/2,"-π/2"),(0,"0"),(pi/2,"π/2"),(pi,"π")]
    layoutlr_x_axis . laxis_title .= "(Preferred) Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Rate"
    --layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Response"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (blue `withOpacity` 0.7)
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    layoutlr_plots
        %= ( ++ [Left $ vlinePlot "" (solidLine 3 $ opaque red) x0])

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 5 (opaque black)
        plot_points_values .= zip (toList mus) (rtrns $ toList rs)



-- Main --

main :: IO ()
main = do

    rs <- realize . generate $ lkl >.>* x0

    void $ goalRenderableToPDF "neural-circuits" "population-response" 500 200 . toRenderable $ lyt rs

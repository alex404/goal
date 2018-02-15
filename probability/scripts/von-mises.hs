{-# LANGUAGE TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Program ---


-- Globals --

nsmps :: Int
nsmps = 1000

mu,kap :: Double
mu = -2
kap = 2

tru :: Source # VonMises
tru = Point $ doubleton mu kap

-- Plot

mn,mx :: Double
(mn,mx) = (-pi,pi)

xs :: [Double]
xs = range mn mx 200

nb :: Int
nb = 25


-- Main --

main :: IO ()
main = do

    smps <- realize . replicateM nsmps $ generate tru

    let lyt = execEC $ do

            histogramLayoutLR nb mn mx

            layoutlr_title .= "Sampling Algorithm vs . Unnormalized von Mises Density"
            layoutlr_left_axis . laxis_title .= "Sample Count"
            layoutlr_left_axis . laxis_override .= axisGridHide

            layoutlr_right_axis . laxis_title .= "Unnormalized Probability Mass"
            layoutlr_right_axis . laxis_override .= axisGridHide

            layoutlr_x_axis . laxis_title .= "Value"
            layoutlr_x_axis . laxis_override .= axisGridHide

            plotLeft . fmap plotBars . liftEC $ do
                histogramPlot nb mn mx [smps]
                plot_bars_titles .= ["Samples"]
                plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

            plotRight . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_title .= "Density"
                plot_lines_values .= [ [(x,unnormalizedDensity (transition tru) x)  | x <- xs] ]

    goalRenderableToSVG "probability" "von-mises" 800 400 $ toRenderable lyt


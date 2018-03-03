{-# LANGUAGE TypeFamilies #-}

--- Imports ---

import Attractor

-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability


--- Globals ---


x0 = 0
stps = 200
pltrng = range mnx mxx stps


--- Main ---


main = do

    rs <- runWithSystemRandom . standardGenerate $ emsn >.>* x0

    let tclyt = toRenderable . execEC $ do

            layout_y_axis . laxis_title .= "Rate"
            layout_x_axis . laxis_title .= "Stimulus"

            plot . liftEC $ do
                plot_lines_style .= dashedLine 2 [6,3] (opaque blue)
                plot_lines_values .= ( zip pltrng <$> transpose
                    (listCoordinates . (gn />) . dualTransition <$> (emsn >$>* pltrng)) )

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip pltrng $ (/fromIntegral nn) . sum <$>
                    (listCoordinates . (gn />) . dualTransition <$> (emsn >$>* pltrng))]

    let rsplyt = toRenderable . execEC $ do

            let posterior = dcdn >.>* rs

            layoutlr_left_axis . laxis_title .= "Probability Density"
            layoutlr_right_axis . laxis_title .= "Response Count"
            layoutlr_x_axis . laxis_title .= "Stimulus"

            layoutlr_plots
                .= [ Left $ vlinePlot "" (solidLine 2 $ opaque black) x0 ]

            plotRight . liftEC $ do
                plot_points_style .= filledCircles 4 (opaque black)
                plot_points_values .= zip rngn rs

            plotLeft . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque red)
                plot_lines_values .= [zip pltrng $ density posterior <$> pltrng]

    let rnbl = gridToRenderable . weights (1,1) $ (tval tclyt ./. tval rsplyt)

    goalRenderableToSVG sbdr "population-code" 800 800 rnbl



{-# LANGUAGE TypeFamilies #-}

--- Imports ---

import Pendulum

-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Globals ---


x0 = (0,0)
stps = 200
chnstps = 10^5
qrng = range mnq mxq stps
dqrng = range mndq mxdq stps


--- Main ---


main = do

    rs <- runWithSystemRandom . standardGenerate $ emsn >.>* x0
    chn <- runWithSystemRandom $ pendulumChain x0

    let qtclyt = toRenderable . execEC $ do

            layout_y_axis . laxis_title .= "Rate"
            layout_x_axis . laxis_title .= "Stimulus"
            layout_title .= "Angle"

            let qtcs' = take nnq . listCoordinates . (gnq />) . dualTransition <$> (emsn >$>* zip qrng (repeat 0))

            plot . liftEC $ do
                plot_lines_style .= dashedLine 2 [6,3] (opaque blue)
                plot_lines_values .= ( zip qrng <$> transpose qtcs' )

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip qrng $ (/fromIntegral nnq) . sum <$> qtcs']

    let dqtclyt = toRenderable . execEC $ do

            layout_y_axis . laxis_title .= "Rate"
            layout_x_axis . laxis_title .= "Stimulus"
            layout_title .= "Angular Velocity"

            let dqtcs' = drop nnq . listCoordinates . (gndq />) . dualTransition <$> (emsn >$>* zip (repeat 0) dqrng)

            plot . liftEC $ do
                plot_lines_style .= dashedLine 2 [6,3] (opaque blue)
                plot_lines_values .= ( zip dqrng <$> transpose dqtcs' )

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip dqrng $ (/fromIntegral nndq) . sum <$> dqtcs']

{-
    let rsplyt = toRenderable . execEC $ do

            let posterior = dcdn >.>* rs

            layoutlr_left_axis . laxis_title .= "Unnormalized Probability Density"
            layoutlr_right_axis . laxis_title .= "Response Count"
            layoutlr_x_axis . laxis_title .= "Stimulus"

            layoutlr_plots
                .= [ Left $ vlinePlot "" (solidLine 2 $ opaque black) x0 ]

            plotRight . liftEC $ do
                plot_points_style .= filledCircles 4 (opaque black)
                plot_points_values .= zip rngq rs

            plotLeft . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque red)
                plot_lines_values .= [zip pltrng [ exp $ sufficientStatistic VonMises q <.> posterior | q <- pltrng] ]
                -}

    let rnbl = gridToRenderable . weights (1,1) $ (tval qtclyt ./. tval dqtclyt)

    goalRenderableToSVG sbdr "population-code" 800 800 rnbl

    let (qs,dqs) = unzip . take chnstps $ streamChain chn

    putStr "Max Angle:"
    print $ maximum qs
    putStr "Min Angle:"
    print $ minimum qs
    putStr "Max Angular Velocity:"
    print $ maximum dqs
    putStr "Min Angular Velocity:"
    print $ minimum dqs



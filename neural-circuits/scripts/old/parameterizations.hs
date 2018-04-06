{-# LANGUAGE ScopedTypeVariables,TypeOperators,DataKinds #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Program ---


-- Globals --

-- Population Code

nrm1,nrm2,nrm3 :: Source # Normal
nrm1 = Point $ doubleton (-3) 1
nrm2 = Point $ doubleton 0 1.5
nrm3 = Point $ doubleton 3 1

nbns :: Int
nbns = 20

mn,mx :: Double
mn = -10
mx = 10

-- Main --

main :: IO ()
main = do

    smps1 :: Vector 20 Double <- realize . replicateMV $ generate nrm1
    smps2 :: Vector 20 Double <- realize . replicateMV $ generate nrm2
    smps3 :: Vector 20 Double <- realize . replicateMV $ generate nrm3

    let smps = joinV smps1 (joinV smps2 smps3)

    let hstlyt = execEC $ do

            goalLayout
            histogramLayout 10 mn mx

            layout_y_axis . laxis_title .= "Sample Count"
            layout_x_axis . laxis_title .= "x"

            plot . fmap plotBars . liftEC $ do
                histogramPlot nbns mn mx [concat $ toList <$> [smps1,smps2,smps3]]
                plot_bars_item_styles .= [(solidFillStyle $ opaque red, Nothing)]
                plot_bars_spacing .= BarsFixGap 2 0

    let nrm :: Source # Normal
        nrm = transition $ sufficientStatisticT smps

    let nrmlyt = execEC $ do

            goalLayout

            let pltxs = range (mn+0.5) (mx-0.5) 500

            layout_y_axis . laxis_title .= "Density"
            layout_x_axis . laxis_title .= "x"

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [zip pltxs $ density nrm <$> pltxs]


    goalRenderableToPDF "neural-circuits/parameterizations" "histogram" 250 200 $ toRenderable hstlyt
    goalRenderableToPDF "neural-circuits/parameterizations" "normal" 250 200 $ toRenderable nrmlyt

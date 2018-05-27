{-# LANGUAGE TypeOperators,DataKinds #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Program ---


-- Globals --

-- Population Code

vr0 = 0.75
mu0 = 2
sp0 :: Source # Normal
sp0 = Point $ doubleton 0.75 2

vr,mn,mx :: Double
vr = 2
mn = -6
mx = 6

mus :: Vector 20 Double
mus = rangeV mn mx

sps :: Vector 20 (Source # Normal)
sps = Point . flip doubleton vr <$> mus

gn :: Double
gn = 5

ppc = normalPopulationEncoder sps gn
hrm = affineToHarmonium (transition sp0) ppc

-- Samples
nstps = 500
x0s :: Vector 500 Double
x0s = rangeV (-3) 3

-- Plot
(pltmn,pltmx) = (mn, mx)
pltxs = range mn mx 200
nb = 20

-- Main --

main = do

    chn <- realize $ bulkGibbsChain (harmoniumTranspose hrm) x0s
    let chnstps = take nstps $ streamChain chn
    let xsmps = snd <$> last chnstps


    let lyt = execEC $ do

            histogramLayoutLR nb mn mx

            layoutlr_title .= "Gibbs Sampling on Population Codes"
            layoutlr_left_axis . laxis_title .= "Sample Count"
            layoutlr_left_axis . laxis_override .= axisGridHide

            layoutlr_right_axis . laxis_title .= "Density"
            layoutlr_right_axis . laxis_override .= axisGridHide

            layoutlr_x_axis . laxis_title .= "Value"
            layoutlr_x_axis . laxis_override .= axisGridHide

            plotLeft . fmap plotBars . liftEC $ do
                histogramPlot nb pltmn pltmx [toList xsmps]
                plot_bars_titles .= ["Samples"]
                plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

            plotRight . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_title .= "Approximate Prior"
                plot_lines_values .= [ zip pltxs $ density sp0 <$> pltxs ]

    goalRenderableToSVG "simulation" "gibbs" 800 400 $ toRenderable (lyt :: LayoutLR Double Int Double)

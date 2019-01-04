{-# LANGUAGE FlexibleContexts,TypeOperators,DataKinds #-}
-- A script which demonstrates how the binomial and poisson distributions
-- approximate each other.

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

--- Globals ---

lmda1 :: Double
lmda1 = 5

rng1,rng2 :: [Int]
rng1 = [0..20]
rng2 = [0..10]


--- Function ---


stringedPlot :: [(Int,Double)] -> Colour Double -> Double -> Double -> EC (Layout Int Double) ()
stringedPlot pnts clr sz1 sz2 = do

    plot . liftEC $ do

        plot_points_style .= filledCircles sz1 (opaque clr)
        plot_points_values .= pnts

    plot . liftEC $ do

        plot_lines_style .= solidLine sz2 (opaque clr)
        plot_lines_values .= [pnts]

labels :: Show x => [x] -> [(x,String)]
labels xs = zip xs $ show <$> xs

--- Script ---

main :: IO ()
main = do

    let lyt1 p = execEC $ do

            goalLayout

            layout_y_axis . laxis_title .= "Probability Mass"
            layout_x_axis . laxis_title .= "Count"
            layout_x_axis . laxis_generate .= const (makeAxis (map show) ([-1,5],[-1,5],[-1,5]))
            layout_y_axis . laxis_generate .= scaledAxis def (0,1)
            layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride (labels [0,2,4])

            let pd :: Source # Poisson
                pd = Point $ S.singleton p
                ppnts = zip rng2 $ density pd <$> rng2

            stringedPlot ppnts red 5 3

            let bd :: Source # Bernoulli
                bd = Point $ S.singleton p
                bpnts = zip rng2 $ density bd <$> [False,True]

            stringedPlot bpnts black 3 2

    let lyt2 = execEC $ do

            goalLayout

            layout_y_axis . laxis_title .= "Probability Mass"
            layout_x_axis . laxis_title .= "Count"
            layout_x_axis . laxis_generate .= const (makeAxis (map show) ([-1,16],[-1,16],[-1,16]))
            layout_y_axis . laxis_generate .= scaledAxis def (0,0.3)
            layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride (labels [0,5,10,15])

            let pd :: Source # Poisson
                pd = Point $ S.singleton lmda1
                ppnts = zip rng1 $ density pd <$> rng1

            let bd1 :: Source # Binomial 10
                bd1 = Point $ S.singleton (lmda1 / 10)

                bd2 :: Source # Binomial 100
                bd2 = Point $ S.singleton (lmda1 / 100)

            let bplt bd = do

                    let bpnts = zip rng1 $ density bd <$> take (binomialTrials bd+1) rng1
                    stringedPlot bpnts black 3 2

            stringedPlot ppnts red 5 3
            bplt bd1
            --bplt 20
            bplt bd2


    goalRenderableToSVG "probability/poisson-binomial" "poisson-bernoulli1" 150 150 . toRenderable $ lyt1 0.1
    goalRenderableToSVG "probability/poisson-binomial" "poisson-bernoulli2" 150 150 . toRenderable $ lyt1 0.2
    goalRenderableToSVG "probability/poisson-binomial" "poisson-bernoulli3" 150 150 . toRenderable $ lyt1 0.4
    goalRenderableToSVG "probability/poisson-binomial" "poisson-binomial" 500 150 $ toRenderable lyt2

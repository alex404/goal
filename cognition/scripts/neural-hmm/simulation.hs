{-# LANGUAGE TypeOperators #-}

--- Imports ---


-- Goal --

import NeuralHMM

import Goal.Core hiding (Down)
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition


--- Globals ---


nstps = 30

main = do

    gef1 <- fromList mg1 . read <$> goalReadFile sbdr "multilayer-perceptronef1"
    alss0 <- read <$> goalReadFile sbdr "descent"
    let [lns,lopts,lef1,lcd1,_,_,lef3,lcd3] = alss0 :: [[Double]]
        alss = [lns,lopts,lef1,lcd1,lef3,lcd3]

    let x0 = Green
        p0 = transitionTo Natural $ Standard # fromList mx [0.33,0.33]
    xchn <- runWithSystemRandom $ markovChain x0
    nmly <- runWithSystemRandom responseMealy
    let gflt = harmoniumEncodedFilter y01 amtx1 bmtx1 gef1
        optflt = parametricFilter (transition . discretePrediction trns . transition) (rectifiedBayesRule emsn (zero mx)) p0

    let optchn = filterChain xchn nmly optflt
        (xs,ns,optblfs) = unzip3 . take nstps $ streamChain optchn
        zs = take nstps . streamChain $ xchn >>> nmly >>> gflt

    let lytrsp = execEC $ do

            goalLayout
            layout_y_axis . laxis_override .= (axisGridHide . axisLabelsOverride [(-3,"Red"),(3,"Blue")])
            layout_y_axis . laxis_generate .= scaledIntAxis defaultIntAxis (-5,5)
            layout_x_axis . laxis_title .= "Step"

            let ks = [1..] :: [Int]

            sequence_ $ do

                (k,n) <- zip ks ns
                (i,spk) <- zip [-4,-3..] n

                return . plot . liftEC $ do
                    plot_points_style .= filledPolygon 3 4 True (opaque black)
                    when (spk /= 0) $ plot_points_values .= [(k,i :: Int)]

    let lytsm = execEC $ do

            goalLayout
            layout_y_axis . laxis_generate .= scaledIntAxis defaultIntAxis (-5,5)
            layout_x_axis . laxis_title .= "Step"
            layout_y_axis . laxis_override .= (axisGridHide . axisLabelsOverride [(-3,"Learned"),(0,"Stimulus"),(3,"Optimal")])

            let ks = [1..] :: [Int]

            sequence_ $ do

                x <- xcat
                let ks' = fst <$> filter ((== x) . snd) (zip ks xs)

                return . plot . liftEC $ do
                    plot_points_style .= filledCircles 2 (opaque $ clrs x)
                    plot_points_values .= zip ks' (repeat 0)

            sequence_ $ do

                (k,optblf) <- zip ks optblfs
                x <- xcat

                return . plot . liftEC $ do
                    plot_points_style .= filledPolygon 3 4 False (clrs x `withOpacity` density optblf x)
                    plot_points_values .= [(k,optpos x)]

            sequence_ $ do

                (k,z) <- zip ks zs
                x <- xcat

                return . plot . liftEC $ do
                    plot_points_style .= filledPolygon 3 4 False (clrs x `withOpacity` density (dcdz1 >.> z) x)
                    plot_points_values .= [(k,estpos x)]

    let lnlyt = execEC $ do

            goalLayout
            layout_y_axis . laxis_generate .= scaledAxis def (0.75,1)
            layout_x_axis . laxis_title .= "Epoch"
            layout_y_axis . laxis_title .= "Average -Log-Likelihood"
            layout_legend .= Nothing --Just (legend_orientation .~ LOCols 3 $ def)

            sequence_ $ do

                (ttl,ln,ls) <- zip3 ttls pltlns alss

                return . plot . liftEC $ do

                    plot_lines_title .= ttl
                    plot_lines_style .= ln
                    plot_lines_values .= [zip [(1 :: Int)..] ls]

    putStrLn "Percent of Optimum Achieved:"
    print $ (last lef1 - last lns) / (last lopts - last lns)

    goalRenderableToPDF sbdr "neural-hmm-descent" 250 250 $ toRenderable lnlyt
    goalRenderableToPDF sbdr "neural-hmm-response-simulation" 230 113 $ toRenderable lytrsp
    goalRenderableToPDF sbdr "neural-hmm-simulation" 250 125 $ toRenderable lytsm

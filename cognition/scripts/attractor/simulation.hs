{-# LANGUAGE Arrows #-}

-- Goal --

import Attractor

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition


--- Program ---


-- Globals --

x0 = 0
nstps = 99
tcstps = 99
ts' = take nstps ts
dtnm = "multilayer-perceptronef1"
sbdr' = sbdr ++ "/tuning-curves"

-- Main --

main :: IO ()
main = do

    bl <- goalDoesFileExist sbdr dtnm
    g <- if bl
              then fromList mg1 . read <$> goalReadFile sbdr dtnm
              else error $ "Script requires a " ++ dtnm ++ " file"
    alss0 <- read <$> goalReadFile sbdr "descent"

    let [lns,lopts,lef1,lcd1,_,_,lef3,lcd3] = alss0 :: [[Double]]
        alss = [lns,lopts,lef1,lcd1,lef3,lcd3]

    xchn <- runWithSystemRandom (attractorChain x0)
    nmly <- runWithSystemRandom attractorResponseMealy

    let gflt = harmoniumEncodedFilter y01 amtx1 bmtx1 g
        schn = proc () -> do
            x <- xchn -< ()
            n <- nmly -< x
            z <- gflt -< n
            kblf <- kflt -< n
            returnA -< (x,n,z,kblf)
        (xs,ns,zs,kblfs) = unzip4 . take tcstps $ streamChain schn
        (mn',mx') = (-4,4)
        smps = range mn' mx' 100
    --tcs <- runWithSystemRandom $ initialPredictedTuningCurves emsn amtx1 smps 100 g
    let bnsz = 0.25
        mnbn = -3
        mxbn = 3
        bns = [mnbn,mnbn+bnsz..mxbn]
        bls = [\x -> x > bn-bnsz/2 && x < bn+bnsz/2 | bn <- bns]
        tcs' = predictionTuningCurves (zip xs zs) bls g

    let lnlyt = execEC $ do

            goalLayout
            layout_y_axis . laxis_generate .= scaledAxis def (0,2)
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Avg. -Log-Likelihood"
            layout_legend .= Nothing

            sequence_ $ do

                (ttl,pltln,ls) <- zip3 ttls pltlns (alss :: [[Double]])

                return . plot . liftEC $ do

                    plot_lines_title .= ttl
                    plot_lines_style .= pltln
                    plot_lines_values .= [zip [(1 :: Int)..] ls]

    let tclyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Position"
            layout_y_axis . laxis_title .= "Activation"
            layout_x_axis . laxis_generate .= scaledAxis def (mn',mx')
            layout_x_axis . laxis_generate .= scaledAxis def (mnbn,mxbn)

            sequence_ $ do

                (k,clr) <- zip [1,101,114,6,8,10] defaultColorSeq

                return . plot . liftEC $ do

                    plot_lines_values .= [zip bns (tcs' !! k)]
                    plot_lines_style .= solidLine 2 clr

{-
    let tclytk k = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Position"
            layout_y_axis . laxis_title .= "Activation"
            layout_x_axis . laxis_generate .= scaledAxis def (mn',mx')

            plot . liftEC $ do

                plot_lines_values .= [zip smps (tcs !! k)]
                plot_lines_style .= solidLine 2 (opaque black)
                -}

    let tclytk' k = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Position"
            layout_y_axis . laxis_title .= "Activation"
            layout_x_axis . laxis_generate .= scaledAxis def (mn',mx')

            plot . liftEC $ do

                plot_lines_values .= [zip bns (tcs' !! k)]
                plot_lines_style .= solidLine 2 (opaque black)

    let smlyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Time"
            layout_y_axis . laxis_title .= "Position"
            layout_legend .= Nothing

            plot . liftEC $ do

                plot_lines_values .= [zip ts' xs]
                plot_lines_style .= solidLine 2 (opaque black)
                plot_lines_title .= "True"

            plot . liftEC $ do

                plot_points_style .= filledCircles 2 (opaque black)
                plot_points_values .= zip ts' (coordinate 0 . transitionTo Standard <$> (dcdn >$>* ns))
                plot_points_title .= "Observations"

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (red `withOpacity` 0.6)
                plot_lines_values .= [zip ts' (coordinate 0 . transitionTo Standard <$> (dcdz1 >$> zs))]
                plot_lines_title .= "Learned"

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (green `withOpacity` 0.6)
                plot_lines_values .= [zip ts' (coordinate 0 . transitionTo Standard <$> kblfs)]
                plot_lines_title .= "Optimal"

    putStrLn "Percent of Optimum Achieved:"
    print $ 100 * (last (alss !! 2) - last (head alss)) / (last (alss !! 1) - last (head alss))
    goalRenderableToPDF sbdr "attractor-tuning-curves" 300 150 $ toRenderable tclyt
    --sequence_ [ goalRenderableToPDF sbdr' ("attractor-tuning-curves" ++ show k) 300 150 $ toRenderable (tclytk k) | k <- [0..499]]
    sequence_ [ goalRenderableToPDF sbdr' ("attractor-tuning-curves" ++ show k ++ "'") 250 125 $ toRenderable (tclytk' k) | k <- [0..199]]
    goalRenderableToPDF sbdr "attractor-simulation" 250 250 $ toRenderable smlyt
    goalRenderableToPDF sbdr "attractor-descent" 250 125 $ toRenderable lnlyt
    print $ maximum xs
    print $ minimum xs

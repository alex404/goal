{-# LANGUAGE Arrows,FlexibleContexts #-}


-- Goal --

import Pendulum

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition

-- Qualified --


--- Program ---


-- Globals --

nstps = 199
tcstps = 199 --tcstps = 200000
axs = [pi,pi/2,0,-pi/2,-pi]
ts' = take nstps ts
dtnm = "multilayer-perceptroncd1"
sbdr' = sbdr ++ "/tuning-curves"
x0 = (2,4)



-- Functions --


-- Main --

main :: IO ()
main = do

    bl <- goalDoesFileExist sbdr dtnm
    g <- if bl
              then fromList mg1 . read <$> goalReadFile sbdr dtnm
              else error $ "Script requires a " ++ dtnm ++ " file"

    alss0 <- read <$> goalReadFile sbdr "descent"

    let [lns,lekfs,_,lef1,lcd1,lef3,lcd3] = alss0 :: [[Double]]
        alss = [lns,lekfs,lef1,lcd1,lef3,lcd3]

    xchn <- runWithSystemRandom (pendulumChain x0)
    nmly <- runWithSystemRandom pendulumResponseMealy

    let gflt = harmoniumEncodedFilter y01 amtx1 bmtx1 g
        schn = proc () -> do
            x <- xchn -< ()
            n <- nmly -< x
            z <- gflt -< n
            kblf <- kflt -< n
            ekblf <- ekflt -< n
            returnA -< (x,n,z,kblf,ekblf)
        (xs,ns,zs,kblfs,ekblfs) = unzip5 . take tcstps $ streamChain schn

    let qbnsz = pi/10
        qbnsz2 = qbnsz/2
        (mnqbn,mxqbn) = (mnq+qbnsz2,mxq-qbnsz2)
        qbns = [mnqbn,mnqbn+qbnsz..mxqbn]
        dqbnsz = 2
        dqbnsz2 = dqbnsz/2
        (mndqbn,mxdqbn) = (mndq+dqbnsz2+6,mxdq-dqbnsz2-6)
        dqbns = [mndqbn,mndqbn+dqbnsz..mxdqbn]
        bns = [ (qbn,dqbn) | dqbn <- dqbns, qbn <- qbns ]
        bls = [\(q,dq) -> q > qbn-qbnsz2 && q < qbn+qbnsz2 && dq > dqbn-dqbnsz2 && dq < dqbn+dqbnsz2 | (qbn,dqbn) <- bns]
        --tcs' = breakEvery (length qbns) <$> predictionTuningCurves (zip xs zs) bls g

    --tcs0 <- runWithSystemRandom $ initialPredictedTuningCurves emsn amtx2 bns 100 g
    --let tcs = breakEvery (length qbns) <$> tcs0

    let lnlyt = execEC $ do

            goalLayout
            layout_legend .= Nothing
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Mean Squared Error"

            sequence_ $ do

                (ttl,pltln,ls) <- zip3 ttls pltlns (alss :: [[Double]])

                return . plot . liftEC $ do

                    plot_lines_title .= ttl
                    plot_lines_style .= pltln
                    plot_lines_values .= [zip [(1 :: Int)..] ls]

{-
    let tclyt tcs k = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Angle"
            layout_y_axis . laxis_title .= "Activation"
            layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
            layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(pi,"π"),(pi/2,"π/2"),(0,"0"),(-pi/2,"-π/2"),(-pi,"-π")]

            sequence_ $ do

                let tcsk = tcs !! k
                (tc,clr) <- zip tcsk . rgbaGradient (0,0,0,1) (1,0,0,1) $ length tcsk

                return . plot . liftEC $ do

                    plot_lines_values .= [zip qbns tc]
                    plot_lines_style .= solidLine 2 clr
                    -}

    let smlytq = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Time"
            layout_y_axis . laxis_title .= "Angle"
            layout_y_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
            layout_y_axis . laxis_override .= axisGridHide . axisLabelsOverride [(pi,"π"),(pi/2,"π/2"),(0,"0"),(-pi/2,"-π/2"),(-pi,"-π")]

            plot . liftEC $ do

                plot_lines_values .= periodicBreaker (zip ts' $ fst <$> xs)
                plot_lines_style .= solidLine 2 (opaque black)

            plot . liftEC $ do

                plot_points_style .= filledCircles 2 (opaque black)
                plot_points_values .= zip ts' (coordinate 0 . transitionTo Standard <$> (dcdn >$>* ns))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (red `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 0 . transitionTo Standard <$> (dcdz1 >$> zs)))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (green `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 0 . transitionTo Standard <$> ekblfs))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (purple `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 0 . transitionTo Standard <$> kblfs))

    let smlytdq = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Time"
            layout_y_axis . laxis_title .= "Angular Velocity"

            plot . liftEC $ do

                plot_lines_values .= [zip ts' $ snd <$> xs]
                plot_lines_style .= solidLine 2 (opaque black)

            plot . liftEC $ do

                plot_points_style .= filledCircles 2 (opaque black)
                plot_points_values .= zip ts' (coordinate 2 . transitionTo Standard <$> (dcdn >$>* ns))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (red `withOpacity` 0.6)
                plot_lines_values .= [zip ts' (coordinate 2 . transitionTo Standard <$> (dcdz1 >$> zs))]

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (green `withOpacity` 0.6)
                plot_lines_values .= [zip ts' (coordinate 2 . transitionTo Standard <$> ekblfs)]

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (purple `withOpacity` 0.6)
                plot_lines_values .= [zip ts' (coordinate 2 . transitionTo Standard <$> kblfs)]

    let prslyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Time"
            layout_y_axis . laxis_title .= "Precision"

            plot . liftEC $ do

                plot_points_style .= filledCircles 2 (opaque black)
                plot_points_values .= zip ts' (coordinate 1 . transitionTo Standard <$> (dcdn >$>* ns))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (red `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 1 . transitionTo Standard <$> (dcdz1 >$> zs)))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (green `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 1 . transitionTo Standard <$> ekblfs))

            plot . liftEC $ do

                plot_lines_style .= solidLine 2 (blue `withOpacity` 0.6)
                plot_lines_values .= periodicBreaker (zip ts' (coordinate 1 . transitionTo Standard <$> kblfs))


    putStrLn "Percent of Optimum Achieved:"
    print $ 100 * (last (alss !! 2) - last (head alss)) / (last (alss !! 1) - last (head alss))
    --sequence_ [ goalRenderableToPDF sbdr' ("pendulum-tuning-curves" ++ show k) 300 150 $ toRenderable (tclyt tcs k) | k <- [0..499]]
    --sequence_ [ goalRenderableToPDF sbdr' ("pendulum-tuning-curves" ++ show k ++ "'") 300 150 $ toRenderable (tclyt tcs' k) | k <- [0..499]]
    goalRenderableToPDF sbdr "pendulum-angle-simulation" 600 150 $ toRenderable smlytq
    goalRenderableToPDF sbdr "pendulum-velocity-simulation" 600 150 $ toRenderable smlytdq
    goalRenderableToPDF sbdr "pendulum-variance" 300 150 $ toRenderable prslyt
    goalRenderableToPDF sbdr "pendulum-descent" 300 300 $ toRenderable lnlyt
    print . minimum $ fst <$> xs
    print . maximum $ fst <$> xs
    print . minimum $ snd <$> xs
    print . maximum $ snd <$> xs

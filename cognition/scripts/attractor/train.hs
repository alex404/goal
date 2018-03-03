{-# LANGUAGE Arrows,FlexibleContexts,DataKinds,TypeFamilies,TypeOperators #-}

--- Imports ---


-- Goal --

import Attractor

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition

--import Control.Parallel.Strategies

--- Program ---


-- Globals --

epss = iterate (/1.25) 0.00005
--cdns = concat . transpose . take 1 $ repeat [1,2..]
--nstpzs = [1..]
nstpzs = (2^) <$> [0..]
trnchnln = 10000
nepch = 20
-- Adam
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
-- Momentum
--mumx = 0.999

-- Histogram
hstchnln = 200000
nbns = 20
llmnbn = -1
llmxbn = 4

-- Filenames
sbdr' = sbdr ++ "/train"
mlpfl = "multilayer-perceptron"
hstfl = "histogram"
dstfl = "descent"

-- Functions --

filterFun l = not (isInfinite l) && not (isNaN l)

epoch (ncs,alss0) (epch,eps,nstpz) = do

    -- Log-Likelihood Histogram
    --x0 <- runWithSystemRandom randomAttractorState
    let x0 = 0
    xchn <- runWithSystemRandom $ attractorChain x0
    nmly <- runWithSystemRandom attractorResponseMealy

    ncrcmlys <- runWithSystemRandom . sequence $ neuralCircuitTrainer eps bt1 bt2 rg nstpz <$> ncs

    let ncrcschn = streamChain $ xchn >>> nmly >>> mapMealy ncrcmlys
        ncrcs' = ncrcschn !! trnchnln
        hrmflts = neuralCircuitFilter <$> ncrcs'

    -- Log-Likelihood Histogram
    --x0' <- runWithSystemRandom randomAttractorState
    let x0' = 0
    xchn' <- runWithSystemRandom $ attractorChain x0'
    nmly' <- runWithSystemRandom attractorResponseMealy

    let hstchn = proc () -> do
            x <- xchn' -< ()
            n <- nmly' -< x
            blfs <- mapMealy (kflt:hrmflts) -< n
            returnA -< negate . log . flip density x <$> (dcdn >.>* n):blfs

        lss = transpose . take hstchnln $ streamChain hstchn

        hstlyt = execEC $ do

            histogramLayout nbns llmnbn llmxbn
            layout_x_axis . laxis_title .= "-Log-Likelihood"
            layout_y_axis . laxis_title .= "Samples (Thousands)"

            plot . fmap plotBars . liftEC $ do

                void $ histogramPlot nbns llmnbn llmxbn lss
                plot_bars_item_styles .= [ (FillStyleSolid $ opaque clr,Nothing) | clr <- pltclrs ]
                plot_bars_titles .= ttls

    -- Plot
    let flnm = show epch ++ hstfl

    goalRenderableToSVG sbdr' flnm 800 400 $ toRenderable hstlyt
    putStrLn $ "Wrote " ++ flnm ++ ".svg"

    -- Averages
    let als1 = [ average $ filter filterFun ls | ls <- lss ]

    sequence_ $ do

        (al,ls,ttl) <- zip3 als1 lss ttls

        return $ do

            putStr $ ttl ++ " average log-likelihood: "
            print al
            putStr $ ttl ++ " number of singular beliefs: "
            print . length $ filter (not . filterFun) ls
            putStrLn ""


    return (increment <$> ncrcs',als1:alss0)

-- Main --

main :: IO ()
main = do

    g1ef0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg1
    g1cd0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg1
    g2ef0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg2
    g2cd0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg2
    g3ef0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg3
    g3cd0 <- runWithSystemRandom $ initialize (Standard # fromList Normal [0,0.00001]) mg3

    let ncrcef1 = HarmoniumCircuit emsn g1ef0 amtx1 bmtx1 dcdy1 dcdz1 y01 Nothing
        ncrccd1 = HarmoniumCircuit emsn g1cd0 amtx1 bmtx1 dcdy1 dcdz1 y01 (Just 1)
        ncrcef2 = HarmoniumCircuit emsn g2ef0 amtx2 bmtx2 dcdy2 dcdz2 y02 Nothing
        ncrccd2 = HarmoniumCircuit emsn g2cd0 amtx2 bmtx2 dcdy2 dcdz2 y02 (Just 1)
        ncrcef3 = HarmoniumCircuit emsn g3ef0 amtx3 bmtx3 dcdy3 dcdz3 y03 Nothing
        ncrccd3 = HarmoniumCircuit emsn g3cd0 amtx3 bmtx3 dcdy3 dcdz3 y03 (Just 1)
        ncs0 = [NC1 ncrcef1,NC1 ncrccd1,NC1 ncrcef2,NC1 ncrccd2,NC1 ncrcef3,NC1 ncrccd3]

    ([NC1 ncrcef1',NC1 ncrccd1',NC1 ncrcef2',NC1 ncrccd2',NC1 ncrcef3',NC1 ncrccd3'],alss0)
        <- foldM epoch (ncs0,replicate 6 []). take nepch $ zip3 [0..] epss nstpzs

    let alss = transpose $ reverse alss0

    let lnlyt = execEC $ do

            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Log-Likelihood"

            sequence_ $ do

                (ttl,clr,ls) <- zip3 ttls pltclrs alss

                return . plot . liftEC $ do

                    plot_lines_title .= ttl
                    plot_lines_style .= solidLine 2 (opaque clr)
                    plot_lines_values .= [zip [(1 :: Int)..] ls]

    goalWriteFile sbdr (mlpfl ++ "ef1") . show $ listCoordinates (neuralNetwork ncrcef1')
    goalWriteFile sbdr (mlpfl ++ "cd1") . show $ listCoordinates (neuralNetwork ncrccd1')
    goalWriteFile sbdr (mlpfl ++ "ef2") . show $ listCoordinates (neuralNetwork ncrcef2')
    goalWriteFile sbdr (mlpfl ++ "cd2") . show $ listCoordinates (neuralNetwork ncrccd2')
    goalWriteFile sbdr (mlpfl ++ "ef3") . show $ listCoordinates (neuralNetwork ncrcef3')
    goalWriteFile sbdr (mlpfl ++ "cd3") . show $ listCoordinates (neuralNetwork ncrccd3')
    goalWriteFile sbdr "descent" $ show alss
    goalRenderableToSVG sbdr dstfl 800 400 $ toRenderable lnlyt


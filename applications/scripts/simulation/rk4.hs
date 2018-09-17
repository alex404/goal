--- Imports ---


-- Goal --

import Goal.Core

import Goal.Simulation

import qualified Data.Vector.Storable as C

--- Script ---


main = do

    -- Generation --

    -- We can simulate sin either as non-autonomous or second order autonomous
    let sin' t _ = C.singleton $ cos t
        vsin' x = C.fromList [x C.! 1,-x C.! 0]
        exp' = id

        t0 = 0
        tf = 10
        dt1 = 2
        dt2 = 1
        dt3 = 0.1
        ts1 = [t0,t0+dt1..tf-dt1]
        ts2 = [t0,t0+dt2..tf-dt2]
        ts3 = [t0,t0+dt3..tf-dt3]

        sx0 = C.fromList [0]
        vx0 = C.fromList [0,1]
        ex0 = C.fromList [1]

        sinFlow = nonAutonomousODE sin' t0 sx0
        sinFlowEuler = nonAutonomousODEEuler sin' t0 sx0
        vsinFlow = autonomousODE vsin' vx0
        expFlow = autonomousODE exp' ex0
        expFlowEuler = autonomousODEEuler exp' ex0

        ssimulator ts mly = zip ts $ (C.! 0) <$> stream ts mly
        vsimulator ts mly = zip ts $ (C.! 0) <$> stream ts mly
        esimulator ts mly = zip ts $ (C.! 0) <$> stream ts mly

    -- Plots --

    -- Sin

    let sinrnbl = toRenderable . execEC $ do

            layout_title .= "Sin Wave (dt = {1,0.1})"

            plot . liftEC $ do
                plot_lines_style .= dashedLine 3 [10,5] (opaque black)
                plot_lines_title .= "True"
                plot_lines_values .= [zip ts3 $ sin <$> ts3]

            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque blue)
                plot_lines_title .= "RK4"
                plot_lines_values .=
                  [ ssimulator ts2 sinFlow, ssimulator ts3 sinFlow ]

            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque red)
                plot_lines_title .= "Euler"
                plot_lines_values .=
                  [ ssimulator ts2 sinFlowEuler, ssimulator ts3 sinFlowEuler ]

            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque purple)
                plot_lines_title .= "RK4 2nd Order"
                plot_lines_values .=
                  [ vsimulator ts2 vsinFlow , vsimulator ts3 vsinFlow ]

    -- Exponential

    let exprnbl = toRenderable . execEC $ do

            layout_title .= "Exponential Function (dt = {2,1,0.1})"

            plot . liftEC $ do
                plot_lines_style .= dashedLine 3 [10,5] (opaque black)
                plot_lines_title .= "True"
                plot_lines_values .= [zip ts3 $ exp <$> ts3]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_title .= "RK4"
                plot_lines_values .= [ esimulator ts1 expFlow, esimulator ts2 expFlow, esimulator ts3 expFlow ]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_title .= "Euler"
                plot_lines_values .=
                    [ esimulator ts1 expFlowEuler, esimulator ts2 expFlowEuler, esimulator ts3 expFlowEuler ]

    -- IO

    goalRenderableToSVG "simulation" "rk4" 1200 600 . gridToRenderable . weights (1,1)
        . tallBeside sinrnbl $ tval exprnbl


--- Extra Functions ---


autonomousODEEuler f' x0 =
    accumulateFunction (0,x0) accumulator
      where accumulator t' (t,x) =
                let dt = t' - t
                    x' = x + stepEuler f' dt x
                 in (x',(t',x'))

nonAutonomousODEEuler f' t0 x0 =
    accumulateFunction (t0,x0) accumulator
      where accumulator t' (t,x) =
                let dt = t' - t
                    x' = x + stepEuler' f' t dt x
                 in (x',(t',x'))



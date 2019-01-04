{-# LANGUAGE FlexibleContexts,Arrows #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Simulation


--- Program ---


-- Globals --

-- Pendulum
m = 1
l = 1
pndl = Pendulum m l
fg = earthGravity
fd = Damping 1
f = (fg,fd)

-- Simulation
qdq0 = fromList (Bundle pndl) [3.14,0]
dt = 0.03
mxt = 20
ts = [0,dt..mxt]
flw = lagrangianFlow f qdq0
tivl = 2.0
fps = round $ recip dt

-- Plot
mnq = -pi
mndq = -5
mxq = pi
mxdq = 5

-- Functions --

vectorFieldLayout :: EC (Layout Double Double) ()
vectorFieldLayout = do

    layout_x_axis . laxis_generate .= scaledAxis def (mnq,mxq)
    layout_x_axis . laxis_title .= "Generalized Position"

    layout_y_axis . laxis_generate .= scaledAxis def (mndq,mxdq)
    layout_y_axis . laxis_title .= "Generalized Velocity"

pathToRenderable tqdqs =

    let phslyt = execEC $ do

            vectorFieldLayout

            layout_title .= "Phase Space"

            plot . liftEC $ do
                plot_lines_title .= "Phase"
                plot_lines_values .= [[(coordinate 0 qdq,coordinate 1 qdq) | (_,qdq) <- tqdqs]]
                plot_lines_style .= solidLine 3 (opaque black)

        enrlyt = execEC $ do

            layout_title .= "Energy"

            layout_x_axis . laxis_title .= "Time"

            layout_y_axis . laxis_generate .= scaledAxis def (0,20)
            layout_y_axis . laxis_title .= "Joules"

            let (tks,tus,tms) = unzip3 [ let (k,u) = mechanicalEnergy fg qdq in ((t,k),(t,u),(t,k+u)) | (t,qdq) <- tqdqs ]

            plot . liftEC $ do
                plot_lines_title .= "Kinetic"
                plot_lines_values .= [tks]
                plot_lines_style .= solidLine 3 (opaque red)

            plot . liftEC $ do
                plot_lines_title .= "Potential"
                plot_lines_values .= [tus]
                plot_lines_style .= solidLine 3 (opaque blue)

            plot . liftEC $ do
                plot_lines_title .= "Total"
                plot_lines_values .= [tms]
                plot_lines_style .= solidLine 3 (opaque purple)

        knmlyt = execEC $ do

            layout_title .= "Kinematics"

            let (_,qdq) = head tqdqs
                [tht] = listCoordinates $ position qdq
                (x,y) = (l * sin tht, -l * cos tht)

            layout_x_axis . laxis_generate .= scaledAxis def (-2,2)

            layout_y_axis . laxis_generate .= scaledAxis def (-2,2)

            plot . liftEC $ do
                plot_lines_values .= [[(0,0),(x,y)]]
                plot_lines_style .= solidLine 3 (opaque black)

            plot . liftEC $ do
                plot_points_values .= [(0,0)]
                plot_points_style .= hollowCircles 6 4 (opaque black)

            plot . liftEC $ do
                plot_points_values .= [(x,y)]
                plot_points_style .= filledCircles 6 (opaque black)


     in toRenderable . weights (1,1) . wideAbove enrlyt $ tval knmlyt .|. tval phslyt

-- Main --

main :: IO ()
main = do

    let mly = proc t -> do
            x <- flw -< t
            txs <- trajectoryWindow tivl -< (t,x)
            returnA -< pathToRenderable txs

    goalRenderablesToAnimation "simulation" "pendulum-simulation" fps 800 600 $ stream ts mly


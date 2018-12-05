{-# LANGUAGE Arrows #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Simulation
import Goal.Probability

import qualified Data.Vector.Storable as C
import qualified Numeric.LinearAlgebra.HMatrix as M


--- Program ---


-- Globals --

dt = 0.02
fps = round $ recip dt
mxt = 20
tivl = 2
t0 = 0
ts = [t0,dt..mxt]
x0 = C.fromList [0]
mu t _ = C.fromList [cos t]
sgma t _ = M.fromLists [[1 + cos (2*t)]]

trajectoryToRenderable ln = toRenderable . execEC $ do

    layout_y_axis . laxis_generate .= scaledAxis def (-6,6)

    plot . liftEC $ do
        plot_lines_title .= "Ito Process"
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [ln]

-- Main --

main :: IO ()
main = do

    flw <- runWithSystemRandom $ itoProcess mu sgma (t0,x0)

    let mly = proc t -> do
            x <- flw -< t
            pth <- trajectoryWindow tivl -< (t,x C.! 0)
            returnA -< trajectoryToRenderable pth

    goalRenderablesToAnimation "simulation" "ito-process"  fps 800 600 $ stream ts mly

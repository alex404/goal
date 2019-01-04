{-# LANGUAGE DataKinds,TypeOperators,FlexibleContexts #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Simulation

import qualified Goal.Core.Vector.Boxed as B

--- Program ---


-- Globals --

-- Pendulum
type Pendulum' = Pendulum (1/1) (1/1)

fg :: Gravity
fg = earthGravity
fd :: Damping
fd = Damping 1

-- Plot
mnq,mndq,mxq,mxdq :: Double
mnq = -pi
mndq = -5
mxq = pi
mxdq = 5

nrng :: Int
nrng = 10

rngq,rngdq :: [Double]
rngq = range mnq mxq nrng
rngdq = range mndq mxdq nrng

qdq0 :: Phase Pendulum' Double
qdq0 = Point $ B.doubleton 0 0

-- Functions --

vectorFieldLayout :: EC (Layout Double Double) ()
vectorFieldLayout = do

    goalLayout

    layout_x_axis . laxis_generate .= scaledAxis def (mnq,mxq)
    layout_x_axis . laxis_title .= "Generalized Position"

    layout_y_axis . laxis_generate .= scaledAxis def (mndq,mxdq)
    layout_y_axis . laxis_title .= "Generalized Velocity"

vectorFieldPlot :: AlphaColour Double -> EC (PlotVectors Double Double) ()
vectorFieldPlot clr = do

    plot_vectors_grid .= [ (x,y) | x <- rngq, y <- rngdq ]
    plot_vectors_style . vector_line_style . line_color .= clr
    plot_vectors_style . vector_head_style . point_color .= clr

pairToPendulum :: (x,x) -> Phase Pendulum' x
pairToPendulum (q,dq) = Point $ B.doubleton q dq

pairField :: (ForceField f Pendulum', RealFloat x) => f -> (x,x) -> (x,x)
pairField f qdq0' =
    let qdq = pairToPendulum qdq0'
     in B.toPair . coordinates $ vectorField f qdq

-- Main --

main :: IO ()
main = do

    let rnbl = toRenderable . execEC $ do

            vectorFieldLayout

            layout_title .= "Pendulum Vector Field"

            plot . fmap plotVectorField . liftEC $ do

                vectorFieldPlot $ opaque red
                plot_vectors_mapf .= pairField (fg,fd)
                plot_vectors_title .= "Damped"

            plot . fmap plotVectorField . liftEC $ do

                vectorFieldPlot $ opaque blue
                plot_vectors_mapf .= pairField fg
                plot_vectors_title .= "Conservative"

    goalRenderableToSVG "simulation" "pendulum-vector-field" 800 800 rnbl


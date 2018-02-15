{-# LANGUAGE Arrows #-}

-- | Provides tools for dynamic plots.
module Goal.Simulation.Plot
    ( -- * Plotting Mealys
      chainWindow
    , trajectoryWindow
    -- * Physics
    , sliceVectorField
    -- * Animation
    , goalRenderablesToAnimation
    ) where


-- Imports --


-- Goal --

import Goal.Core
import Goal.Geometry
--import Goal.Simulation.Mealy
import Goal.Simulation.Physics

--import System.Process

--- Processes ---


-- | When this 'Mealy' is fed a 'stream' of outputs from a 'Chain', it returns
-- lists of buffered output, which allows a window of outputs to be dynamically
-- plotted.
chainWindow
    :: Int -- ^ Buffer Size
    -> Mealy x [x] -- ^ Mealy
chainWindow n = accumulateMealy [] $ proc (x,xs) -> do
    let xs' = take n $ x:xs
    returnA -< (xs',xs')

-- | When this 'Mealy' is fed a 'stream' of outputs from a 'Flow, it returns
-- lists of buffered output, where the the "tail" of the output falls within the
-- given time interval.
trajectoryWindow
    :: Double -- ^ Time Interval
    -> Mealy (Double,x) [(Double,x)] -- ^ Mealy
trajectoryWindow tivl = accumulateMealy [] $ proc ((t,x),xts) -> do
    let xts' = takeWhile (\(t',_) -> t - t' < tivl) $ (t,x):xts
    returnA -< (xts',xts')

--- Physics ---


-- | Projects a vector field onto a single degree of freedom and expresses the
-- result in terms of the coordinates, allowing us to plot slices of high
-- dimensional vector fields.
sliceVectorField :: (Riemannian Generalized m, ForceField f m, RealFloat x)
    => Double -- ^ Scales the resulting vector field
    -> Int -- ^ Selects the degree of freedom
    -> f -- ^ The force field
    -> (Phase m x) -- ^ 'Anchor' configuration
    -> (Double,Double) -- ^ Input position, velocity
    -> (Double,Double) -- ^ Output (scaled) velocity, acceleration
sliceVectorField scl i f qdq =
    pointToPair scl i . vectorField f . pairToPoint i qdq

pairToPoint
    :: Manifold m
    => Int
    -> Point Directional (PhaseSpace m) x
    -> (Double, Double)
    -> Point Directional (PhaseSpace m) x
pairToPoint n qdq (x,dx) =
    let (hxs,_:txs) = splitAt n . listCoordinates $ position qdq
        (hdxs,_:tdxs) = splitAt n . listCoordinates $ detachTangentVector qdq
     in fromList (manifold qdq) $ hxs ++ x:txs ++ hdxs ++ dx:tdxs

pointToPair :: Manifold m => Double -> Int -> Point c m x -> (Double, Double)
pointToPair scl n dqddq =
     (scl * coordinate n dqddq, scl * coordinate (n + div 2 (dimension $ manifold dqddq)) dqddq)



--- Animations ---


goalRenderablesToAnimation
    :: String -- ^ Subdirectory
    -> String -- ^ Filename
    -> Int -- ^ Frames per second
    -> Int -- ^ x-pixels
    -> Int -- ^ y-pixels
    -> [Renderable a] -- ^ Renderables
    -> IO () -- ^ Write animation
goalRenderablesToAnimation sbdr flnm0 fps nx ny rnbls = do
    let sbdr' = sbdr ++ "/" ++ flnm0
    sbpth <- goalCreateSubdirectory sbdr'
    sequence_ [ goalRenderableToPNG sbdr' flnm nx ny rnbl | (flnm,rnbl) <- zip (show <$> [0..]) rnbls ]
    callProcess "ffmpeg"
        ["-y", "-r", show fps, "-s", show nx ++ "x" ++ show ny, "-i", sbpth ++ "/%d.png", sbpth ++ ".mp4"]
    goalRemoveSubdirectory sbdr'

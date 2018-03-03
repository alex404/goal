{-# LANGUAGE TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition

-- Qualified --

import qualified Data.Vector.Storable as C
import qualified Numeric.LinearAlgebra.HMatrix as M


--- Globals ---


nstps = 49

-- Simulation --

dt = 0.1
ts = [0,dt..]

-- Pendulum --

m = 1
l = 1
pndl = Pendulum m l

fg = earthGravity
Gravity fgx = fg
fdx = 0.1
fd = Damping fdx
f = (fg,fd)

mnq = -3
mxq = 3
mndq = -6
mxdq = 6

sgx = 1
sgma qdq = fromList (Tensor (Tangent $ velocity qdq) (Tangent $ velocity qdq)) [sgx]

sgy11 = 1
sgy22 = 2

pendulumTransition :: Coordinates -> RandST r Coordinates
pendulumTransition qdq =
    fmap coordinates . langevinTransition dt f sgma $ fromCoordinates (Bundle pndl) qdq

pendulumChain :: Coordinates -> RandST r (Chain Coordinates)
pendulumChain x0 = chain x0 pendulumTransition

randomPendulumState :: RandST s Coordinates
randomPendulumState = do
    q <- generate $ Standard # fromList (Uniform mnq mxq) []
    dq <- generate $ Standard # fromList (Uniform mndq mxdq) []
    return $ C.fromList [q,dq]

pendulumResponseMealy :: RandST r (Mealy Coordinates Coordinates)
pendulumResponseMealy = accumulateRandomFunction0 $
    \x -> generate $ Standard # fromList (MultivariateNormal 2) (C.toList x ++ [sgy11,0,0,sgy22])

-- Filter --

ekfProcessStep :: Coordinates -> Coordinates
ekfProcessStep cs =
    let ddq xs = coordinates . vectorField f $ fromCoordinates (Bundle pndl) xs
     in cs + stepRK4 ddq dt cs

ekfProcessJacobian :: Coordinates -> M.Matrix Double
ekfProcessJacobian cs =
    let [q,_] = C.toList cs
     in M.fromLists [[1,dt],[-dt*fgx/l*cos q, 1 - dt*fdx]]

ekfProcessCovariance :: M.Matrix Double
ekfProcessCovariance =
     M.fromLists [[0,0],[0,dt*sgx]]

ekfObservationStep :: Coordinates -> Coordinates
ekfObservationStep cs = cs

ekfObservationJacobian :: Coordinates -> M.Matrix Double
ekfObservationJacobian _ =
     M.fromLists [[1,0],[0,1]]

ekfObservationCovariance :: M.Matrix Double
ekfObservationCovariance =
     M.fromLists [[sgy11,0],[0,sgy22]]


-- Functions --

coordinateLayout :: Int -> [Coordinates] -> [Coordinates] -> [Standard :#: MultivariateNormal] -> Layout Double Double
coordinateLayout n xs ys ekfs = execEC $ do

    goalLayout

    layout_x_axis . laxis_title .= "Time"
    layout_y_axis . laxis_title .= if n == 0 then "Position" else "Velocity"

    plot . liftEC $ do

        plot_lines_values .= [zip ts ((C.! n) <$> xs)]
        plot_lines_style .= solidLine 3 (opaque black)

    plot . liftEC $ do

        plot_points_style .= filledCircles 3 (opaque black)
        plot_points_values .=
            zip ts ((C.! n) <$> ys)

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip ts (coordinate n <$> ekfs)]

--- Main ---


main :: IO ()
main = do

    x0 <- runWithSystemRandom randomPendulumState
    let p0 = Standard # fromList (MultivariateNormal 2) [x0 C.! 0, x0 C.! 1, 1, 0, 0, 1]

    xchn <- runWithSystemRandom (pendulumChain x0)
    nmly <- runWithSystemRandom pendulumResponseMealy
    let ekf = parametricFilter
                (extendedKalmanPrediction ekfProcessStep ekfProcessJacobian ekfProcessCovariance)
                (extendedKalmanInference ekfObservationStep ekfObservationJacobian ekfObservationCovariance) p0
        ekfchn = filterChain xchn nmly ekf
        (xs,ys,ekfs) = unzip3 . take nstps $ streamChain ekfchn
        lyt0 = coordinateLayout 0 xs ys ekfs
        lyt1 = coordinateLayout 1 xs ys ekfs

    void . goalRenderableToPDF "cognition/extended-kalman" "ekf-position" 500 200 $ toRenderable lyt0
    void . goalRenderableToPDF "cognition/extended-kalman" "ekf-velocity" 500 200 $ toRenderable lyt1

-- | Basic definitions for simulations of dynamical systems.
module Goal.Simulation.Flow (
    -- * Flows
      Flow
    -- ** Deterministic
    , autonomousODE
    , nonAutonomousODE
    , lagrangianFlow
    -- ** Stochastic
    , itoTransition
    , itoProcess
    , langevinTransition
    , langevinFlow
    -- * Integral Curves
    , stepEuler
    , stepRK4
    , stepEuler'
    , stepRK4'
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Simulation.Physics
import Goal.Simulation.Mealy

-- Qualified --

import qualified Numeric.LinearAlgebra.HMatrix as M


--- Flow ---


-- | A 'Flow' is a 'Mealy' automata with a certain semantics. In particular, 'stream'ing a 'Flow' takes an ordered set of times and returns the state of the process specified by the 'Flow' at each time.
type Flow x = Mealy Double x

-- Euclidean --

-- | Creates a 'Flow' which represents a deterministic, autonomous ordinary differential equation. Time is assumed to start at 0.
autonomousODE
    :: (Coordinates -> Coordinates) -- ^ Differential Equation
    -> Coordinates -- ^ Initial State
    -> Flow Coordinates -- ^ Differential Process
autonomousODE f' x0 =
    accumulateFunction (0,x0) accumulator
      where accumulator t' (t,x)
                | t' == t = (x,(t,x))
                | otherwise =
                      let dt = t' - t
                          x' = x + stepRK4 f' dt x
                       in (x',(t',x'))

-- | Creates a 'Flow' which represents a deterministic, non-autonomous ordinary differential equation.
nonAutonomousODE
    :: (Double -> Coordinates -> Coordinates) -- ^ Differential Equation
    -> Double -- ^ Initial Time
    -> Coordinates -- ^ Initial State
    -> Flow Coordinates -- ^ Differential Process
nonAutonomousODE f' t0 x0 =
    accumulateFunction (t0,x0) accumulator
      where accumulator t' (t,x)
                | t' == t = (x,(t,x))
                | otherwise =
                      let dt = t' - t
                          x' = x + stepRK4' f' t dt x
                       in (x',(t',x'))

-- Mechanical --

-- | Creates a 'Flow' which represents the evolution of a mechanical system. Time is assumed to start at 0.
lagrangianFlow :: (Riemannian Generalized m, ForceField f m)
    => f -- ^ Force field
    -> (Partials :#: PhaseSpace m) -- ^ Initial State
    -> Flow (Partials :#: PhaseSpace m) -- ^ Flow
lagrangianFlow f p0 =
    accumulateFunction (0,p0) accumulator
      where accumulator t' (t,qdq)
                | t' == t = (qdq,(t,qdq))
                | otherwise =
                      let dt = t' - t
                          ddq xs = coordinates . vectorField f $ fromCoordinates (manifold qdq) xs
                          qdq' = qdq <+> fromCoordinates (manifold qdq) (stepRK4 ddq dt $ coordinates qdq)
                       in (qdq', (t',qdq'))

--- Stochastic ---

-- | The transition function of a non-autonomous Ito process on a Euclidean space.
itoTransition
    :: Double -- ^ Time step size
    -> (Double -> Coordinates -> Coordinates) -- ^ The drift
    -> (Double -> Coordinates -> M.Matrix Double) -- ^ The diffusion
    -> (Double,Coordinates) -- ^ The current (time,state)
    -> RandST s Coordinates -- ^ The Ito transition
itoTransition dt mu sgma (t,x) = do
    let x' = x + stepRK4' mu t dt x
    generate $ Standard # joinMultivariateNormal x' (M.scale (sqrt dt) $ sgma t x)

-- | Returns an Ito 'Flow' on Euclidean space.
itoProcess
    :: (Double -> Coordinates -> Coordinates) -- ^ The drift
    -> (Double -> Coordinates -> M.Matrix Double) -- ^ The diffusion
    -> (Double,Coordinates) -- ^ The initial (time,state)
    -> RandST s (Flow Coordinates) -- ^ The Ito transition
itoProcess mu0 sgma0 (t0,x0) =
    accumulateRandomFunction (t0,x0) (accumulator mu0 sgma0)
      where accumulator mu sgma t' (t,x)
                | t' == t = return (x,(t,x))
                | otherwise = do
                      let dt = t' - t
                      x' <- itoTransition dt mu sgma (t,x)
                      return (x',(t',x'))

-- | The transition function for an autonomous, stochastic mechanical system.
langevinTransition :: (Riemannian Generalized m, ForceField f m)
    => Double -- ^ Time step
    -> f -- ^ The force field
    -> ( Partials :#: PhaseSpace m
       -> Function Differentials Differentials
       :#: Tensor (GeneralizedAcceleration m) (GeneralizedAcceleration m) ) -- ^ The second order noise function
    -> (Partials :#: PhaseSpace m) -- ^ The current state
    -> RandST s (Partials :#: PhaseSpace m) -- ^ The Langevin transition
langevinTransition dt f sgma qdq = do
    let flx p = matrixSquareRoot $ sgma p
        dq = bundleToTangent qdq
    nrms <- replicateM (dimension $ manifold dq) . generate $ Standard # fromList Normal [0,1]
    let lng p = sqrt dt /> (flx p >.> fromList (Tangent $ bundleToTangent p) nrms)
        ddq xs = coordinates . vectorField (f, lng) $ fromCoordinates (manifold qdq) xs
        qdq' = qdq <+> fromCoordinates (manifold qdq) (stepRK4 ddq dt $ coordinates qdq)
    return qdq'

-- | Returns a realization of a stochastic mechanical system. Time is assumed to start at 0.
langevinFlow :: (Riemannian Generalized m, ForceField f m)
    => f -- ^ The force field
    -> ( Partials :#: PhaseSpace m
       -> Function Differentials Differentials
       :#: Tensor (GeneralizedAcceleration m) (GeneralizedAcceleration m) ) -- ^ The second order noise function
    -> (Partials :#: PhaseSpace m) -- ^ The initial state
    -> RandST s (Flow (Partials :#: PhaseSpace m)) -- ^ The Ito transition
langevinFlow f0 sgma0 qdq0 =
    accumulateRandomFunction (0,qdq0) (accumulator f0 sgma0)
      where accumulator f sgma t' (t,qdq)
                | t' == t = return (qdq,(t,qdq))
                | otherwise = do
                      let dt = t' - t
                      qdq' <- langevinTransition dt f sgma qdq
                      return (qdq',(t',qdq'))

--- Integral Curves ---


-- | Returns the difference from an autonomous Euler step.
stepEuler
    :: (Coordinates -> Coordinates) -- ^ The derivative of the function to simulate
    -> Double -- ^ The time step 'dt'
    -> Coordinates -- ^ The state of the system at the current time
    -> Coordinates -- ^ The difference to the state of the system at time 't' + 'dt'
stepEuler f' dt x = realToFrac dt * f' x

-- | Returns the difference from an autonomous RK4 step.
stepRK4
    :: (Coordinates -> Coordinates) -- ^ The derivative of the function to simulate
    -> Double -- ^ The time step 'dt'
    -> Coordinates -- ^ The state of the system at the current time
    -> Coordinates -- ^ The difference to the state of the system at time 't' + 'dt'
stepRK4 f' dt x =
    let k1 = realToFrac dt * f' x
        k2 = realToFrac dt * f' (x + k1 / 2)
        k3 = realToFrac dt * f' (x + k2 / 2)
        k4 = realToFrac dt * f' (x + k3)
    in (k1 + 2 * k2 + 2 * k3 + k4) / 6

-- | Returns the difference from a non-autonomous Euler step.
stepEuler'
    :: (Double -> Coordinates -> Coordinates) -- ^ The derivative of the function to simulate
    -> Double -- ^ The current time 't'
    -> Double -- ^ The time step 'dt'
    -> Coordinates -- ^ The state of the system at the current time
    -> Coordinates -- ^ The difference to the state of the system at time 't' + 'dt'
stepEuler' f' t dt x = realToFrac dt * f' t x

-- | Returns the difference from a non-autonomous RK4 step.
stepRK4'
    :: (Double -> Coordinates -> Coordinates) -- ^ The derivative of the function to simulate
    -> Double -- ^ The current time 't'
    -> Double -- ^ The time step 'dt'
    -> Coordinates -- ^ The state of the system at the current time
    -> Coordinates -- ^ The difference to the state of the system at time 't' + 'dt'
stepRK4' f' t dt x =
    let k1 = realToFrac dt * f' t x
        k2 = realToFrac dt * f' (t + dt / 2) (x + k1 / 2)
        k3 = realToFrac dt * f' (t + dt / 2) (x + k2 / 2)
        k4 = realToFrac dt * f' (t + dt) (x + k3)
    in (k1 + 2 * k2 + 2 * k3 + k4) / 6

{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts #-}

module Pendulum where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Simulation
import Goal.Probability
import Goal.Cognition

import qualified Data.Vector.Storable as C
import qualified Numeric.LinearAlgebra.HMatrix as H

--- Globals ---


-- Directory --

sbdr = "pendulum"
ttls = ["Observation","EKF","KF","RCR-EF","RCR-CD","ASM-EF","ASM-CD","NIV-EF","NIV-CD"]
pltclrs = [black,green,blue,darkblue,red,darkred,violet,darkviolet]
pltlns = [sld black, sld green, sld' red, dshd red, sld' blue, dshd blue, sld' purple, dshd purple]
        where sld = solidLine 2 . opaque
              sld' = solidLine 2 . (`withOpacity` 0.6)
              dshd = dashedLine 2 [4,4] . opaque


-- Simulation --


-- Simulation --

dt = 0.02
ts = [0,dt..]

-- Pendulum --

m = 1
l = 1
pndl = Pendulum m l
fg = earthGravity
Gravity fgx = fg
fdx = 0.1
fd = Damping fdx

sgmans = 1
sgma qdq = fromList (Tensor (Tangent $ velocity qdq) (Tangent $ velocity qdq)) [sgmans]
mnq = -pi
mndq = -12
mxq = pi
mxdq = 12

thrsh = pi

f = (fg,fd)

-- Population Codes --

nnq = 10
qprsn = recip 2
rngq = tail $ range mnq mxq (nnq+1)
tcsq = [ fromList VonMises [q,qprsn] | q <- rngq ]

nndq = 10
dqvr = 4
rngdq = range mndq mxdq nndq
tcsdq = [ fromList Normal [dq,dqvr] | dq <- rngdq ]

gnq = dt * 100
gndq = dt * 100
--emsn = vonMisesPopulationEncoder sps gn
emsn = vonMisesNormalPopulationEncoder' tcsq tcsdq gnq gndq
dcdn = matrixTranspose . snd $ splitAffine emsn
Affine mn _ = manifold emsn

-- Filtering Population --

mz = Replicated Poisson 20
thtz = 10

(dcdz1,amtx1) = orthogonalCode mz thtz emsn
dcdz2 = dcdz1
dcdz3 = dcdn
amtx2 = amtx1
amtx3 = matrixIdentity mn

-- Prediction Population --

dcdy1 = dcdz1
dcdy2 = dcdn
dcdy3 = dcdn

bmtx1 = Function Mixture Mixture # matrixIdentity mz
bmtx2 = amtx1
bmtx3 = Function Mixture Mixture # matrixIdentity mn

y01 = Mixture # zero mz
y02 = Mixture # zero mn
y03 = Mixture # zero mn

-- Neural Network --

mg1 = ThreeLayerPerceptron mz (Replicated Bernoulli 500) mz
mg2 = ThreeLayerPerceptron mn (Replicated Bernoulli 500) mz
mg3 = ThreeLayerPerceptron mn (Replicated Bernoulli 500) mn
--mg3 = NeuralField (Replicated Bernoulli 500) mn dt


--- Neural Networks ---

data NeuralCircuits
    = NC1 ( HarmoniumCircuit (ThreeLayerPerceptron (Replicated Poisson) (Replicated Bernoulli) (Replicated Poisson))
            (Replicated Poisson) (VonMises, Normal) )
    | NC2 ( HarmoniumCircuit (NeuralField (Replicated Bernoulli) (Replicated Poisson))
            (Replicated Poisson) (VonMises, Normal) )

--- Neural Circuits ---

neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC1 ncrc) = (NC1 ^<<) <$> harmoniumCircuitTrainer eps bt1 bt2 rg nstpz ncrc
neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC2 ncrc) = (NC2 ^<<) <$> harmoniumCircuitTrainer eps bt1 bt2 rg nstpz ncrc

neuralCircuitFilter (NC1 ncrc) = harmoniumCircuitFilter ncrc
neuralCircuitFilter (NC2 ncrc) = harmoniumCircuitFilter ncrc

increment (NC1 ncrc) = NC1 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }
increment (NC2 ncrc) = NC2 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }

-- Kalman Filtering --

ek0 = transition $ Standard # fromList (VonMises,Normal) [0.1,0.1,0.1,10]
--ek0' = transition $ Standard # fromList (MultivariateNormal 2) [0.1,0.1,10,0,0,10]
ekflt = parametricFilter approximateExtendedKalmanPrediction (rectifiedBayesRule emsn) ek0
--ekflt' = parametricFilter
    --(transition . extendedKalmanPrediction ekfProcessStep ekfProcessJacobian ekfProcessCovariance . transition) approximateKalmanInference ek0'
kflt = parametricFilter approximateKalmanPrediction (rectifiedBayesRule emsn) ek0

kfProcessStep :: Coordinates -> Coordinates
kfProcessStep cs =
    let [q,dq] = C.toList cs
     in cs + C.fromList [dt * dq, -dt*(m*fgx*l*q + fdx*dq)]

kfProcessJacobian :: H.Matrix Double
kfProcessJacobian =
     H.fromLists [[1,dt],[-dt*fgx*l*m, 1 - dt*fdx]]

ekfProcessStep :: Coordinates -> Coordinates
ekfProcessStep cs =
    let ddq xs = coordinates . vectorField f $ fromCoordinates (Bundle pndl) xs
     in cs + stepRK4 ddq dt cs

ekfProcessJacobian :: Coordinates -> H.Matrix Double
ekfProcessJacobian cs =
    let [q,_] = C.toList cs
     in H.fromLists [[1,dt],[-dt*fgx/l*cos q, 1 - dt*fdx]]

ekfProcessCovariance :: H.Matrix Double
ekfProcessCovariance =
     H.fromLists [[0,0],[0,dt*sgmans]]

approximateKalmanPrediction :: Natural :#: (VonMises, Normal) -> Natural :#: (VonMises, Normal)
approximateKalmanPrediction p0 =
    let [q,qp,dq,dqsd] = listCoordinates $ transitionTo Standard p0
        p = fromList (MultivariateNormal 2) [q,dq,recip qp,0,0,dqsd]
        [q',dq',qsd',_,_,dqsd'] = listCoordinates $ kalmanPrediction kfProcessJacobian ekfProcessCovariance p
     in transitionTo Natural $ Standard # fromList (VonMises,Normal) [q',recip qsd', dq', dqsd']

approximateExtendedKalmanPrediction :: Natural :#: (VonMises, Normal) -> Natural :#: (VonMises, Normal)
approximateExtendedKalmanPrediction p0 =
    let [q,qp,dq,dqsd] = listCoordinates $ transitionTo Standard p0
        p = fromList (MultivariateNormal 2) [q,dq,recip qp,0,0,dqsd]
        [q',dq',qsd',_,_,dqsd'] = listCoordinates $ extendedKalmanPrediction ekfProcessStep ekfProcessJacobian ekfProcessCovariance p
     in transitionTo Natural $ Standard # fromList (VonMises,Normal) [q',recip qsd', dq', dqsd']

{-
approximateKalmanInference :: [Int] -> Natural :#: MultivariateNormal -> Natural :#: MultivariateNormal
approximateKalmanInference n p0 =
    let [qn,qpn,dqn,dqvrn] = listCoordinates . transitionTo Standard $ dcdn >.>* n
        q0 = coordinate 0 $ transitionTo Natural p0
        mvnp = transitionTo Natural $ Standard # fromList (MultivariateNormal 2) [qn + 2*pi*fromIntegral (revolutions q0),dqn,1/qpn,0,0,dqvrn]
     in p0 <+> mvnp
     -}

--- Functions ---

squaredError :: Natural :#: (VonMises, Normal) -> (Double,Double) -> Double
squaredError p (q,dq) =
    let [qht,_,dqht,_] = listCoordinates $ transitionTo Standard p
        q' = q - fromIntegral (revolutions q) * 2 * pi
     in (minimum [(q'-qht + pi)^2,(q'-qht)^2,(q'-qht - pi)^2] + (dq-dqht)^2)
     --in (q'-qht)^2

{-
squaredError' :: Natural :#: MultivariateNormal -> (Double,Double) -> Double
squaredError' p (q,dq) =
    let (qht:dqht:_) = listCoordinates $ transitionTo Standard p
        q' = q - fromIntegral (revolutions q) * 2 * pi
     in (minimum [(q'-qht + pi)^2,(q'-qht)^2,(q'-qht - pi)^2] + (dq-dqht)^2)
     --in (q'-qht)^2
     -}


periodicBreaker :: [(Double,Double)] -> [[(Double,Double)]]
periodicBreaker [] = []
periodicBreaker txs =
    let (txs0,txs') = span (\((_,x),(_,x')) -> abs (x - x') < 2) . zip txs $ tail txs
        txs0' = if null txs0 then [] else fst (head txs0):(snd <$> txs0)
     in txs0' : periodicBreaker (snd <$> txs')

pendulumTransition :: (Double,Double) -> RandST r (Double,Double)
pendulumTransition (q,dq) = do
    qdq' <- langevinTransition dt f sgma $ fromList (Bundle pndl) [q,dq]
    let [q',dq'] = listCoordinates $ periodic [True,False] qdq'
    return (q',dq')

pendulumChain :: (Double,Double) -> RandST r (Chain (Double,Double))
pendulumChain x0 = chain x0 pendulumTransition

randomPendulumState :: RandST s (Double,Double)
randomPendulumState = do
    q <- generate $ Standard # fromList (Uniform mnq mxq) []
    dq <- generate $ Standard # fromList (Uniform mndq mxdq) []
    return (q,dq)

pendulumResponseMealy :: RandST r (Mealy (Double,Double) [Int])
pendulumResponseMealy = accumulateRandomFunction0 $
    \x -> standardGenerate $ emsn >.>* x


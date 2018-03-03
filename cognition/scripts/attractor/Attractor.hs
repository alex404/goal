{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts #-}

module Attractor where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition


--- Globals ---


sbdr = "attractor"
ttls = ["Observation","True","RCR-EF","RCR-CD","ASM-EF","ASM-CD","NIV-EF","NIV-CD"]
pltclrs = [black,green,blue,darkblue,red,darkred,violet,darkviolet]
pltlns = [sld black, sld green, sld' red, dshd red, sld' blue, dshd blue, sld' purple, dshd purple]
        where sld = solidLine 2 . opaque
              sld' = solidLine 2 . (`withOpacity` 0.6)
              dshd = dashedLine 2 [4,4] . opaque


-- Simulation --

dt = 0.02
ts = [0,dt..]

-- 1-D Linear Dynamical System --

a = -1
b = 1
b2 = b^2

-- Discrete-Time
xdnm = 1+dt*a
b2stp = b2*dt

-- Population Codes --

mnx = -6
mxx = 6

-- Observation Population
nn = 10
vrn = 2
mnn = mnx
mxn = mxx
rngn = range mnn mxn nn
sps = [ fromList Normal [x,vrn] | x <- rngn ]
gn = 100 * dt
emsn = normalPopulationEncoder sps gn
Affine mn _ = manifold emsn
dcdn = matrixTranspose (snd $ splitAffine emsn)

-- Filtering Population --

mz = Replicated Poisson 10
thtz = 10
(dcdz1,amtx1) = orthogonalCode mz thtz emsn
amtx2 = amtx1
amtx3 = matrixIdentity mn
dcdz2 = dcdz1
dcdz3 = dcdn

-- Prediction Population --

dcdy1 = dcdz1
dcdy2 = dcdn
dcdy3 = dcdn

bmtx1 = Function Mixture Mixture # matrixIdentity mz
bmtx2 = amtx1
bmtx3 = Function Mixture Mixture # matrixIdentity mn

y01 = Mixture # zero mz
y02 = Mixture # zero mn
y03 = y02

-- Neural Network --

--mg1 = NeuralField (Replicated Bernoulli 500) mz dt
mg1 = ThreeLayerPerceptron mz (Replicated Bernoulli 200) mz
mg2 = ThreeLayerPerceptron mn (Replicated Bernoulli 200) mz
mg3 = ThreeLayerPerceptron mn (Replicated Bernoulli 200) mn
--mg3 = NeuralField (Replicated Bernoulli 500) mn dt


--- Kalman Filter ---

kf0 = transition $ Standard # fromList Normal [0,10]
kflt = parametricFilter (transition . kalmanPrediction1D xdnm b2stp . transition) (rectifiedBayesRule emsn (zero Normal)) kf0

attractorTransition :: Double -> RandST r Double
attractorTransition x = generate $ Standard # fromList Normal [xdnm*x, b2stp]

attractorChain :: Double -> RandST r (Chain Double)
attractorChain x0 = chain x0 attractorTransition

attractorResponseMealy :: RandST r (Mealy Double [Int])
attractorResponseMealy = accumulateRandomFunction0 $
    \x -> standardGenerate $ emsn >.>* x

randomAttractorState :: RandST s Double
randomAttractorState = generate $ Standard # fromList (Uniform (mnx/2) (mxx/2)) []


--- Neural Circuits ---

data NeuralCircuits
    = NC1 ( HarmoniumCircuit (ThreeLayerPerceptron (Replicated Poisson) (Replicated Bernoulli) (Replicated Poisson))
            (Replicated Poisson) Normal )
    | NC2 ( HarmoniumCircuit (NeuralField (Replicated Bernoulli) (Replicated Poisson))
            (Replicated Poisson) Normal )
    | NC3 ( HarmoniumCircuit (NeuralField (Replicated Bernoulli) (Replicated MeanNormal))
            (Replicated Poisson) Normal )

neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC1 ncrc) = (NC1 ^<<) <$> harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz ncrc
neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC2 ncrc) = (NC2 ^<<) <$> harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz ncrc
neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC3 ncrc) = (NC3 ^<<) <$> harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz ncrc

neuralCircuitFilter (NC1 ncrc) = harmoniumCircuitFilter ncrc
neuralCircuitFilter (NC2 ncrc) = harmoniumCircuitFilter ncrc
neuralCircuitFilter (NC3 ncrc) = harmoniumCircuitFilter ncrc

increment (NC1 ncrc) = NC1 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }
increment (NC2 ncrc) = NC2 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }
increment (NC3 ncrc) = NC3 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }

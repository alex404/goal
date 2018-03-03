{-# LANGUAGE TypeOperators,TypeFamilies #-}

module NeuralHMM where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition


--- Markov Model ---

sbdr = "neural-hmm"
ttls = ["Observation","True","RCR-EF","RCR-CD","ASM-EF","ASM-CD","NIV-EF","NIV-CD"]
pltclrs = [black,green,blue,darkblue,red,darkred,violet,darkviolet]
pltlns = [sld black, sld green, sld' red, dshd red, sld' blue, dshd blue, sld' purple, dshd purple]
        where sld = solidLine 2 . opaque
              sld' = solidLine 2 . (`withOpacity` 0.6)
              dshd = dashedLine 2 [8,8] . opaque

-- Types --

data Colours = Red | Green | Blue deriving (Eq, Show)

xcat = [Red, Green, Blue]
mx = Categorical xcat

-- Functions --

trns :: Colours -> Standard :#: Categorical [Colours]
trns Red = fromList mx [0.8,0.15]
trns Green = fromList mx [0.25,0.5]
trns Blue = fromList mx [0.05,0.15]

clrs :: Colours -> Colour Double
clrs Red = red
clrs Green = green
clrs Blue = blue

optpos :: Colours -> Int
optpos Red = 2
optpos Green = 3
optpos Blue = 4

estpos :: Colours -> Int
estpos Red = -4
estpos Green = -3
estpos Blue = -2


--- Neural Circuit ---


-- Observation Population --

gn = 0.4
shft = -5

nn = 10
mn = Replicated Poisson nn

bs = (+shft) . (*gn) . fromIntegral <$> [0..nn-1]
mtx1 = zipWith subtract bs $ reverse bs
mtx2 = zipWith subtract bs . replicate nn . log . average $ exp <$> bs

b = fromList (Replicated Poisson nn) bs
mtx = fromList (Tensor mn mx) . concat . transpose $ [mtx1,mtx2]

emsn :: Function Mixture Natural :#: Affine (Replicated Poisson) (Categorical [Colours])
emsn = joinAffine b mtx

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
y03 = Mixture # zero mn

-- Neural Network --

mg1 = ThreeLayerPerceptron mz (Replicated Bernoulli 100) mz
mg2 = ThreeLayerPerceptron mn (Replicated Bernoulli 100) mz
mg3 = ThreeLayerPerceptron mn (Replicated Bernoulli 100) mn


--- Functions ---


markovChain :: Colours -> RandST s (Chain Colours)
markovChain x0 = chain x0 (generate . trns)

responseMealy :: RandST s (Mealy Colours [Int])
responseMealy = accumulateRandomFunction0 (standardGenerate . (emsn >.>*))


--- Neural Networks ---

data NeuralCircuits
    = NC1 ( HarmoniumCircuit (ThreeLayerPerceptron (Replicated Poisson) (Replicated Bernoulli) (Replicated Poisson))
            (Replicated Poisson) (Categorical [Colours]) )

neuralCircuitTrainer eps bt1 bt2 rg nstpz (NC1 ncrc) = (NC1 ^<<) <$> harmoniumCircuitTrainer' eps bt1 bt2 rg nstpz ncrc
neuralCircuitFilter (NC1 ncrc) = harmoniumCircuitFilter ncrc

increment (NC1 ncrc) = NC1 $ ncrc { maybeCDN = (+1) <$> maybeCDN ncrc }

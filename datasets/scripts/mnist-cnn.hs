{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Datasets.MNIST


-- Globals ---


-- Network --

type NeuralNetwork = Categorical Int 10 <+< Convolutional 100 2 1 MNISTHeight MNISTWidth 1 Bernoulli Poisson

-- Data --


-- Training --

nepchs :: Int
nepchs = 1000

eps :: Double
eps = -0.05

-- Momentum
mxmu :: Double
mxmu = 0.999

mu :: Int -> Double
mu = defaultMomentumSchedule mxmu

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8


-- Main --


main :: IO ()
main = undefined

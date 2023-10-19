#!/usr/bin/env cabal
{- cabal:
build-depends: base
-}
{- project:
packages: /home/alex404/code/goal/core
-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}

--- Imports ---

import Goal.Core
import Goal.Geometry

--- Globals ---

-- Functions --

f :: Point # Euclidean 2 -> Double
f p =
    let [x, y] = listCoordinates p
     in square x + square y + square (x - y)

-- Gradient Descent --

p0 :: Cartesian # Euclidean 2
p0 = fromTuple (-4, 2)

bnd :: Double
bnd = 0.0001

cauchify :: [Cartesian # Euclidean 2] -> [Cartesian # Euclidean 2]
cauchify = cauchySequence euclideanDistance bnd

eps, mtm :: Double
eps = -0.05
mtm = 0.9

path :: GradientPursuit -> [Cartesian # Euclidean 2]
path gp = cauchify $ gradientSequence (differential f) eps gp p0

grds, mtms, adms :: [Cartesian # Euclidean 2]
grds = path Classic
mtms = path $ defaultMomentumPursuit mtm
adms = path defaultAdamPursuit

-- Plot --

ldpth :: String
ldpth = "."

rng :: [Double]
rng = range (-4) 4 100

isosmps :: [(Double, Double, Double)]
isosmps = do
    x <- rng
    y <- rng
    return (x, y, f $ B.fromTuple (x, y))

isonm, grdnm, mtmnm, admnm :: String
isonm = "isosamples"
grdnm = "gradient-descent"
mtmnm = "momentum"
admnm = "adam"

--- Main ---

main :: IO ()
main = do
    putStrLn "Gradient Descent Steps:"
    print $ length grds - 1

    putStrLn "Momentum Steps:"
    print $ length mtms - 1

    putStrLn "Adam Steps:"
    print $ length adms - 1

    goalExport ldpth isonm isosmps
    goalExport ldpth grdnm $ listCoordinates <$> grds
    goalExport ldpth mtmnm $ listCoordinates <$> mtms
    goalExport ldpth admnm $ listCoordinates <$> adms

    runGnuplot ldpth "gradient-descent"

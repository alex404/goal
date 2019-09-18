{-# LANGUAGE DataKinds,TypeOperators #-}


--- Imports ---


import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

--- Globals ---


-- Functions --

f :: RealFrac x => B.Vector 2 x -> x
f xs =
    let (x,y) = G.toPair xs
        two :: Int
        two = 2
     in x^two + y^two + (x-y)^two

-- Plot --

niso :: Int
niso = 10

cntrf :: Double -> Double -> Double
cntrf x y = f $ G.doubleton x y

rng :: [Double]
rng = range (-4) 4 400

-- Gradient Descent --

p0 :: Cartesian # Euclidean 2
p0 = Point $ G.doubleton (-4) 2

bnd,eps :: Double
bnd = 0.0001
eps = -0.05

cauchify :: [Cartesian # Euclidean 2] -> [Cartesian # Euclidean 2]
cauchify = cauchySequence euclideanDistance bnd

grds,mtms,adms :: [Cartesian # Euclidean 2]
grds = cauchify $ gradientSequence (differential f) eps Classic p0
mtms = cauchify $ gradientSequence (differential f) eps (defaultMomentumPursuit 0.9) p0
adms = cauchify $ gradientSequence (differential f) eps defaultAdamPursuit p0

-- Plot --

ldpth :: String
ldpth = "."

isosmps :: [(Double, Double, Double)]
isosmps = do
    x <- rng
    y <- rng
    return (x,y,f $ B.doubleton x y)

isonm :: String
isonm = "isosamples"

grdnm :: String
grdnm = "gradient-ascent"

mtmnm :: String
mtmnm = "momentum"

admnm :: String
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

    goalExport ldpth grdnm $ listCoordinates <$> grds
    goalExport ldpth mtmnm $ listCoordinates <$> mtms
    goalExport ldpth admnm $ listCoordinates <$> adms

    runGnuplot ldpth "gradient-ascent"


{-# LANGUAGE TypeOperators, TypeFamilies, FlexibleContexts #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

--- Globals ---

-- True Normal --

sp1 :: Source # Normal
sp1 = Point $ S.doubleton 2 3

-- Gradient Descent --

meps,neps,geps :: Double
meps = -0.1
neps = -0.01
geps = -0.1

mbnd,nbnd,gbnd :: Double
mbnd = 1e-10
nbnd = 1e-10
gbnd = 1e-10

sp0 :: Source # Normal
sp0 = Point $ S.doubleton 0.5 1.5

np0 :: Natural # Normal
np0 = transition sp0

mp0 :: Mean # Normal
mp0 = transition sp0

-- Plot --

mnmu,mxmu,mnvr,mxvr :: Double
mnmu = 0
mxmu = 4
mnvr = 1
mxvr = 5

murng,vrrng :: (Double,Double,Int)
murng = (mnmu,mxmu,1000)
vrrng = (mnvr,mxvr,1000)

niso :: Int
niso = 20

naturalDifferentials
     :: Point Natural Normal
     -> CotangentPair Natural Normal
naturalDifferentials q = joinTangentPair q $ crossEntropyDifferential (transition sp1) q

-- Functions --

-- Layout --

main :: IO ()
main = do

    --let mps = cauchySequence relativeEntropy mbnd $ vanillaGradientSequence meps meanDifferentials mp0
        --nmps = take stps $ gradientSequence meps mixtureDifferentials mp0

    let nps = cauchySequence (transition2 relativeEntropy) nbnd
            $ vanillaGradientSequence naturalDifferentials neps Classic np0
        gps = cauchySequence (transition2 relativeEntropy) gbnd
            $ gradientSequence naturalDifferentials geps Classic np0

    --putStrLn "Mean Coordinate Descent Steps:"
    --print $ length mps

    putStrLn "Natural Coordinate Descent Steps:"
    print $ length nps

    putStrLn "Geometric Gradient Descent Steps:"
    print $ length gps

    goalExport ldpth grdnm $ listCoordinates <$> grds
    goalExport ldpth mtmnm $ listCoordinates <$> mtms
    goalExport ldpth admnm $ listCoordinates <$> adms

    runGnuplot ldpth "coordinate-ascent"

#! stack runghc

{-# LANGUAGE GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Data.List as L
import qualified Goal.Core.Vector.Storable as S


--- Globals ---


-- Prior

psrx :: Source # Normal
psrx = fromTuple (0,10)

prx :: Natural # Normal
prx = toNatural psrx

-- Emission Distribution

escl,evr,eshft :: Double
escl = 1
evr = 4
eshft = 1

enzx :: Natural # Tensor NormalMean NormalMean
enzx = singleton $ escl / evr

enrz :: Natural # Normal
enrz = fromTuple (eshft,-1/(2*evr))

efzx :: Natural # Affine Tensor NormalMean Normal NormalMean
efzx = join enrz enzx

-- Transition Distribution

tscl,tvr,tshft :: Double
tscl = 0.5
tvr = 4
tshft = 2

tnzx :: Natural # Tensor NormalMean NormalMean
tnzx = singleton $ tscl / tvr

tnrz :: Natural # Normal
tnrz = fromTuple (tshft,-1/(2*tvr))

tfzx :: Natural # Affine Tensor NormalMean Normal NormalMean
tfzx = join tnrz tnzx

-- Latent Process

ltnt :: Natural # LatentProcess Tensor Tensor NormalMean NormalMean Normal Normal
ltnt = joinLatentProcess prx efzx tfzx

-- Conjugation Curve Plotting

xsmps :: [Double]
xsmps = range (-10) 10 1000

ys :: [Double]
ys = potential <$> efzx >$>* xsmps

sx :: Double -> S.Vector 3 Double
sx x = S.fromTuple (1,x, x**2)

rho0 :: Double
rprms :: Natural # Normal
(rho0,rprms) = conjugationParameters efzx

bts :: S.Vector 3 Double
bts =
    let (rho1,rho2) = S.toPair $ coordinates rprms
     in S.fromTuple (rho0,rho1,rho2)

yhts :: [Double]
yhts = S.dotMap bts $ sx <$> xsmps

-- Simulation

nstps :: Int
nstps = 20


--- Main ---


main :: IO ()
main = do

    zxpth <- realize $ sampleLatentProcess nstps ltnt

    let (zpth,xpth) = unzip zxpth
        flts = conjugatedFiltering tfzx efzx prx zpth
        (mus,vrs) = unzip $ S.toPair . coordinates . toSource <$> flts
        sds = sqrt <$> vrs

    goalExport "." "conjugation-curve" $ zip3 xsmps ys yhts
    runGnuplot "." "conjugation-curve"

    goalExport "." "kalman-filter" $ L.zip5 [0 :: Int ..] xpth zpth mus sds
    runGnuplot "." "kalman-filter"

    putStrLn "Average SD:"
    print $ average sds

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
psrx = fromTuple (0,2)

prx :: Natural # Normal
prx = toNatural psrx

-- Emission Distribution

enzx :: Natural # Tensor NormalMean NormalMean
enzx = 1

esrz :: Source # Normal
esrz = fromTuple (0,1)

efzx :: Natural # Affine Tensor NormalMean Normal NormalMean
efzx = join (toNatural esrz) enzx

-- Transition Distribution

tnzx :: Natural # Tensor NormalMean NormalMean
tnzx = 1

tsrz :: Source # Normal
tsrz = fromTuple (0,1)

tfzx :: Natural # Affine Tensor NormalMean Normal NormalMean
tfzx = join (toNatural tsrz) tnzx

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
        flts = conjugatedFiltering efzx tfzx prx zpth
        (mus,vrs) = unzip $ S.toPair . coordinates . toSource <$> flts
        sds = sqrt <$> vrs

    goalExport "." "conjugation-curve" $ zip3 xsmps ys yhts
    runGnuplot "." "conjugation-curve"

    goalExport "." "kalman-filter" $ L.zip5 [0 :: Int ..] xpth zpth mus sds
    runGnuplot "." "kalman-filter"

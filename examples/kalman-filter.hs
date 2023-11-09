{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Control.Monad (replicateM)
import Data.List qualified as DL
import Numeric (showFFloat)

--- Globals ---

-- Prior

prx, prx0 :: Source # Normal
prx = fromTuple (0, 1)
prx0 = fromTuple (0, 10)

nprx, nprx0 :: Natural # Normal
nprx = transition prx
nprx0 = transition prx0

-- Emission Distribution

enrm, enrm0 :: Source # Normal
enrm = fromTuple (0, 2)
enrm0 = fromTuple (0, 10)

escl, escl0 :: Source # Tensor (StandardNormal 1) (StandardNormal 1)
escl = 1
escl0 = 1

efzx, efzx0 :: Source # LinearModel L.PositiveDefinite 1 1
efzx = join enrm escl
efzx0 = join enrm0 escl0

nefzx, nefzx0 :: Natural # LinearModel L.PositiveDefinite 1 1
nefzx = transition efzx
nefzx0 = transition efzx0

-- Transition Distribution

tnrm, tnrm0 :: Source # Normal
tnrm = fromTuple (2, 4)
tnrm0 = fromTuple (0, 10)

tscl, tscl0 :: Source # Tensor (StandardNormal 1) (StandardNormal 1)
tscl = 0.5
tscl0 = 1

tfzx, tfzx0 :: Source # LinearModel L.PositiveDefinite 1 1
tfzx = join tnrm tscl
tfzx0 = join tnrm0 tscl0

ntfzx, ntfzx0 :: Natural # LinearModel L.PositiveDefinite 1 1
ntfzx = transition tfzx
ntfzx0 = transition tfzx0

-- Latent Process

kflt, kflt0 :: Natural # KalmanFilter 1 1
kflt = joinLatentProcess nprx nefzx ntfzx
kflt0 = joinLatentProcess nprx0 nefzx0 ntfzx0

-- Conjugation Curve Plotting

xsmps :: [S.Vector 1 Double]
xsmps = range (-10) 10 1000

ys :: [Double]
ys = potential <$> nefzx >$>* xsmps

sx :: Double -> S.Vector 3 Double
sx x = S.fromTuple (1, x, x ** 2)

rho0 :: Double
rprms :: Natural # Normal
(rho0, rprms) = conjugationParameters nefzx

bts :: S.Vector 3 Double
bts =
    let (rho1, rho2) = S.toPair $ coordinates rprms
     in S.fromTuple (rho0, rho1, rho2)

yhts :: [Double]
yhts = S.dotMap bts $ sx . S.head <$> xsmps

-- Simulation

nstps :: Int
nstps = 20

nepchs :: Int
nepchs = 50

nsmps :: Int
nsmps = 250

--- Main ---

main :: IO ()
main = do
    xzpth <- realize $ sampleLatentProcess nstps kflt

    let xpth, zpth :: [S.Vector 1 Double]
        (xpth, zpth) = unzip xzpth

        flts = conjugatedFiltering kflt xpth
        (mus, vrs) = unzip $ S.toPair . coordinates . toSource <$> flts
        sds = sqrt <$> vrs

    xzpths <- realize . replicateM nsmps $ sampleLatentProcess nstps kflt

    let xpths = map fst <$> xzpths

    let em = latentProcessExpectationMaximization xpths

        kflts = take nepchs $ iterate em kflt0
        kflt1 = last kflts

    let llsf kflt' = average $ logObservableDensities kflt' xpths
        lls = zip [0 :: Int ..] $ llsf <$> kflts

    mapM_ print lls

    let flts0 = conjugatedFiltering kflt0 zpth
        (mus0, vrs0) = unzip $ S.toPair . coordinates . toSource <$> flts0
        sds0 = sqrt <$> vrs0

    let flts1 = conjugatedFiltering kflt1 zpth
        (mus1, vrs1) = unzip $ S.toPair . coordinates . toSource <$> flts1
        sds1 = sqrt <$> vrs1

    let xpth0 = S.head <$> xpth
        zpth0 = S.head <$> zpth

    let mse ms = sqrt . sum $ square <$> zipWith subtract ms xpth0
        shower x = showFFloat (Just 2) x ""

    putStrLn
        $ concat
            [ "Optimal: LL: "
            , shower $ llsf kflt
            , ", Initial: LL: "
            , shower $ llsf kflt0
            , ", Learned: LL: "
            , shower $ llsf kflt1
            ]

    putStrLn "Sample Latent Path Statistics:"
    putStrLn $ concat ["Optimal MSE/Average SD: ", shower $ mse mus, "/", shower $ average sds]
    putStrLn $ concat ["Initial MSE/Average SD: ", shower $ mse mus0, "/", shower $ average sds0]
    putStrLn $ concat ["Average MSE/Average SD: ", shower $ mse mus1, "/", shower $ average sds1]

    let json =
            toJSON
                [ "conjugation-curve" .= zip3 xsmps ys yhts
                , "kalman-filter" .= DL.transpose [xpth0, zpth0, mus, sds, mus1, sds1]
                ]

    --- Process data
    rsltsfl <- resultsFilePath "factor-analysis.json"
    exportJSON rsltsfl json

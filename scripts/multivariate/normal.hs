#! stack runghc

{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

type TrueNormal = IsotropicNormal 2
type FitNormal = IsotropicNormal 2

nsmps :: Int
nsmps = 1000

mux,muy,vrx,cvr,vry :: Double
mux = 1
muy = -1
vrx = 3
cvr = 1
vry = 2

tru :: Source # TrueNormal
tru = fromTuple (mux,muy,vrx)

mn,mx :: Double
mn = -6
mx = 6

nrng :: Int
nrng = 100

--- Main ---


main :: IO ()
main = do

    smps <- realize $ sample nsmps tru

    let mtru :: Mean # TrueNormal
        mtru = toMean tru

    let mtru' :: Mean # TrueNormal
        mtru' = averageSufficientStatistic smps

    let mfit :: Mean # FitNormal
        mfit = averageSufficientStatistic smps

    let nfit :: Natural # FitNormal
        nfit = toNatural mfit

    --let nnrm :: Natural # IsotropicNormal 2
    --    nnrm = mle smps

    synthsmps <- realize $ sample 1000000 nfit

    let mrefit :: Mean # FitNormal
        mrefit = averageSufficientStatistic synthsmps

    let nrefit :: Natural # FitNormal
        nrefit = toNatural mrefit


    let ds1 = densities tru smps
    let ds2 = densities (toNatural tru) smps

    putStrLn "Densities"
    print . average $ square <$> zipWith (-) ds1 ds2
    putStrLn "True"
    print tru
    putStrLn "True Means"
    print mtru
    putStrLn "True Sample Means"
    print mtru'
    putStrLn "Mean Fit"
    print mfit
    putStrLn "Natural Fit"
    print nfit
    putStrLn "Isotransformed Fit"
    print . toNatural $ toMean nfit
    putStrLn "Resampled Fit"
    print $ nrefit
    putStrLn "Resampled Means"
    print $ mrefit

    let dsmps nrm = do
            x <- range mn mx 100
            y <- range mn mx 100
            return (x,y,density nrm $ S.doubleton x y)

        trups = dsmps tru
        lrnps = dsmps (toNatural tru)

    let ldpth = "normal"
        smpnm = "samples"
        trunm = "true-lines"
        lrnnm = "learned-lines"

    goalExport ldpth smpnm $ S.toList <$> smps

    goalExport ldpth trunm trups

    goalExport ldpth lrnnm lrnps

    runGnuplotWithVariables  ldpth "multivariate"
        [("xmn",show mn),("xmx",show mx),("ymn",show mn),("ymx",show mx)]


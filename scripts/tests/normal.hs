#! stack runghc

{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

type TrueNormal = FullNormal 2
type FitNormal = FullNormal 2

nsmps :: Int
nsmps = 20

mux,muy,vrx,cvr,vry :: Double
mux = 1
muy = -1
vrx = 3
cvr = 1
vry = 2

tru :: Source # TrueNormal
tru = fromTuple (mux,muy,vrx,cvr,vry)

mn,mx :: Double
mn = -6
mx = 6

nrng :: Int
nrng = 100

--- Main ---


main :: IO ()
main = do

    print tru
    let (mu,sgma) = split tru
        invsgma = inverse sgma
        nnrm :: Natural # TrueNormal
        nnrm = join (breakPoint $ invsgma >.> mu) . breakPoint $ (-0.5) * invsgma

    print . toSource $ toNatural tru

    smps <- realize $ sample nsmps tru

    --print smps

    let mfit :: Mean # FitNormal
        mfit = averageSufficientStatistic smps

    let nfit :: Natural # FitNormal
        nfit = toNatural mfit

    let dsmps nrm = do
            x <- range mn mx 5
            y <- range mn mx 5
            return (density nrm $ S.doubleton x y)

    print . sqrt . sum $ square <$> zipWith (-) (dsmps tru) (dsmps $ toNatural tru)
    --print . snd . split $ tru
    --print . snd . split . toSource $ toNatural tru
    --print . potential $ toNatural tru
    --print . dualPotential $ transition tru

    --let (nmu,nsgma) = split $ toNatural tru
    --    (insgma,lndt,_) = inverseLogDeterminant . negate $ 2 * nsgma
    --    ptn = 0.5 * (nmu <.> (insgma >.> nmu)) -0.5 * lndt

    --print ptn
    --let tru' :: Natural # MVNCovariance (MVNMean 2) (MVNMean 2)
    --    tru' =  fromTensor . toTensor . snd . split $ toNatural tru
    --print $ toNatural tru'
    --print $ unnormalizedLogDensities (toNatural tru) smps
    --let ldpth = "normal"
    --    smpnm = "samples"
    --    trunm = "true-lines"
    --    lrnnm = "learned-lines"

    --goalExport ldpth smpnm $ S.toList <$> smps

    --goalExport ldpth trunm trups

    --goalExport ldpth lrnnm lrnps

    --runGnuplotWithVariables  ldpth "multivariate"
    --    [("xmn",show mn),("xmx",show mx),("ymn",show mn),("ymx",show mx)]

    --let mtru :: Mean # TrueNormal
    --    mtru = toMean tru

    --let mtru' :: Mean # TrueNormal
    --    mtru' = averageSufficientStatistic smps

    ----let nnrm :: Natural # IsotropicNormal 2
    ----    nnrm = mle smps

    --synthsmps <- realize $ sample 1000000 nfit

    --let mrefit :: Mean # FitNormal
    --    mrefit = averageSufficientStatistic synthsmps

    --let nrefit :: Natural # FitNormal
    --    nrefit = toNatural mrefit


    --let ds1 = densities tru smps
    --let ds2 = densities (toNatural tru) smps

    --putStrLn "Densities"
    --print . average $ square <$> zipWith (-) ds1 ds2
    --putStrLn "True"
    --print tru
    --putStrLn "True Means"
    --print mtru
    --putStrLn "True Sample Means"
    --print mtru'
    --putStrLn "Mean Fit"
    --print mfit
    --putStrLn "Natural Fit"
    --print nfit
    --putStrLn "Isotransformed Fit"
    --print . toNatural $ toMean nfit
    --putStrLn "Resampled Fit"
    --print $ nrefit
    --putStrLn "Resampled Means"
    --print $ mrefit



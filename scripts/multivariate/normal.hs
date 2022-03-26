#! stack runghc

{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

nsmps :: Int
nsmps = 10

mux,muy,vrx,cvr,vry :: Double
mux = 1
muy = -1
vrx = 2
cvr = 1
vry = 3

tru :: Source # SymmetricNormal 2
tru = fromTuple (mux,muy,vrx,cvr,vry)

mn,mx :: Double
mn = -6
mx = 6

nrng :: Int
nrng = 100

type FitNormal = IsotropicNormal 2

--- Main ---


main :: IO ()
main = do

    smps <- realize $ sample nsmps tru

    --let nnrm :: Natural # SymmetricNormal 2
    --    nnrm = mle smps

    let nnrm :: Natural # FitNormal
        nnrm = mle smps

    --let nnrm :: Natural # IsotropicNormal 2
    --    nnrm = mle smps

    synthsmps <- realize $ sample 1000000 nnrm

    let mnrm :: Mean # FitNormal
        mnrm = averageSufficientStatistic synthsmps


    putStrLn "Fit"
    print nnrm
    putStrLn "Isotransformed Fit"
    print . toNatural $ toMean nnrm
    putStrLn "Resampled Fit"
    print $ toNatural mnrm
    putStrLn "Means"
    print $ toMean nnrm
    putStrLn "Sampled Means"
    print mnrm

    let dsmps nrm = do
            x <- range mn mx 100
            y <- range mn mx 100
            return (x,y,density nrm $ S.doubleton x y)

        trups = dsmps tru
        lrnps = dsmps nnrm

    let ldpth = "normal"
        smpnm = "samples"
        trunm = "true-lines"
        lrnnm = "learned-lines"

    goalExport ldpth smpnm $ S.toList <$> smps

    goalExport ldpth trunm trups

    goalExport ldpth lrnnm lrnps

    runGnuplotWithVariables  ldpth "multivariate"
        [("xmn",show mn),("xmx",show mx),("ymn",show mn),("ymx",show mx)]


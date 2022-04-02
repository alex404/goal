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

    smps <- realize $ sample nsmps tru

    let mfit :: Mean # FitNormal
        mfit = averageSufficientStatistic smps

    let nfit :: Natural # FitNormal
        nfit = toNatural mfit

    let dsmps nrm = do
            x <- range mn mx 100
            y <- range mn mx 100
            return (x,y,density nrm $ S.doubleton x y)

        trups = dsmps tru
        lrnps = dsmps nfit

    let ldpth = "normal"
        smpnm = "samples"
        trunm = "true-lines"
        lrnnm = "learned-lines"

    goalExport ldpth smpnm $ S.toList <$> smps

    goalExport ldpth trunm trups

    goalExport ldpth lrnnm lrnps

    runGnuplotWithVariables  ldpth "multivariate"
        [("xmn",show mn),("xmx",show mx),("ymn",show mn),("ymx",show mx)]

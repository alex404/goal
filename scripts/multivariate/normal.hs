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

tru :: Source # MultivariateNormal 2
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

    --let nnrm :: Natural # MultivariateNormal 2
    --    nnrm = mle smps

    let nnrm :: Natural # IsotropicNormal 2
        nnrm = mle smps

    print nnrm
    print $ toMean nnrm
    print . toNatural $ toMean nnrm

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


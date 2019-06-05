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
vrx = 1
cvr = 0.7
vry = 2

tru :: Source # MultivariateNormal 2
tru = Point $ S.fromTuple (mux,muy,vrx,cvr,vry)

mn,mx :: Double
mn = -5
mx = 5

nrng :: Int
nrng = 100

expmnt :: Experiment
expmnt = Experiment "probability" "multivariate"

sbexpmnt :: Maybe Analysis
sbexpmnt = Just $ Analysis "isolines" "normal"

--- Main ---


main :: IO ()
main = do

    smps <- realize $ sample nsmps tru

    let mnrm :: Mean # MultivariateNormal 2
        mnrm = sufficientStatisticT smps
        nnrm = toNatural mnrm

    let dsmps nrm = do
            x <- range mn mx 100
            y <- range mn mx 100
            return (x,y,density nrm $ S.doubleton x y)

        trups = dsmps tru
        lrnps = dsmps nnrm

    goalExport True expmnt sbexpmnt $ S.toList <$> smps

    goalExport False expmnt sbexpmnt trups

    goalExport False expmnt sbexpmnt lrnps

    runGnuplot expmnt sbexpmnt defaultGnuplotOptions "multivariate.gpi"

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

-- Misc

import Control.Monad (replicateM)

--- Globals ---

--- Distribution

type Neurons = 4
type Boltzmann' = Boltzmann Neurons

bltztru :: Natural # Boltzmann'
bltztru =
    -fromTuple
        ( 0.1
        , 0.5
        , -0.25
        , -0.5
        , 0.2
        , 0.25
        , -0.5
        , 1
        , 0.5
        , -0.1
        )

bltz0 :: Natural # Boltzmann'
bltz0 =
    fromTuple
        ( 0
        , 0
        , 0
        , 0
        , 0
        , 0
        , 0
        , 0
        , 0
        , 0
        )

--- Simulation

smpsz :: Int
smpsz = 10000

ncycs :: Int
ncycs = 1000

eps :: Double
eps = 0.05

nepchs :: Int
nepchs = 500

--- Plotting

--- Functions

fitBoltzmann :: Sample (Boltzmann Neurons) -> [Natural # Boltzmann Neurons]
fitBoltzmann blss =
    vanillaGradientSequence (logLikelihoodDifferential blss) eps defaultAdamPursuit bltz0

--- Main ---

main :: IO ()
main = do
    smp <- replicateM smpsz . realize $ gibbsBoltzmann ncycs bltztru

    let bltzs = take nepchs $ fitBoltzmann smp
        dlls = logLikelihood smp <$> bltzs

    print dlls

    mapM_ print . zip (listCoordinates bltztru) . listCoordinates $ last bltzs

--     drch = last drchs
--     dxrng = range dmn dmx pltsmps
--     dyrng = range dmn dmx pltsmps
--
-- let ddnss = do
--         y <- dyrng
--         x <- dxrng
--         return (density2d dtru (x, y), density2d drch (x, y))
--
-- let (dtdns, dldns) = unzip ddnss
--
-- let djson =
--         toJSON
--             [ "log-likelihood" .= dlls
--             , "xrange" .= dxrng
--             , "yrange" .= dyrng
--             , "true-density" .= dtdns
--             , "learned-density" .= dldns
--             , "samples" .= map S.toList dsmps
--             ]
--
--     --- Process data
--     mvnfl <- resultsFilePath "multivariate-normal.json"
--     drchfl <- resultsFilePath "multivariate-dirichlet.json"
--
--     exportJSON mvnfl njson
--     exportJSON drchfl djson

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

--- Misc

import Data.Aeson (ToJSON)

--- Globals ---

nsmps :: Int
nsmps = 20

-- Binomial --

ttlB :: String
ttlB = "binomial"

truB :: Source # Binomial 5
truB = singleton 0.3

rngB :: [Int]
rngB = pointSampleSpace truB

-- Categorical --

ttlC :: String
ttlC = "categorical"

truC :: Source # Categorical 4
truC = fromTuple (0.1, 0.4, 0.1, 0.2)

rngC :: [Int]
rngC = pointSampleSpace truC

-- Poisson --

ttlP :: String
ttlP = "poisson"

truP :: Source # Poisson
truP = singleton 5

rngP :: [Int]
rngP = [0 .. 20]

-- Normal --

ttlN :: String
ttlN = "normal"

truN :: Source # Normal
truN = fromTuple (2, 0.7)

rngN :: [Double]
rngN = range (-3) 7 100

-- Layout --

simulateDistribution ::
    forall m.
    ( Transition Mean Source m
    , DuallyFlatExponentialFamily m
    , AbsolutelyContinuous Source m
    , Generative Source m
    , AbsolutelyContinuous Natural m
    , ToJSON (SamplePoint m)
    ) =>
    String ->
    Sample m ->
    Source # m ->
    IO ()
simulateDistribution ttl rng tru = do
    smps <- realize $ sample nsmps tru

    let mmle :: Mean # m
        mmle = averageSufficientStatistic smps
        smle = toSource mmle
        nmle = toNatural mmle

    let tdns = density tru <$> rng
        sdns = density smle <$> rng
        ndns = density nmle <$> rng

    -- Create a data structure for the combined JSON output
    let jsonData =
            toJSON
                [ "title" .= ttl
                , "samples" .= smps
                , "range" .= rng
                , "true-density" .= tdns
                , "source-density" .= sdns
                , "natural-density" .= ndns
                ]

    -- Export data as JSON
    rsltfl <- resultsFilePath $ concat ["univariate-", ttl, ".json"]
    exportJSON rsltfl jsonData

main :: IO ()
main = do
    simulateDistribution ttlB rngB truB
    simulateDistribution ttlC rngC truC
    simulateDistribution ttlP rngP truP
    simulateDistribution ttlN rngN truN

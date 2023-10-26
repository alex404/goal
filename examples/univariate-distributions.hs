{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

--- Globals ---

nsmps :: Int
nsmps = 20

-- Binomial --

ttlB :: String
ttlB = "binomial"

mnB, mxB :: Double
(mnB, mxB) = (0, 5)

bnsB :: Int
bnsB = 11

truB :: Source # Binomial 5
truB = singleton 0.3

rngB :: [Int]
rngB = pointSampleSpace truB

-- Categorical --

ttlC :: String
ttlC = "categorical"

mnC, mxC :: Double
(mnC, mxC) = (0, 4)

bnsC :: Int
bnsC = 5

truC :: Source # Categorical 4
truC = fromTuple (0.1, 0.4, 0.1, 0.2)

rngC :: [Int]
rngC = pointSampleSpace truC

-- Poisson --

ttlP :: String
ttlP = "poisson"

mnP, mxP :: Double
(mnP, mxP) = (0, 20)

bnsP :: Int
bnsP = 20

truP :: Source # Poisson
truP = singleton 5

rngP :: [Int]
rngP = [0 .. 20]

-- Normal --

ttlN :: String
ttlN = "normal"

mnN, mxN :: Double
(mnN, mxN) = (-3, 7)

bnsN :: Int
bnsN = 20

truN :: Source # Normal
truN = fromTuple (2, 0.7)

rngN :: [Double]
rngN = range (-3) 7 100

-- Layout --

generateLayout ::
    forall m.
    ( Transition Mean Source m
    , LegendreExponentialFamily m
    , AbsolutelyContinuous Source m
    , Generative Source m
    , AbsolutelyContinuous Natural m
    ) =>
    String ->
    Int ->
    Double ->
    Double ->
    Sample m ->
    Source # m ->
    IO ()
generateLayout ttl nb mn mx rng tru = do
    smps <- realize $ sample nsmps tru

    let mmle :: Mean # m
        mmle = averageSufficientStatistic smps
        smle = toSource mmle
        nmle = toNatural mmle

    let tdns = density tru <$> rng
        sdns = density smle <$> rng
        ndns = density nmle <$> rng

    let msmps, mrng :: [Mean # m]
        msmps = sufficientStatistic <$> smps
        dsmps = head . listCoordinates <$> msmps
        mrng = sufficientStatistic <$> rng
        drng = head . listCoordinates <$> mrng

    let (bns, _, [wghts]) = histograms nb (Just (mn, mx)) [dsmps]

    -- Create a data structure for the combined JSON output
    let jsonData =
            toJSON
                [ "title" .= ttl
                , "range" .= drng
                , "true-density" .= tdns
                , "source-density" .= sdns
                , "natural-density" .= ndns
                , "histogram-bins" .= bns
                , "histogram-weights" .= wghts
                ]
    -- Export data as JSON

    rsltfl <- resultsFilePath $ concat ["univariate-", ttl, ".json"]
    exportJSON rsltfl jsonData

main :: IO ()
main = do
    generateLayout ttlB bnsB mnB mxB rngB truB
    generateLayout ttlC bnsC mnC mxC rngC truC
    generateLayout ttlP bnsP mnP mxP rngP truP
    generateLayout ttlN bnsN mnN mxN rngN truN

    rsltsdr <- resultsFilePath ""
    runPythonScriptWithArg "univariate-distributions.py" rsltdr

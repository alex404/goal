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

smpns :: [Int]
smpns = [20, 200, 2000]

eps :: Double
eps = 0.05

nepchs :: Int
nepchs = 500

--- Von Mises

ttlV :: String
ttlV = "von-mises"

truV :: Source # VonMises
truV = fromTuple (2, 2)

mnV, mxV :: Double
(mnV, mxV) = (0, 2 * pi)

rngV :: [Double]
rngV = range mnV mxV 200

fit0V :: Source # VonMises
fit0V = fromTuple (0, 1)

--- CoMPoisson

ttlC :: String
ttlC = "com-poisson"

truC :: Source # CoMPoisson
truC = fromTuple (5, 0.5)

mnC, mxC :: Int
(mnC, mxC) = (0, 20)

rngC :: [Int]
rngC = [mnC .. mxC]

fit0C :: Source # CoMPoisson
fit0C = fromTuple (1, 1)

--- Gamma

ttlG :: String
ttlG = "gamma"

truG :: Source # Gamma
truG = fromTuple (2, 1 / 2)

mnG, mxG :: Double
(mnG, mxG) = (0, 20)

rngG :: [Double]
rngG = tail $ range mnG mxG 200

fit0G :: Source # Gamma
fit0G = fromTuple (1, 1)

--- Main ---

simulateDistribution ::
    forall m.
    ( LegendreExponentialFamily m
    , Generative Natural m
    , AbsolutelyContinuous Natural m
    , Transition Source Natural m
    , Transition Natural Source m
    , ToJSON (SamplePoint m)
    ) =>
    String ->
    Sample m ->
    Source # m ->
    Source # m ->
    IO ()
simulateDistribution ttl rng stru fit0 = do
    --- Training

    let ntru = toNatural stru

    smps <- mapM realize $ (`sample` ntru) <$> smpns
    let smp = last smps

    let fits =
            take nepchs
                . vanillaGradientSequence (exponentialFamilyLogLikelihoodDifferential smp) eps defaultAdamPursuit
                $ toNatural fit0

    let nfit :: Natural # m
        nfit = last fits

    let tdns = density ntru <$> rng
        ndns = density nfit <$> rng

    --- Sufficient Statistics

    let mtru = toMean ntru
        mavgs :: [Mean # m]
        mavgs = averageSufficientStatistic <$> smps

    --- JSON

    let jsonData =
            toJSON
                [ "title" .= ttl
                , "samples" .= smps
                , "true-statistics" .= listCoordinates mtru
                , "sample-statistics" .= (listCoordinates <$> mavgs)
                , "range" .= rng
                , "true-density" .= tdns
                , "natural-density" .= ndns
                ]

    rsltfl <- resultsFilePath $ concat ["univariate-", ttl, ".json"]
    exportJSON rsltfl jsonData

main :: IO ()
main = do
    simulateDistribution ttlV rngV truV fit0V
    simulateDistribution ttlC rngC truC fit0C
    simulateDistribution ttlG rngG truG fit0G

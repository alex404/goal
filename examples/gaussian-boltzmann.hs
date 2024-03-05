{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Globals ---

--- Model

type Neurons = 9
type ObservationSpace = 2

unibnd :: Double
unibnd = 4

--- Plotting

pltsmps :: Int
pltsmps = 100

--- Training

eps :: Double
eps = 0.03

nstps, nits :: Int
nstps = 200
nits = 20

gp :: GradientPursuit
gp = defaultAdamPursuit

emStep ::
    (KnownNat n, KnownNat k, 1 <= k) =>
    [S.Vector n Double] ->
    (Int, Natural # GaussianBoltzmannHarmonium L.PositiveDefinite n k) ->
    IO (Int, Natural # GaussianBoltzmannHarmonium L.PositiveDefinite n k)
emStep xs (k, gbhrm) = do
    let gbhrms = expectationMaximizationAscent eps gp xs gbhrm
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show $ logLikelihood xs gbhrm
            ]
    return (k + 1, gbhrms !! nstps)

--- Main ---

main :: IO ()
main = do
    --- Random Ground-Truth Gaussian-Boltzmann machine
    bltztru :: Natural # Boltzmann Neurons <-
        realize $ uniformInitialize (-unibnd, unibnd)

    shftstru :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli) <-
        realize $ uniformInitialize (-unibnd, unibnd)
    let mvntru = standardNormal
        lmdltru = join mvntru shftstru
        gbhrmtru :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace Neurons
        gbhrmtru = joinConjugatedHarmonium lmdltru bltztru

    --- Random Initial Gaussian-Boltzmann machine
    bltz0 :: Natural # Boltzmann Neurons <-
        realize $ uniformInitialize (-unibnd, unibnd)
    shfts0 :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli) <-
        realize $ uniformInitialize (-unibnd, unibnd)

    let mvn0 = standardNormal
        lmdl0 = join mvn0 shfts0
        gbhrm0 = joinConjugatedHarmonium lmdl0 bltz0

    --- Sample from the Gaussian-Boltzmann machine
    xzsmps :: [(S.Vector ObservationSpace Double, S.Vector Neurons Bool)] <- realize $ sample 1000 gbhrmtru
    let xsmps = fst <$> xzsmps

    --- Calculate plot centres
    let smps = fst <$> xzsmps
        xs = flip S.index 0 <$> smps
        ys = flip S.index 1 <$> smps
        mux = average xs
        muy = average ys
        mnx = mux - 10
        mxx = mux + 10
        mny = muy - 10
        mxy = muy + 10

    let xrng = range mnx mxx pltsmps
        yrng = range mny mxy pltsmps
        xyplts0 = [(x, y) | y <- yrng, x <- xrng]
        xyplts = S.fromTuple <$> xyplts0

    kgbhrms <- iterateM nits (emStep xsmps) (0, gbhrm0)

    let gbhrms = snd <$> kgbhrms

    --- Density

    let trudns = observableDensities gbhrmtru xyplts
        dns0 = observableDensities gbhrm0 xyplts
        lrndns = observableDensities (last gbhrms) xyplts

    let json =
            toJSON
                [ "xys" .= xyplts0
                , "true-density" .= trudns
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                ]
    --- Process data
    flnm <- resultsFilePath "gaussian-boltzmann.json"

    exportJSON flnm json

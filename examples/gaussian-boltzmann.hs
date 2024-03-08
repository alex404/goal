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

--- Initialisation

unibnd :: Double
unibnd = 4

initializeGaussianBoltzmannHarmonium ::
    Random (Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace Neurons)
initializeGaussianBoltzmannHarmonium = do
    bltz :: Natural # Boltzmann Neurons <- uniformInitialize (-unibnd, unibnd)
    shfts :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli) <-
        uniformInitialize (-unibnd, unibnd)
    let mvn = standardNormal
        lmdl = join mvn shfts
    return $ joinConjugatedHarmonium lmdl bltz

bltztru :: Natural # Boltzmann Neurons
bltztru = join (fromTuple (1, 1, 1, 1, 0, 1, 1, 1, 1)) $ -10

mubnd :: Double
mubnd = 4

shftstru :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli)
shftstru =
    fromRows $
        S.fromTuple
            ( fromTuple (-mubnd, -mubnd, -mubnd, 0, 0, 0, mubnd, mubnd, mubnd)
            , fromTuple (-mubnd, 0, mubnd, -mubnd, 0, mubnd, -mubnd, 0, mubnd)
            )

smvntru :: Source # FullNormal ObservationSpace
smvntru = standardNormal / 2

mvntru :: Natural # FullNormal ObservationSpace
mvntru = toNatural smvntru

lmdltru :: Natural # BoltzmannLinearModel L.PositiveDefinite ObservationSpace Neurons
lmdltru = join mvntru shftstru

gbhrmtru :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace Neurons
gbhrmtru = joinConjugatedHarmonium lmdltru bltztru

--- Plotting

pltsmps :: Int
pltsmps = 100

--- Training

eps :: Double
eps = 0.03

nstps, nepchs :: Int
nstps = 200
nepchs = 10

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
    -- gbhrm0 <- realize initializeGaussianBoltzmannHarmonium
    -- let (bss, intrs) = split . snd $ splitConjugatedHarmonium gbhrm0
    --     bsscrds = listCoordinates bss
    --     intrscrds = listCoordinates intrs
    -- print ("Random Biases:" :: String)
    -- print $ roundSD 3 <$> bsscrds
    -- print ("Random Interactions:" :: String)
    -- print $ roundSD 3 <$> intrscrds

    --- Sample from the Gaussian-Boltzmann machine
    -- xzsmps :: [(S.Vector ObservationSpace Double, S.Vector Neurons Bool)] <- realize $ sample 1000 gbhrmtru
    -- let xsmps = fst <$> xzsmps

    --- Calculate plot centres
    -- let smps = fst <$> xzsmps
    --     xs = flip S.index 0 <$> smps
    --     ys = flip S.index 1 <$> smps
    --     mux = average xs
    --     muy = average ys
    --     mnx = mux - 10
    --     mxx = mux + 10
    --     mny = muy - 10
    --     mxy = muy + 10

    let mnx = -10
        mxx = 10
        mny = -10
        mxy = 10

    let xrng = range mnx mxx pltsmps
        yrng = range mny mxy pltsmps
        xyplts0 = [(x, y) | y <- yrng, x <- xrng]
        xyplts = S.fromTuple <$> xyplts0

    -- kgbhrms <- iterateM nepchs (emStep xsmps) (0, gbhrm0)

    -- let gbhrms = snd <$> kgbhrms

    --- Density

    let trudns = observableDensities gbhrmtru xyplts
    -- dns0 = observableDensities gbhrm0 xyplts
    -- lrndns = observableDensities (last gbhrms) xyplts

    let json =
            toJSON
                [ "xys" .= xyplts0
                , "true-density" .= trudns
                , "initial-density" .= trudns
                , "learned-density" .= trudns
                -- , "initial-density" .= dns0
                -- , "learned-density" .= lrndns
                ]
    --- Process data
    flnm <- resultsFilePath "gaussian-boltzmann.json"

    exportJSON flnm json

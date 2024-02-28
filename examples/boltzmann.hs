{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Globals ---

--- Prior

type Neurons = 9
type ObservationSpace = 2

--- Plotting

pltsmps :: Int
pltsmps = 200

--- Main ---

main :: IO ()
main = do
    --- Random Boltzmann machine
    bltz :: Natural # Boltzmann Neurons <- realize $ uniformInitialize (-5, 5)

    let shfts :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli)
        shfts =
            fromRows $
                S.fromTuple
                    ( fromTuple (2, 2, 2, 0, 0, 0, -2, -2, -2)
                    , fromTuple (2, 0, -2, 2, 0, -2, 2, 0, -2)
                    )

    --- Linear Model
    let mvn :: Natural # FullNormal ObservationSpace
        mvn = standardNormal

        lmdl :: Natural # BoltzmannLinearModel L.PositiveDefinite ObservationSpace Neurons
        lmdl = join mvn shfts

    --- Harmonium
    let hrm :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace Neurons
        hrm = joinConjugatedHarmonium lmdl bltz

        rprms :: (Double, Natural # Boltzmann Neurons)
        rprms = linearBoltzmannHarmoniumConjugationParameters lmdl

    --- Calculate statistics
    xzsmps :: [(S.Vector ObservationSpace Double, S.Vector Neurons Bool)] <- realize $ sample 100 hrm
    let smps = fst <$> xzsmps
        xs = flip S.index 0 <$> smps
        ys = flip S.index 1 <$> smps
        (mux, sdx) = estimateMeanSD xs
        (muy, sdy) = estimateMeanSD ys
        mnx = mux - 5
        mxx = mux + 5
        mny = muy - 5
        mxy = muy + 5

        xrng = range mnx mxx pltsmps
        yrng = range mny mxy pltsmps

    --- Lines
    let dnss = exp <$> logConjugatedDensities rprms hrm [S.fromTuple (x, y) | y <- yrng, x <- xrng]

    let json =
            toJSON
                [ "xrange" .= xrng
                , "yrange" .= yrng
                , "density" .= dnss
                ]
    --- Process data
    flnm <- resultsFilePath "boltzmann.json"

    exportJSON flnm json

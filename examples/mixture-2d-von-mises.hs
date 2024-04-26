{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeOperators #-}

--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Data.List qualified as L
import Goal.Core.Vector.Storable qualified as S

--- Globals ---

-- Distributions --

vm1, vm2, vm3 :: Source # (VonMises, VonMises)
vm1 = fromTuple (2.5, 3, 4.5, 2)
vm2 = fromTuple (1.5, 2, 1.5, 3)
vm3 = fromTuple (4.5, 2.5, 1.5, 2.5)

wghts :: Source # Categorical 2
wghts = fromTuple (0.33, 0.33)

sxz :: Source # Mixture (VonMises, VonMises) 2
sxz = joinSourceMixture (S.fromTuple (vm1, vm2, vm3)) wghts

trunxz :: Natural # Mixture (VonMises, VonMises) 2
trunxz = transition sxz

-- Training --

nsmps :: Int
nsmps = 100

eps :: Double
eps = 5e-3

nepchs :: Int
nepchs = 100

nbtch, nmlt :: Int
nbtch = 10
nmlt = div nsmps nbtch

nitrs :: Int
nitrs = 200

emGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    [Natural # Mixture (VonMises, VonMises) 2]
emGP zs nxz0 =
    take nitrs
        . flip iterate nxz0
        $ \nxz -> expectationMaximizationAscent eps defaultAdamPursuit zs nxz !! nepchs

mlGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    [Natural # Mixture (VonMises, VonMises) 2]
mlGP zs nxz0 =
    take nitrs
        . takeEvery nepchs
        $ vanillaGradientSequence (logLikelihoodDifferential zs) eps defaultAdamPursuit nxz0

emMCGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    Random [Natural # Mixture (VonMises, VonMises) 2]
emMCGP xs =
    iterateM
        nitrs
        (iterateChain (nmlt * nepchs) . stochasticMonteCarloExpectationMaximizationAscent eps defaultAdamPursuit xs nbtch nbtch)

mlMCGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    Random [Natural # Mixture (VonMises, VonMises) 2]
mlMCGP xs nxz0 = do
    nxzs <-
        streamChain (nmlt * nepchs * nitrs) $
            stochasticMonteCarloMaximumLikelihoodAscent eps defaultAdamPursuit xs nbtch nbtch nxz0
    return $ takeEvery (nmlt * nepchs) nxzs

-- Plotting --

ellipse :: Source # (VonMises, VonMises) -> [(Double, Double)]
ellipse vm =
    let [mux, kpx, muy, kpy] = listCoordinates vm
        f t = (mux + sqrt (recip kpx) * cos t, muy + sqrt (recip kpy) * sin t)
     in f <$> range 0 (2 * pi) 100

mixtureModelToConfidenceCSV :: Natural # Mixture (VonMises, VonMises) 2 -> [[(Double, Double)]]
mixtureModelToConfidenceCSV nxz =
    let cmps = S.toList . fst $ splitNaturalMixture nxz
     in ellipse . toSource <$> cmps

--- Main ---

main :: IO ()
main = do
    xzs <- realize $ sample nsmps trunxz
    nxz0 <- realize $ uniformInitialize (-0.2, 0.2)

    let xs = fst <$> xzs

    let emnxzs = emGP xs nxz0
        mlnxzs = mlGP xs nxz0

    semnxzs <- realize $ emMCGP xs nxz0
    smlnxzs <- realize $ mlMCGP xs nxz0

    let truces = repeat . negate $ logLikelihood xs trunxz
        emces = negate . logLikelihood xs <$> emnxzs
        mlces = negate . logLikelihood xs <$> mlnxzs
        semces = negate . logLikelihood xs <$> semnxzs
        smlces = negate . logLikelihood xs <$> smlnxzs

    let cecsvs = L.zip5 truces emces mlces semces smlces

    let [trucnfs, emcnfs, mlcnfs, semcnfs, smlcnfs] =
            mixtureModelToConfidenceCSV <$> [trunxz, last emnxzs, last mlnxzs, last semnxzs, last smlnxzs]

    let json =
            toJSON
                [ "cross-entropy" .= cecsvs
                , "true-confidence" .= trucnfs
                , "em-confidence" .= emcnfs
                , "ml-confidence" .= mlcnfs
                , "sem-confidence" .= semcnfs
                , "sml-confidence" .= smlcnfs
                , "samples" .= xs
                ]

    jsnpth <- resultsFilePath "mixture-2d-von-mises.json"

    exportJSON jsnpth json

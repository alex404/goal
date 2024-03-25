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
vm1 = fromTuple (2.5, 2, 4.5, 1.5)
vm2 = fromTuple (1.5, 1.5, 1.5, 2.5)
vm3 = fromTuple (4.5, 3, 1.5, 1.5)

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
eps = 0.01

nstps :: Int
nstps = 200

nbtch :: Int
nbtch = 10

nepchs :: Int
nepchs = 100

emGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    [Natural # Mixture (VonMises, VonMises) 2]
emGP zs nxz0 =
    take nepchs
        . flip iterate nxz0
        $ \nxz -> expectationMaximizationAscent eps defaultAdamPursuit zs nxz !! nstps

mlGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    [Natural # Mixture (VonMises, VonMises) 2]
mlGP zs nxz0 =
    take nepchs
        . takeEvery nstps
        $ vanillaGradientSequence (logLikelihoodDifferential zs) eps defaultAdamPursuit nxz0

emSGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    Random [Natural # Mixture (VonMises, VonMises) 2]
emSGP xs =
    iterateM
        nepchs
        (iterateChain nstps . stochasticConjugatedEMAscent eps defaultAdamPursuit xs nbtch)

mlSGP ::
    -- | Observations
    Sample (VonMises, VonMises) ->
    Natural # Mixture (VonMises, VonMises) 2 ->
    Random [Natural # Mixture (VonMises, VonMises) 2]
mlSGP xs nxz0 = do
    nxzs <- streamChain (nstps * nepchs) $ stochasticConjugatedMLAscent eps defaultAdamPursuit xs nbtch nxz0
    return $ takeEvery nstps nxzs

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

    semnxzs <- realize $ emSGP xs nxz0
    smlnxzs <- realize $ mlSGP xs nxz0

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

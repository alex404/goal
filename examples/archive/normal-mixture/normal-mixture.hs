#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S

import qualified Data.List as L


--- Globals ---

-- Manifolds --

-- Mixture Distributions --

mu1,mu2,mu3 :: Source # MVNMean 2
mu1 = fromTuple (0, 1)
mu2 = fromTuple (-1, -1)
mu3 = fromTuple (1, -1)

cvr1,cvr2,cvr3 :: Source # MVNCovariance (MVNMean 2) (MVNMean 2)
cvr1 = fromTuple (0.5,0.1,0.3)
cvr2 = fromTuple (0.2,-0.1,0.4)
cvr3 = fromTuple (0.3,0.1,0.5)

nrm1,nrm2,nrm3 :: Source # FullNormal 2
nrm1 = join mu1 cvr1
nrm2 = join mu2 cvr2
nrm3 = join mu3 cvr3

nrms :: S.Vector 3 (Natural # FullNormal 2)
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.33
mix2 = 0.33

wghts :: Source # Categorical 2
wghts = fromTuple (mix1,mix2)

truhrm :: Natural # Mixture (FullNormal 2) 2
truhrm = joinNaturalMixture nrms $ toNatural wghts

-- Mixture Distributions --

identity :: Source # Diagonal (MVNMean 2) (MVNMean 2)
identity = 1.5

nrm1',nrm2',nrm3' :: Source # FullNormal 2
nrm1' = join (fromTuple (-0.2,0)) . fromTensor $ toTensor identity
nrm2' = join (fromTuple (0,0)) . fromTensor $ toTensor identity
nrm3' = join (fromTuple (0.2,0)) . fromTensor $ toTensor identity

nrms' :: S.Vector 3 (Natural # FullNormal 2)
nrms' = S.fromTuple (toNatural nrm1',toNatural nrm2',toNatural nrm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Categorical 2
wghts' = fromTuple (mix1',mix2')

hrm0 :: Natural # Mixture (FullNormal 2) 2
hrm0 = joinNaturalMixture nrms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 100

nepchs :: Int
nepchs = 101

-- Functions --


emGP
    :: Sample (FullNormal 2) -- ^ Observations
    -> Int
    -> Natural # Mixture (FullNormal 2) 2
    -> Natural # Mixture (FullNormal 2) 2
emGP zs n nhrm =
    expectationMaximizationAscent 0.05 defaultAdamPursuit zs nhrm !! n

llGP
    :: Sample (FullNormal 2) -- ^ Observations
    -> Natural # Mixture (FullNormal 2) 2
    -> [Natural # Mixture (FullNormal 2) 2]
llGP zs nhrm = vanillaGradientSequence (logLikelihoodDifferential zs) 0.02 defaultAdamPursuit nhrm

emSGP
    :: Int
    -> Sample (FullNormal 2) -- ^ Observations
    -> Natural # Mixture (FullNormal 2) 2
    -> Random [Natural # Mixture (FullNormal 2) 2]
emSGP n xs nhrm =  iterateM n
    (iterateChain 1000 . stochasticConjugatedEMAscent 1e-4 defaultAdamPursuit xs 100) nhrm


--emCD
--    :: Sample (Normal,Normal) -- ^ Observations
--    -> Natural # Mixture (Normal,Normal) 2
--    -> Random (Natural # Mixture (Normal,Normal) 2)
--emCD zs nhrm = do
--    iterateChain 100 $ gibbsExpectationMaximization eps 2 nsmps defaultAdamPursuit zs nhrm

--filterCat :: [SamplePoint (Mixture (Normal,Normal) 2)] -> Int -> [SamplePoint (Normal,Normal)]
--filterCat cxys n = fst <$> filter ((== n) . snd) cxys

ellipse :: Source # FullNormal 2 -> [(Double,Double)]
ellipse nrm =
    let (xs,ys) = unzip $ bivariateNormalConfidenceEllipse 100 1 nrm
     in zip xs ys

--ellipse :: Source # (Normal, Normal) -> [ConfidenceEllipse]
--ellipse nrm =
--    let [mux,vrx,muy,vry] = listCoordinates nrm
--        f t = ConfidenceEllipse (mux + sqrt vrx * cos t) (muy + sqrt vry * sin t)
--     in f <$> range 0 (2*pi) 100


mixtureModelToConfidenceCSV :: Natural # Mixture (FullNormal 2) 2 -> [[(Double,Double)]]
mixtureModelToConfidenceCSV hrm =
    let cmps = S.toList . fst $ splitNaturalMixture hrm
     in ellipse . toSource <$> cmps


--mixtureModelToMeanCSV :: Natural # Mixture (FullNormal 2) 2 -> [(Double,Double)]
--mixtureModelToMeanCSV hrm =
--    let cmps = S.toList . fst $ splitNaturalMixture hrm
--     in S.toPair . coordinates . fst . split . toSource <$> cmps


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps truhrm

    let xys = fst <$> cxys

    let emhrms = take nepchs $ iterate (expectationMaximization xys) hrm0
        --llgahrms = take nepchs . takeEvery 100 $ llGP xys hrm0
        emgahrms = take nepchs $ iterate (emGP xys 100) hrm0

    --llsahrms0 <- realize . streamChain (200 * nepchs)
    --    $ stochasticConjugatedMLAscent 3e-3 defaultAdamPursuit xys 100 hrm0
    llsahrms0 <- realize $ emSGP nepchs xys hrm0
    let llgahrms = llsahrms0

    let trunlls = repeat . average $ negate . logObservableDensity truhrm <$> xys
        emnlls = [ average $ negate . logObservableDensity hrm <$> xys | hrm <- emhrms ]
        llganlls = [ average $ negate . logObservableDensity hrm <$> xys | hrm <- llgahrms ]
        emganlls = [ average $ negate . logObservableDensity hrm <$> xys | hrm <- emgahrms ]

    let cedcsvs = L.zip4 trunlls emnlls llganlls emganlls

    let emhrm = last emhrms
        llgahrm = last llgahrms
        emgahrm = last emgahrms
        [trucnfs,cnfs0,emcnfs,llgacnfs,emgacnfs] = mixtureModelToConfidenceCSV
            <$> [truhrm,hrm0,emhrm,llgahrm,emgahrm]

    --let [trumn,emmn,llgamn,emgamn] = mixtureModelToMeanCSV <$> [truhrm,emhrm,llgahrm,emgahrm]

    let ldpth = "data"

    goalExport ldpth "log-likelihood" cedcsvs

    goalExportLines ldpth "true-confidence" trucnfs
    goalExportLines ldpth "confidence0" cnfs0
    goalExportLines ldpth "em-confidence" emcnfs
    goalExportLines ldpth "llga-confidence" llgacnfs
    goalExportLines ldpth "emga-confidence" emgacnfs

    --goalExport ldpth "true-mean" trumn
    --goalExport ldpth "em-mean" emmn
    --goalExport ldpth "llga-mean" llgamn
    --goalExport ldpth "emga-mean" emgamn

    goalExport ldpth "samples" $ S.toPair <$> xys

    runGnuplot ldpth "log-likelihood-ascent"
    runGnuplot ldpth "mixture-components"

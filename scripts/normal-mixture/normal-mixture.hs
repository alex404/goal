#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

-- Manifolds --

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 1.5
muy1 = 4.5
mux2 = 3
muy2 = 3
mux3 = 4.5
muy3 = 1.5

vrx1,vrx2,vrx3,vry1,vry2,vry3 :: Double
vrx1 = 0.5
vry1 = 0.5
vrx2 = 0.5
vry2 = 0.5
vrx3 = 0.5
vry3 = 0.5

nrm1,nrm2,nrm3 :: Source # (Normal,Normal)
nrm1 = fromTuple (mux1, vrx1, muy1, vry1)
nrm2 = fromTuple (mux2, vrx2, muy2, vry2)
nrm3 = fromTuple (mux3, vrx3, muy3, vry3)

nrms :: S.Vector 3 (Natural # (Normal,Normal))
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Categorical 2
wghts = fromTuple (mix1,mix2)

truhrm :: Natural # Mixture (Normal,Normal) 2
truhrm = joinNaturalMixture nrms $ toNatural wghts

-- Mixture Distributions --

nrm1',nrm2',nrm3' :: Source # (Normal,Normal)
nrm1' = fromTuple (2, 0.5, pi, 0.5)
nrm2' = fromTuple (3, 0.5, pi, 0.5)
nrm3' = fromTuple (4, 0.5, pi, 0.5)

nrms' :: S.Vector 3 (Natural # (Normal,Normal))
nrms' = S.fromTuple (toNatural nrm1',toNatural nrm2',toNatural nrm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Categorical 2
wghts' = fromTuple (mix1',mix2')

hrm0 :: Natural # Mixture (Normal,Normal) 2
hrm0 = joinNaturalMixture nrms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 50

eps :: Double
eps = 0.01

bnd :: Double
bnd = 1e-2

admmlt :: Int
admmlt = 10

-- EM
nepchs :: Int
nepchs = 200

-- Functions --


emGD
    :: Sample (Normal,Normal) -- ^ Observations
    -> Natural # Mixture (Normal,Normal) 2
    -> Natural # Mixture (Normal,Normal) 2
emGD zs nhrm =
    expectationMaximizationAscent eps defaultAdamPursuit zs nhrm !! 100

emCD
    :: Sample (Normal,Normal) -- ^ Observations
    -> Natural # Mixture (Normal,Normal) 2
    -> Random r (Natural # Mixture (Normal,Normal) 2)
emCD zs nhrm = do
    iterateChain 100 $ gibbsExpectationMaximization eps 2 nsmps defaultAdamPursuit zs nhrm

filterCat :: [SamplePoint (Mixture (Normal,Normal) 2)] -> Int -> [SamplePoint (Normal,Normal)]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

ellipse :: Source # (Normal, Normal) -> [ConfidenceEllipse]
ellipse nrm0 =
    let [mux,vrx,muy,vry] = listCoordinates nrm0
        mu = S.fromTuple (mux,muy)
        cvr = S.diagonalMatrix (S.fromTuple (vrx,vry))
        nrm = joinMultivariateNormal mu cvr
        (xs,ys) = unzip $ bivariateNormalConfidenceEllipse 100 1 nrm
     in zipWith ConfidenceEllipse xs ys

mixtureModelToConfidenceCSV :: Natural # Mixture (Normal,Normal) 2 -> [[ConfidenceEllipse]]
mixtureModelToConfidenceCSV hrm =
    let cmps = S.toList . fst $ splitNaturalMixture hrm
     in ellipse . toSource <$> cmps


mixtureModelToMeanCSV :: Natural # Mixture (Normal,Normal) 2 -> [NormalMeans]
mixtureModelToMeanCSV hrm =
    let cmps = S.toList . fst $ splitNaturalMixture hrm
     in [ NormalMeans mu1 mu2 | [mu1,_,mu2,_] <- listCoordinates . toSource <$> cmps ]

-- CSV --

data CrossEntropyDescent = CrossEntropyDescent
    { trueCrossEntropy :: Double
    , emCrossEntropy :: Double
    , cdCrossEntropy :: Double }
    deriving (Generic, Show)

instance FromNamedRecord CrossEntropyDescent
instance ToNamedRecord CrossEntropyDescent
instance DefaultOrdered CrossEntropyDescent
instance NFData CrossEntropyDescent

-- Harmonium Number: 0 <=> True, 1 <=> Expectation-Maximization, 2 <=> Adam Descent
data ConfidenceEllipse = ConfidenceEllipse
    { precisionX :: Double
    , precisionY :: Double
    } deriving (Generic, Show)

instance FromNamedRecord ConfidenceEllipse
instance ToNamedRecord ConfidenceEllipse
instance DefaultOrdered ConfidenceEllipse
instance NFData ConfidenceEllipse

data TrainingSamples = TrainingSamples
    { xVal :: Double
    , yVal :: Double
    } deriving (Generic, Show)

instance FromNamedRecord TrainingSamples
instance ToNamedRecord TrainingSamples
instance DefaultOrdered TrainingSamples
instance NFData TrainingSamples

data NormalMeans = NormalMeans
    { xMean :: Double
    , yMean :: Double
    } deriving (Generic, Show)

instance FromNamedRecord NormalMeans
instance ToNamedRecord NormalMeans
instance DefaultOrdered NormalMeans
instance NFData NormalMeans


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps truhrm

    let xys = hHead <$> cxys

    let emhrms = take nepchs $ iterate (emGD xys) hrm0

    cdhrms <- realize $ iterateM nepchs (emCD xys) hrm0

    let trunlls = repeat . average $ negate . log . mixtureDensity truhrm <$> xys
        emnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- emhrms ]
        cdnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- cdhrms ]

    let cedcsvs = zipWith3 CrossEntropyDescent trunlls emnlls cdnlls

    let emhrm1 = last emhrms
        cdhrm1 = last cdhrms
        [trucnfs,emcnfs,cdcnfs] = mixtureModelToConfidenceCSV <$> [truhrm,emhrm1,cdhrm1]

    let [trumn,emmn,cdmn] = mixtureModelToMeanCSV <$> [truhrm,emhrm1,cdhrm1]

    let xycsv = [ TrainingSamples x y | (x,y) <- xys ]

    let ldpth = "data"

    goalExportNamed ldpth "log-likelihood" cedcsvs

    goalExportNamedLines ldpth "true-confidence" trucnfs
    goalExportNamedLines ldpth "em-confidence" emcnfs
    goalExportNamedLines ldpth "cd-confidence" cdcnfs

    goalExportNamed ldpth "true-mean" trumn
    goalExportNamed ldpth "em-mean" emmn
    goalExportNamed ldpth "cd-mean" cdmn

    goalExportNamed ldpth "samples" xycsv

    runGnuplot ldpth "log-likelihood-ascent"
    runGnuplot ldpth "mixture-components"

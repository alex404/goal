#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

-- VonMises Components --

vm1,vm2,vm3 :: Source # (VonMises,VonMises)
vm1 = fromTuple (1.5, 2, 4.5, 2)
vm2 = fromTuple (3, 2, 3, 2)
vm3 = fromTuple (4.5, 2, 1.5, 2)

vms :: S.Vector 3 (Natural # (VonMises,VonMises))
vms = S.fromTuple (toNatural vm1,toNatural vm2,toNatural vm3)

wghts :: Source # Categorical 2
wghts = fromTuple (0.25,0.25)

truhrm :: Natural # Mixture (VonMises,VonMises) 2
truhrm = joinNaturalMixture vms $ toNatural wghts

-- Mixture Distributions --

vm1',vm2',vm3' :: Source # (VonMises,VonMises)
vm1' = fromTuple (2, 0.5, pi, 0.5)
vm2' = fromTuple (3, 0.5, pi, 0.5)
vm3' = fromTuple (4, 0.5, pi, 0.5)

vms' :: S.Vector 3 (Natural # (VonMises,VonMises))
vms' = S.fromTuple (toNatural vm1',toNatural vm2',toNatural vm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Categorical 2
wghts' = fromTuple (mix1',mix2')

hrm0 :: Natural # Mixture (VonMises,VonMises) 2
hrm0 = joinNaturalMixture vms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 50

eps :: Double
eps = 0.05

bnd :: Double
bnd = 1e-5

admmlt :: Int
admmlt = 10

-- EM
nepchs :: Int
nepchs = 200

-- Functions --


vonMisesEM
    :: Sample (VonMises,VonMises) -- ^ Observations
    -> Natural # Mixture (VonMises,VonMises) 2
    -> Natural # Mixture (VonMises,VonMises) 2
vonMisesEM zs nhrm =
    cauchyLimit euclideanDistance bnd $ expectationMaximizationAscent eps defaultAdamPursuit zs nhrm

filterCat :: [SamplePoint (Mixture (VonMises,VonMises) 2)] -> Int -> [SamplePoint (VonMises,VonMises)]
filterCat cxys n = fst <$> filter ((== n) . snd) cxys

ellipse :: Source # (VonMises, VonMises) -> [ConfidenceEllipse]
ellipse vm =
    let [mux,kpx,muy,kpy] = listCoordinates vm
        f t = ConfidenceEllipse (mux + sqrt (recip kpx) * cos t) (muy + sqrt (recip kpy) * sin t)
     in f <$> range 0 (2*pi) 100

mixtureModelToConfidenceCSV :: Natural # Mixture (VonMises,VonMises) 2 -> [[ConfidenceEllipse]]
mixtureModelToConfidenceCSV hrm =
    let cmps = S.toList . fst $ splitNaturalMixture hrm
     in ellipse . toSource <$> cmps

mixtureModelToMeanCSV :: Natural # Mixture (VonMises,VonMises) 2 -> [VonMisesMeans]
mixtureModelToMeanCSV hrm =
    let cmps = S.toList . fst $ splitNaturalMixture hrm
     in [ VonMisesMeans mu1 mu2 | [mu1,_,mu2,_] <- listCoordinates . toSource <$> cmps ]

-- CSV --

data CrossEntropyDescent = CrossEntropyDescent
    { trueCrossEntropy :: Double
    , emCrossEntropy :: Double
    , sgdCrossEntropy :: Double }
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

data VonMisesMeans = VonMisesMeans
    { xMean :: Double
    , yMean :: Double
    } deriving (Generic, Show)

instance FromNamedRecord VonMisesMeans
instance ToNamedRecord VonMisesMeans
instance DefaultOrdered VonMisesMeans
instance NFData VonMisesMeans


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps truhrm

    let xys = fst <$> cxys

    let emhrms = take nepchs $ iterate (vonMisesEM xys) hrm0
        --itrhrm1 :: Natural # Mixture VonMises 2
        --(itrhrm1,itrlls) = iterativeMixtureOptimization
        --    (div nepchs 3) eps defaultAdamPursuit Nothing 0.1 xys $ toNatural vm2'

    let admhrms = take nepchs . takeEvery admmlt
            $ vanillaGradientSequence (logLikelihoodDifferential xys) eps defaultAdamPursuit hrm0

    let trunlls = repeat . average $ negate . logObservableDensity truhrm <$> xys
        emnlls = [ average $ negate . logObservableDensity hrm <$> xys | hrm <- emhrms ]
        admnlls = [ average $ negate . logObservableDensity hrm <$> xys | hrm <- admhrms ]
        --itrnlls = negate <$> concat itrlls

    let cedcsvs = zipWith3 CrossEntropyDescent trunlls emnlls admnlls

    let emhrm1 = last emhrms
        admhrm1 = last admhrms
        [trucnfs,emcnfs,admcnfs] = mixtureModelToConfidenceCSV <$> [truhrm,emhrm1,admhrm1]

    let [trumn,emmn,admmn] = mixtureModelToMeanCSV <$> [truhrm,emhrm1,admhrm1]

    let xycsv = [ TrainingSamples x y | (x,y) <- xys ]

    let ldpth = "data"

    goalExportNamed ldpth "log-likelihood" cedcsvs

    goalExportNamedLines ldpth "true-confidence" trucnfs
    goalExportNamedLines ldpth "em-confidence" emcnfs
    goalExportNamedLines ldpth "adam-confidence" admcnfs

    goalExportNamed ldpth "true-mean" trumn
    goalExportNamed ldpth "em-mean" emmn
    goalExportNamed ldpth "adam-mean" admmn

    goalExportNamed ldpth "samples" xycsv

    runGnuplot ldpth "cross-entropy-descent"
    runGnuplot ldpth "mixture-components"

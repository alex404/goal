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

wghts :: Source # Categorical 2
wghts = fromTuple (0.25,0.25)

strumxmdl :: Source # Mixture (VonMises,VonMises) 2
strumxmdl = joinSourceMixture (S.fromTuple (vm1,vm2,vm3)) wghts

trumxmdl :: Natural # Mixture (VonMises,VonMises) 2
trumxmdl = transition strumxmdl

-- Initialization Distribution
mxmdlint :: Source # Normal
mxmdlint = fromTuple (0,0.1)

-- Mixture Distributions --

--vm1',vm2',vm3' :: Source # (VonMises,VonMises)
--vm1' = fromTuple (2, 0.5, pi, 0.5)
--vm2' = fromTuple (3, 0.5, pi, 0.5)
--vm3' = fromTuple (4, 0.5, pi, 0.5)
--
--vms' :: S.Vector 3 (Natural # (VonMises,VonMises))
--vms' = S.fromTuple (toNatural vm1',toNatural vm2',toNatural vm3')
--
--mix1',mix2' :: Double
--mix1' = 0.33
--mix2' = 0.33
--
--wghts' :: Source # Categorical 2
--wghts' = fromTuple (mix1',mix2')
--
--hrm0 :: Natural # Mixture (VonMises,VonMises) 2
--hrm0 = joinNaturalMixture vms' $ toNatural wghts'

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
vonMisesEM zs nmxmdl =
    cauchyLimit euclideanDistance bnd $ expectationMaximizationAscent eps defaultAdamPursuit zs nmxmdl

ellipse :: Source # (VonMises, VonMises) -> [(Double,Double)]
ellipse vm =
    let [mux,kpx,muy,kpy] = listCoordinates vm
        f t = (mux + sqrt (recip kpx) * cos t, muy + sqrt (recip kpy) * sin t)
     in f <$> range 0 (2*pi) 100

mixtureModelToConfidenceCSV :: Natural # Mixture (VonMises,VonMises) 2 -> [[(Double,Double)]]
mixtureModelToConfidenceCSV mxmdl =
    let cmps = S.toList . fst $ splitNaturalMixture mxmdl
     in ellipse . toSource <$> cmps

mixtureModelToMeanCSV
    :: Natural # Mixture (VonMises,VonMises) 2 -> [(Double,Double)]
mixtureModelToMeanCSV mxmdl =
    let cmps = S.toList . fst $ splitNaturalMixture mxmdl
     in [ (mu1,mu2) | [mu1,_,mu2,_] <- listCoordinates . toSource <$> cmps ]


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps trumxmdl
    mxmdl0 <- realize $ initialize mxmdlint

    let xys = fst <$> cxys

    let emmxmdls = take nepchs $ iterate (vonMisesEM xys) mxmdl0
        --itrmxmdl1 :: Natural # Mixture VonMises 2
        --(itrmxmdl1,itrlls) = iterativeMixtureOptimization
        --    (div nepchs 3) eps defaultAdamPursuit Nothing 0.1 xys $ toNatural vm2'

    let admmxmdls = take nepchs . takeEvery admmlt
            $ vanillaGradientSequence (logLikelihoodDifferential xys) eps defaultAdamPursuit mxmdl0

    let trunlls = repeat . average $ negate . logObservableDensity trumxmdl <$> xys
        emnlls = [ average $ negate . logObservableDensity mxmdl <$> xys | mxmdl <- emmxmdls ]
        admnlls = [ average $ negate . logObservableDensity mxmdl <$> xys | mxmdl <- admmxmdls ]
        --itrnlls = negate <$> concat itrlls

    let cedcsvs = zip3 trunlls emnlls admnlls

    let emmxmdl1 = last emmxmdls
        admmxmdl1 = last admmxmdls
        [trucnfs,emcnfs,admcnfs] = mixtureModelToConfidenceCSV <$> [trumxmdl,emmxmdl1,admmxmdl1]

    let [trumn,emmn,admmn] = mixtureModelToMeanCSV <$> [trumxmdl,emmxmdl1,admmxmdl1]

    let ldpth = "data"

    goalExport ldpth "log-likelihood" cedcsvs

    print $ length trucnfs

    goalExportLines ldpth "true-confidence" trucnfs
    goalExportLines ldpth "em-confidence" emcnfs
    goalExportLines ldpth "adam-confidence" admcnfs

    goalExport ldpth "true-mean" trumn
    goalExport ldpth "em-mean" emmn
    goalExport ldpth "adam-mean" admmn

    goalExport ldpth "samples" xys

    runGnuplot ldpth "cross-entropy-descent"
    runGnuplot ldpth "mixture-components"

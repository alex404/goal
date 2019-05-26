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

kpx1,kpx2,kpx3,kpy1,kpy2,kpy3 :: Double
kpx1 = 2
kpy1 = 2
kpx2 = 2
kpy2 = 2
kpx3 = 2
kpy3 = 2

vm1,vm2,vm3 :: Source # (VonMises,VonMises)
vm1 = Point $ S.fromTuple (mux1, kpx1, muy1, kpy1)
vm2 = Point $ S.fromTuple (mux2, kpx2, muy2, kpy2)
vm3 = Point $ S.fromTuple (mux3, kpx3, muy3, kpy3)

vms :: S.Vector 3 (Natural # (VonMises,VonMises))
vms = S.fromTuple (toNatural vm1,toNatural vm2,toNatural vm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Categorical Int 2
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Mixture (VonMises,VonMises) 2
truhrm = joinMixture vms $ toNatural wghts

-- Mixture Distributions --

vm1',vm2',vm3' :: Source # (VonMises,VonMises)
vm1' = Point $ S.fromTuple (2, 0.5, pi, 0.5)
vm2' = Point $ S.fromTuple (3, 0.5, pi, 0.5)
vm3' = Point $ S.fromTuple (4, 0.5, pi, 0.5)

vms' :: S.Vector 3 (Natural # (VonMises,VonMises))
vms' = S.fromTuple (toNatural vm1',toNatural vm2',toNatural vm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Categorical Int 2
wghts' = Point $ S.doubleton mix1' mix2'

hrm0 :: Natural # Mixture (VonMises,VonMises) 2
hrm0 = joinMixture vms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 50

eps :: Double
eps = -0.05

bnd :: Double
bnd = 1e-5

admmlt :: Int
admmlt = 10

-- EM
nepchs :: Int
nepchs = 200

-- Functions --


--vonMisesEM
--    :: Sample (VonMises,VonMises) -- ^ Observations
--    -> (Natural # Mixture (VonMises,VonMises) 2, S.Vector 3 (Natural # (VonMises,VonMises)))
--    -> (Natural # Mixture (VonMises,VonMises) 2, S.Vector 3 (Natural # (VonMises,VonMises)))
--vonMisesEM zs (hrm,nzs0) =
--    let zs' = hSingleton <$> zs
--        (cats,mzs) = deepMixtureExpectationStep zs' $ transposeHarmonium hrm
--        nzs = S.zipWith cauchyFun (S.map fromOneHarmonium mzs) nzs0
--     in (joinMixture nzs cats, nzs)
--    where diffFun mz nz = joinTangentPair nz $ crossEntropyDifferential mz nz
--          cauchyFun mz nz = cauchyLimit euclideanDistance bnd
--              $ vanillaGradientSequence (diffFun mz) eps defaultAdamPursuit nz

vonMisesEM
    :: Sample (VonMises,VonMises) -- ^ Observations
    -> Natural # Mixture (VonMises,VonMises) 2
    -> Natural # Mixture (VonMises,VonMises) 2
vonMisesEM zs nhrm =
    cauchyFun $ harmoniumEmpiricalExpectations zs nhrm
    where diffFun mhrm' = pairTangentFunction $ crossEntropyDifferential mhrm'
          cauchyFun mhrm' = cauchyLimit euclideanDistance bnd
              $ vanillaGradientSequence (diffFun mhrm') eps defaultAdamPursuit nhrm

filterCat :: [SamplePoint (Mixture (VonMises,VonMises) 2)] -> Int -> [SamplePoint (VonMises,VonMises)]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

ellipse :: Source # (VonMises, VonMises) -> [ConfidenceEllipse]
ellipse vm =
    let [mux,kpx,muy,kpy] = listCoordinates vm
        f t = ConfidenceEllipse (mux + sqrt (recip kpx) * cos t) (muy + sqrt (recip kpy) * sin t)
     in f <$> range 0 (2*pi) 100

mixtureModelToConfidenceCSV :: Natural # Mixture (VonMises,VonMises) 2 -> [[ConfidenceEllipse]]
mixtureModelToConfidenceCSV hrm =
    let cmps = S.toList . fst $ splitMixture hrm
     in ellipse . toSource <$> cmps

mixtureModelToMeanCSV :: Natural # Mixture (VonMises,VonMises) 2 -> [VonMisesMeans]
mixtureModelToMeanCSV hrm =
    let cmps = S.toList . fst $ splitMixture hrm
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

expmnt :: Experiment
expmnt = Experiment "probability" "von-mises-mixture"


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps truhrm

    let xys = hHead <$> cxys

    let emhrms = take nepchs $ iterate (vonMisesEM xys) hrm0
        --itrhrm1 :: Natural # Mixture VonMises 2
        --(itrhrm1,itrlls) = iterativeMixtureOptimization
        --    (div nepchs 3) eps defaultAdamPursuit Nothing 0.1 xys $ toNatural vm2'

    let sgd = pairTangentFunction $ stochasticMixtureDifferential xys
        admhrms = take nepchs . takeEvery admmlt
            $ vanillaGradientSequence sgd eps defaultAdamPursuit hrm0

    let trunlls = repeat . average $ negate . log . mixtureDensity truhrm <$> xys
        emnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- emhrms ]
        admnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- admhrms ]
        --itrnlls = negate <$> concat itrlls

    let cedcsvs = zipWith3 CrossEntropyDescent trunlls emnlls admnlls

    let emhrm1 = last emhrms
        admhrm1 = last admhrms
        (cnfcsv:cnfcsvs) = concat $ mixtureModelToConfidenceCSV <$> [truhrm,emhrm1,admhrm1]

    let mncsvs = mixtureModelToMeanCSV <$> [truhrm,emhrm1,admhrm1]

    let xycsv = [ TrainingSamples x y | (x,y) <- xys ]

    goalExportNamed True expmnt Nothing cedcsvs

    goalExportNamed False expmnt Nothing cnfcsv

    mapM_ (goalExportNamed False expmnt Nothing) cnfcsvs

    mapM_ (goalExportNamed False expmnt Nothing) mncsvs

    goalExportNamed False expmnt Nothing xycsv

    runGnuplot expmnt Nothing defaultGnuplotOptions "cross-entropy-descent.gpi"
    runGnuplot expmnt Nothing defaultGnuplotOptions "mixture-components.gpi"

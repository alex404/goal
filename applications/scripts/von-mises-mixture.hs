{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import Data.List


--- Globals ---

prjdr :: String
prjdr = "probability/von-mises-mixture"

-- Manifolds --

type Latent = Categorical Int 2
type Observable = (VonMises, VonMises)
type Harmonium' = Harmonium Tensor Observable Latent

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 1.5
muy1 = 1.5
mux2 = 3
muy2 = 3
mux3 = 4.5
muy3 = 4.5

kpx1,kpx2,kpx3,kpy1,kpy2,kpy3 :: Double
kpx1 = 2
kpy1 = 2
kpx2 = 2
kpy2 = 2
kpx3 = 2
kpy3 = 2

vm1,vm2,vm3 :: Source # Observable
vm1 = Point $ S.fromTuple (mux1, kpx1, muy1, kpy1)
vm2 = Point $ S.fromTuple (mux2, kpx2, muy2, kpy2)
vm3 = Point $ S.fromTuple (mux3, kpx3, muy3, kpy3)

vms :: S.Vector 3 (Natural # Observable)
vms = S.fromTuple (toNatural vm1,toNatural vm2,toNatural vm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Harmonium Tensor Observable Latent
truhrm = buildMixtureModel vms $ toNatural wghts

-- Mixture Distributions --

vm1',vm2',vm3' :: Source # Observable
vm1' = Point $ S.fromTuple (2, 0.5, pi, 0.5)
vm2' = Point $ S.fromTuple (3, 0.5, pi, 0.5)
vm3' = Point $ S.fromTuple (4, 0.5, pi, 0.5)

vms' :: S.Vector 3 (Natural # Observable)
vms' = S.fromTuple (toNatural vm1',toNatural vm2',toNatural vm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Latent
wghts' = Point $ S.doubleton mix1' mix2'

hrm0 :: Natural # Harmonium Tensor Observable Latent
hrm0 = buildMixtureModel vms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 200

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

vonMisesEM
    :: Sample Observable -- ^ Observations
    -> (Natural # Harmonium', S.Vector 3 (Natural # Observable))
    -> (Natural # Harmonium', S.Vector 3 (Natural # Observable))
vonMisesEM zs (hrm,nzs0) =
    let zs' = hSingleton <$> zs
        (cats,mzs) = deepMixtureModelExpectationStep zs' $ transposeHarmonium hrm
        nzs = S.zipWith cauchyFun (S.map fromOneHarmonium mzs) nzs0
     in (buildMixtureModel nzs cats, nzs)
    where diffFun mz nz = joinTangentPair nz $ crossEntropyDifferential mz nz
          cauchyFun mz nz = cauchyLimit euclideanDistance bnd
              $ vanillaGradientSequence (diffFun mz) eps defaultAdamPursuit nz

filterCat :: [SamplePoint Harmonium'] -> Int -> [SamplePoint Observable]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

ellipse :: Source # (VonMises, VonMises) -> [ConfidenceEllipse]
ellipse vm =
    let [mux,kpx,muy,kpy] = listCoordinates vm
        f t = ConfidenceEllipse (mux + sqrt (recip kpx) * cos t) (muy + sqrt (recip kpy) * sin t)
     in f <$> range 0 (2*pi) 100

mixtureModelToConfidenceCSV :: Natural # Harmonium Tensor Observable Latent -> [[ConfidenceEllipse]]
mixtureModelToConfidenceCSV hrm =
    let cmps = S.toList . fst $ splitMixtureModel hrm
     in ellipse . toSource <$> cmps

mixtureModelToMeanCSV :: Natural # Harmonium Tensor Observable Latent -> [VonMisesMeans]
mixtureModelToMeanCSV hrm =
    let cmps = S.toList . fst $ splitMixtureModel hrm
     in [ VonMisesMeans mu1 mu2 | [mu1,_,mu2,_] <- listCoordinates . toSource <$> cmps ]

-- CSV --

data CrossEntropyDescent = CrossEntropyDescent
    { trueCrossEntropy :: Double
    , emCrossEntropy :: Double
    , itrCrossEntropy :: Double
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

    let xys = hHead <$> cxys

    let emhrms = take nepchs $ fst <$> iterate (vonMisesEM xys) (hrm0, vms')
        itrhrm1 :: Natural # Harmonium'
        (itrhrm1,itrlls) = iterativeMixtureModelOptimization
            (div nepchs 3) eps defaultAdamPursuit Nothing (0.1) xys $ toNatural vm2'

    let sgd hrm = joinTangentPair hrm $ stochasticMixtureModelDifferential xys hrm
        admhrms = take nepchs . takeEvery admmlt
            $ vanillaGradientSequence sgd eps defaultAdamPursuit hrm0

    let emhrm1 = last emhrms
        admhrm1 = last admhrms
        (cnfcsv:cnfcsvs) = concat $ mixtureModelToConfidenceCSV <$> [truhrm,emhrm1,itrhrm1,admhrm1]

    goalWriteNamedAnalysis "probability" "von-mises-mixture" "mixture-components" Nothing cnfcsv

    mapM_ (goalAppendNamedAnalysis "probability" "von-mises-mixture" "mixture-components" Nothing) cnfcsvs

    let mncsvs = mixtureModelToMeanCSV <$> [truhrm,emhrm1,admhrm1]

    mapM_ (goalAppendNamedAnalysis "probability" "von-mises-mixture" "mixture-components" Nothing) mncsvs

    let xycsv = [ TrainingSamples x y | (x,y) <- xys ]

    goalAppendNamedAnalysis "probability" "von-mises-mixture" "mixture-components" Nothing xycsv

    let trunlls = repeat . average $ negate . log . mixtureDensity truhrm <$> xys
        emnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- emhrms ]
        admnlls = [ average $ negate . log . mixtureDensity hrm <$> xys | hrm <- admhrms ]
        itrnlls = negate <$> concat itrlls

    let cedcsvs = zipWith4 CrossEntropyDescent trunlls emnlls itrnlls admnlls

    goalWriteNamedAnalysis "probability" "von-mises-mixture" "cross-entropy-descent" Nothing cedcsvs


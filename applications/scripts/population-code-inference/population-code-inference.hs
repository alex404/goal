{-# LANGUAGE DeriveGeneric,DataKinds,TypeOperators,TypeFamilies,FlexibleContexts #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import Goal.Geometry
import Goal.Probability

-- Qualified --

import qualified Data.List as L

--- Program ---


-- CSV --

expmnt :: Experiment
expmnt = Experiment "probability" "population-code-inference"

data PopulationCodeInference = PopulationCodeInference
    { radians :: Double
    , conjugatePosterior :: Double
    , numericalPosterior :: Double }
    deriving (Generic, Show)

instance FromNamedRecord PopulationCodeInference
instance ToNamedRecord PopulationCodeInference
instance DefaultOrdered PopulationCodeInference
instance NFData PopulationCodeInference


-- Globals --

x0 :: Double
x0 = pi+1

mcts :: Source # Categorical Int 2
mcts = Point $ S.doubleton 0.5 0.2

sp0 :: Source # VonMises
sp0 = Point $ S.doubleton pi 10

sp1 :: Source # VonMises
sp1 = Point $ S.doubleton 1 10

sp2 :: Source # VonMises
sp2 = Point $ S.doubleton 5 10

prr :: Natural # Harmonium Tensor VonMises (Categorical Int 2)
prr = joinMixtureModel (S.map toNatural $ S.fromTuple (sp0,sp1,sp2)) (toNatural mcts)

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

type NNeurons = 6

type Neurons = R NNeurons Poisson

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mnx mxx

xsmps :: [Double]
xsmps = init $ range mnx mxx 50

pltsmps :: [Double]
pltsmps = init $ range mnx mxx 500

kp :: Double
kp = 2

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map (Point . flip S.doubleton kp) mus

gn0 :: Double
gn0 = log 0.5

lkl0 :: Mean #> Natural # R NNeurons Poisson <* VonMises
lkl0 = joinVonMisesPopulationEncoder (Left gn0) $ S.map toNatural sps

rho0 :: Double
rprms0 :: Natural # VonMises
(rho0,rprms0) = regressConjugationParameters lkl0 xsmps

rprms1,rprms2 :: Natural # VonMises
rprms1 = Point $ S.doubleton 1 0
rprms2 = Point $ S.doubleton 2 0

lkl1,lkl2 :: Mean #> Natural # R NNeurons Poisson <* VonMises
lkl1 = conjugatePopulationEncoder rho0 rprms1 xsmps lkl0
lkl2 = conjugatePopulationEncoder rho0 rprms2 xsmps lkl0

numericalPosteriorFunction :: [Response NNeurons] -> Double -> Double
numericalPosteriorFunction zs =
    let xs = range mnx mxx 100
     in fst $ numericalRecursiveBayesianInference 1e-6 mnx mxx xs [lkl0,lkl1,lkl2] zs (mixtureDensity prr)


-- Main --


main :: IO ()
main = do

    z0 <- realize . samplePoint $ lkl0 >.>* x0
    z1 <- realize . samplePoint $ lkl1 >.>* x0
    z2 <- realize . samplePoint $ lkl2 >.>* x0

    let zs = [z0,z1,z2]
        zcsvs = L.transpose $ S.toList mus : fmap (map fromIntegral . B.toList) zs

    let pst1 = conjugatedBayesRule rprms0 lkl0 z0 prr
        pst2 = conjugatedBayesRule rprms1 lkl1 z1 pst1
        pst3 = conjugatedBayesRule rprms2 lkl2 z2 pst2

        pst0' = numericalPosteriorFunction []
        pst1' = numericalPosteriorFunction $ take 1 zs
        pst2' = numericalPosteriorFunction $ take 2 zs
        pst3' = numericalPosteriorFunction zs

        pstcsvs pst pst' = zipWith3 PopulationCodeInference pltsmps (mixtureDensity pst <$> pltsmps) (pst' <$> pltsmps)

    goalExport True expmnt Nothing zcsvs

    goalExportNamed False expmnt Nothing $ pstcsvs prr pst0'
    goalExportNamed False expmnt Nothing $ pstcsvs pst1 pst1'
    goalExportNamed False expmnt Nothing $ pstcsvs pst2 pst2'
    goalExportNamed False expmnt Nothing $ pstcsvs pst3 pst3'

    runGnuplot expmnt Nothing defaultGnuplotOptions "samples.gpi"
    runGnuplot expmnt Nothing defaultGnuplotOptions "posteriors.gpi"

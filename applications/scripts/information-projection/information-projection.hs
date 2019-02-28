{-# LANGUAGE
    ScopedTypeVariables,
    TypeFamilies,
    TypeOperators,
    FlexibleContexts,
    DataKinds,
    DeriveGeneric,
    Arrows
    #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Data.List

--- Globals ---


-- Manifolds --

type Latent = Categorical Int 2

-- Mixture Distributions --

mu1,mu2,mu3 :: Double
mu1 = -4
mu2 = 0.5
mu3 = 3

nrm1,nrm2,nrm3 :: Source # Normal
nrm1 = Point $ S.doubleton mu1 0.6
nrm2 = Point $ S.doubleton mu2 0.4
nrm3 = Point $ S.doubleton mu3 0.5

nrms :: S.Vector 3 (Natural # Normal)
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.3
mix2 = 0.3

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

hrm :: Natural # Harmonium Tensor Normal Latent
hrm = joinMixtureModel nrms $ toNatural wghts

-- Training --

sx0 :: Source # Normal
sx0 = Point $ S.doubleton 0 20

nx0 :: Natural # Normal
nx0 = transition sx0

ceeps,ipeps :: Double
ceeps = -0.005
ipeps = -0.01

nbtch :: Int
nbtch = 50

nipsmps :: Int
nipsmps = 10

-- Plotting --

pltsmps :: [Double]
pltsmps = range (-7) 7 100

-- CSV --

data InformationProjection = InformationProjection
    { xAxis :: Double
    , trueDensity :: Double
    , ceDensity :: Double
    , ipDensity :: Double }
    deriving (Generic, Show)

instance FromNamedRecord InformationProjection
instance ToNamedRecord InformationProjection
instance DefaultOrdered InformationProjection
instance NFData InformationProjection

expmnt :: Experiment
expmnt = Experiment "probability" "information-projection"

--- Main ---


main :: IO ()
main = do

    tcxs <- realize $ sample nbtch hrm
    let txs = hHead <$> tcxs

    let cechn = chainCircuit nx0 $ proc nx -> do
            dnx <- arr $ stochasticCrossEntropyDifferential txs -< nx
            gradientCircuit ceeps defaultAdamPursuit -< breakPoint $ joinTangentPair nx dnx

    let ipchn = chainCircuit nx0 $ proc nx -> do
            dnx <- arrM $ harmoniumInformationProjectionDifferential nipsmps (transposeHarmonium hrm) -< nx
            gradientCircuit ipeps defaultAdamPursuit -< breakPoint $ joinTangentPair nx dnx

    cenx <- realize $ iterateChain 1000 cechn
    ipnx <- realize $ iterateChain 1000 ipchn

    let trusmps = mixtureDensity hrm <$> pltsmps
        cesmps = density cenx <$> pltsmps
        ipsmps = density ipnx <$> pltsmps

    let csv = zipWith4 InformationProjection pltsmps trusmps cesmps ipsmps

    goalExportNamed True expmnt Nothing csv
    runGnuplot expmnt Nothing defaultGnuplotOptions "information-projection.gpi"

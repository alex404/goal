#! stack runghc

{-# LANGUAGE
    DeriveGeneric,
    UndecidableInstances,
    StandaloneDeriving,
    ScopedTypeVariables,
    DataKinds,
    TypeOperators,
    TypeFamilies,
    FlexibleContexts
    #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

-- Unqualified --

import qualified Data.List as L


--- Project ---


data Univariate m = Univariate
    { sampleValue :: SamplePoint m
    , trueDensity :: Double
    , sourceFit :: Double
    , exponentialFamilyFit :: Double }

deriving instance Show (SamplePoint m) => Show (Univariate m)
deriving instance Generic (Univariate m)

instance ToField (SamplePoint m) => ToNamedRecord (Univariate m) where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered (Univariate m) where
    headerOrder = goalCSVOrder


--- Globals ---

nsmps :: Int
nsmps = 20

-- Binomial --

ttlB :: String
ttlB = "binomial"

mnB,mxB :: Double
(mnB,mxB) = (0,5)

bnsB :: Int
bnsB = 11

truB :: Source # Binomial 5
truB = singleton 0.3

rngB :: [Int]
rngB = pointSampleSpace truB

-- Categorical --

ttlC :: String
ttlC = "categorical"

mnC,mxC :: Double
(mnC,mxC) = (0,4)

bnsC :: Int
bnsC = 5

truC :: Source # Categorical 4
truC = fromTuple (0.1,0.4,0.1,0.2)

rngC :: [Int]
rngC = pointSampleSpace truC

-- Poisson --

ttlP :: String
ttlP = "poisson"

mnP,mxP :: Double
(mnP,mxP) = (0,20)

bnsP :: Int
bnsP = 20

truP :: Source # Poisson
truP = singleton 5

rngP :: [Int]
rngP = [0..20]

-- Normal --

ttlN :: String
ttlN = "normal"

mnN,mxN :: Double
(mnN,mxN) = (-3,7)

bnsN :: Int
bnsN = 20

truN :: Source # Normal
truN = fromTuple (2,0.7)

rngN :: [Double]
rngN = range (-3) 7 100

-- Layout --

generateLayout
    :: forall m
    . ( Transition Source Mean m, Transition Source Natural m, LegendreExponentialFamily m
      , MaximumLikelihood Source m, AbsolutelyContinuous Source m, Generative Source m
      , AbsolutelyContinuous Natural m, ToField (SamplePoint m), Real (SamplePoint m) )
    => String
    -> Int
    -> Double
    -> Double
    -> Sample m
    -> Source # m
    -> IO ()
generateLayout ttl nb mn mx rng tru = do

    smps <- realize $ sample nsmps tru

    let smle :: Source # m
        smle = mle smps
        nmle = toNatural smle

    let tdns = density tru <$> rng
        sdns = density smle <$> rng
        ndns = density nmle <$> rng

    let univs :: [Univariate m]
        univs = L.zipWith4 Univariate rng tdns sdns ndns

    let (bns,_,[wghts]) = histograms nb (Just (mn,mx)) [realToFrac <$> smps]

    goalExportNamed ttl "densities" univs

    goalExport ttl "histogram" $ zip bns wghts

    runGnuplot ttl "univariate"

main :: IO ()
main = do

    generateLayout ttlB bnsB mnB mxB rngB truB
    generateLayout ttlC bnsC mnC mxC rngC truC
    generateLayout ttlP bnsP mnP mxP rngP truP
    generateLayout ttlN bnsN mnN mxN rngN truN

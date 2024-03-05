{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Aeson (decode)
import Data.ByteString.Lazy qualified as B
import Data.Maybe (fromJust)
import Data.Proxy (Proxy (..))

--- Globals ---

--- Data loading

-- Function to decode the JSON file into the MNISTSizeData type
decodeMNISTSize :: FilePath -> IO (Maybe [[Int]])
decodeMNISTSize filePath = do
    content <- B.readFile filePath
    -- Directly parse the content into a Map String [[Int]]
    let mnistData = decode content :: Maybe [[Int]]
    return mnistData

--- Type Synonyms

type ObservationSpace = 2
type MNISTSizeRows = 28
type MNISTSize = 784
type BoltzmannNeurons = 9

--- Linear Population

mumn, mumx, stp :: Double
mumn = -5
mumx = 5
stp = (mumx - mumn) / (fromIntegral (natValInt $ Proxy @MNISTSizeRows) - 1)

mu0s :: S.Vector MNISTSizeRows Double
mu0s = S.enumFromStepN mumn stp

mus :: S.Vector MNISTSize (Source # StandardNormal 2)
mus = S.concatMap (\x -> S.map (\y -> fromTuple (x, y)) mu0s) mu0s

var :: Double
var = 0.1

tcs :: S.Vector MNISTSize (Source # FullNormal 2)
tcs = S.map (`join` fromTuple (var, 0, var)) mus

gn :: Double
gn = 0.5

gns :: Source # Replicated MNISTSize Poisson
gns = Point $ S.replicate gn

ppc :: Natural # Replicated MNISTSize Poisson <* FullNormal 2
ppc = joinPopulationCode (toNatural gns) (S.map toNatural tcs)

--- Conjugation Parameter Regression

regmn, regmx :: Double
regmn = mumn + 1
regmx = mumx - 1
regres :: Int
regres = 200

regxys :: Sample (FullNormal 2)
regxys = [S.fromTuple (x, y) | x <- range regmn regmx regres, y <- range regmn regmx regres]

chixyht :: Double
rhoxyht :: Natural # FullNormal 2
(chixyht, rhoxyht) = conjugationParameterRegression regxys ppc

--- Initialization

unibnd :: Double
unibnd = 4

--- Training

eps :: Double
eps = 3e-3

nstps, nsmps, nepchs :: Int
nstps = 2000
nsmps = 100
nepchs = 5

gp :: GradientPursuit
gp = defaultAdamPursuit

stochasticEMStep ::
    [S.Vector MNISTSize Int] ->
    (Int, Natural # GaussianBoltzmannPopulationCode L.PositiveDefinite MNISTSize ObservationSpace BoltzmannNeurons) ->
    IO (Int, Natural # GaussianBoltzmannPopulationCode L.PositiveDefinite MNISTSize ObservationSpace BoltzmannNeurons)
stochasticEMStep imgs (k, gbppc) = do
    gbppc' <- realize . iterateChain nstps $ ppcExpectationMaximization imgs eps nsmps gp rhoxyht gbppc
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show $ ppcLogLikelihood (chixyht, rhoxyht) imgs gbppc'
            ]
    return (k + 1, gbppc')

--- Plotting

pltsmps :: Int
pltsmps = 100

pltmn, pltmx :: Double
pltmn = mumn - 2
pltmx = mumx + 2

pltxys :: Sample (DiagonalNormal 2)
pltxys = [S.fromTuple (x, y) | x <- range pltmn pltmx pltsmps, y <- range pltmn pltmx pltsmps]

--- Main ---

main :: IO ()
main = do
    --- Load MNIST data
    jsnpth <- dataFilePath "mnist-compressed.json"
    mnistData <- decodeMNISTSize jsnpth
    imgs0 <- case mnistData of
        Just mnsts -> return mnsts
        Nothing -> error "Failed to decode MNISTSize JSON data."

    let imgs :: [S.Vector MNISTSize Int]
        imgs = fromJust . S.fromList <$> imgs0

    print $ head imgs

    --- Initialize Gaussian-Boltzmann machine
    bltz0 :: Natural # Boltzmann BoltzmannNeurons <- realize $ uniformInitialize (-unibnd, unibnd)

    shfts0 :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated BoltzmannNeurons Bernoulli) <-
        realize $ uniformInitialize (-unibnd, unibnd)

    --- Linear Model
    let mvn0 :: Natural # FullNormal ObservationSpace
        mvn0 = standardNormal

        lmdl0 :: Natural # BoltzmannLinearModel L.PositiveDefinite ObservationSpace BoltzmannNeurons
        lmdl0 = join mvn0 shfts0

    --- Harmonium
    let gbhrm0 :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace BoltzmannNeurons
        gbhrm0 = joinConjugatedHarmonium lmdl0 bltz0

    --- GBPPC

    let gbppc0 ::
            Natural
                # GaussianBoltzmannPopulationCode
                    L.PositiveDefinite
                    MNISTSize
                    ObservationSpace
                    BoltzmannNeurons
        gbppc0 = approximateJoinConjugatedHarmonium rhoxyht ppc gbhrm0

    --- Training
    -- putStrLn $ "Initial Log-Likelihood: " ++ show (ppcLogLikelihood (chixyht, rhoxyht) imgs gbppc0)
    kgbppcs <- iterateM nepchs (stochasticEMStep imgs) (1, gbppc0)

    let gbppcs = snd <$> kgbppcs
        gbhrms = snd . approximateSplitConjugatedHarmonium rhoxyht <$> gbppcs
        lrngbhrm = last gbhrms

    --- Lines
    let dns0 = observableDensities gbhrm0 pltxys
        lrndns = observableDensities lrngbhrm pltxys

    let json =
            toJSON
                [ "xys" .= pltxys
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                ]

    --- Process data
    flnm <- resultsFilePath "mnist-ppc.json"

    exportJSON flnm json

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

type N = 28
type NN = N * N
type BN = 9
type OS = 2
type Neurons n = Replicated n Poisson

--- Linear Population

mumn, mumx, stp :: Double
mumn = -6
mumx = 6
stp = (mumx - mumn) / (fromIntegral (natValInt $ Proxy @N) - 1)

mu0s :: S.Vector N Double
mu0s = S.enumFromStepN mumn stp

mus :: S.Vector NN (Source # StandardNormal 2)
mus = S.concatMap (\x -> S.map (\y -> fromTuple (x, y)) mu0s) mu0s

var :: Double
var = 0.15

tcs :: S.Vector NN (Source # FullNormal 2)
tcs = S.map (`join` fromTuple (var, 0, var)) mus

gns :: Source # Neurons NN
gns = Point $ S.replicate 8

lkl :: Natural # Neurons NN <* FullNormal 2
lkl = joinPopulationCode (toNatural gns) (S.map toNatural tcs)

--- Regression

regbfr :: Double
regbfr = 1

regmn, regmx :: Double
regmn = mumn + regbfr
regmx = mumx - regbfr

regres :: Int
regres = 200

regxys :: Sample (FullNormal 2)
regxys = [S.fromTuple (x, y) | x <- range regmn regmx regres, y <- range regmn regmx regres]

chixyht :: Double
rhoxyht :: Natural # FullNormal 2
(chixyht, rhoxyht) = conjugationParameterRegression regxys lkl

regptns, regptnsht, regptndffs :: [Double]
regptns = potential <$> lkl >$>* regxys
regptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> ptnpltxys
regptndffs = zipWith (-) regptns regptnsht

--- Initialization

unibnd :: Double
unibnd = 0.001

initializeGaussianBoltzmannHarmonium ::
    Random (Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN)
initializeGaussianBoltzmannHarmonium = do
    bltz :: Natural # Boltzmann BN <- uniformInitialize (-unibnd, unibnd)
    shfts :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli) <-
        uniformInitialize (-unibnd, unibnd)

    --- Linear Model
    let mvn = standardNormal
        lmdl = join mvn shfts

    return $ joinConjugatedHarmonium lmdl bltz

bltzhnd :: Natural # Boltzmann BN
bltzhnd = join (fromTuple (1, 1, 1, 1, 0, 1, 1, 1, 1)) $ -10

mubnd :: Double
mubnd = 4

shftshnd :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli)
shftshnd =
    fromRows $
        S.fromTuple
            ( fromTuple (-mubnd, -mubnd, -mubnd, 0, 0, 0, mubnd, mubnd, mubnd)
            , fromTuple (-mubnd, 0, mubnd, -mubnd, 0, mubnd, -mubnd, 0, mubnd)
            )

smvnhnd :: Source # FullNormal OS
smvnhnd = standardNormal / 2

mvnhnd :: Natural # FullNormal OS
mvnhnd = toNatural smvnhnd

lmdlhnd :: Natural # BoltzmannLinearModel L.PositiveDefinite OS BN
lmdlhnd = join mvnhnd shftshnd

gbhrmhnd :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN
gbhrmhnd = joinConjugatedHarmonium lmdlhnd bltzhnd

--- Training

eps :: Double
eps = 3e-3

nstps, nepchs :: Int
nstps = 1000
nepchs = 10

nbtch :: Int
nbtch = 20

gp :: GradientPursuit
gp = defaultAdamPursuit

loggingEMStep ::
    (LinearSubspace y (FullNormal 2), LegendreExponentialFamily y, Generative Natural y) =>
    [Natural # FullNormal 2] ->
    Sample (Replicated NN Poisson) ->
    (Int, Natural # y) ->
    IO (Int, Natural # y)
loggingEMStep ny0s imgs (k, nltnt) = do
    let nltnt' = ppcExpectationMaximizationAscent eps gp ny0s nltnt !! nstps
    -- nltnt' <- realize . iterateChain nstps $ ppcStochasticMaximumLikelihood eps gp nbtch ny0s nltnt
    let ppc' = approximateJoinConjugatedHarmonium rhoxyht lkl nltnt'
    putStrLn
        . concat
        $ [ "\nIteration: "
          , show k
          , "\nLog-Likelihood: "
          , show $ ppcLogLikelihood (chixyht, rhoxyht) imgs ppc'
          ]
    return (k + 1, nltnt')

--- Plotting

pltres :: Int
pltres = 100

ptnpltmn, ptnpltmx :: Double
ptnpltmn = regmn - 2
ptnpltmx = regmx + 2

ptnpltxys :: Sample (FullNormal 2)
ptnpltxys = [S.fromTuple (x, y) | x <- range ptnpltmn ptnpltmx pltres, y <- range ptnpltmn ptnpltmx pltres]

pltptns, pltptnsht, pltptndffs :: [Double]
pltptns = potential <$> lkl >$>* ptnpltxys
pltptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> ptnpltxys
pltptndffs = zipWith (-) pltptns pltptnsht

dnspltmn, dnspltmx :: Double
dnspltmn = regmn
dnspltmx = regmx

dnspltxys :: Sample (FullNormal 2)
dnspltxys = [S.fromTuple (x, y) | x <- range dnspltmn dnspltmx pltres, y <- range dnspltmn dnspltmx pltres]

-- Note, these copies seem required to avoid segfaulting!?!?! I should probably report this as a bug.
dnspltxys2 :: Sample (FullNormal 2)
dnspltxys2 = [S.fromTuple (x, y) | x <- range dnspltmn dnspltmx pltres, y <- range dnspltmn dnspltmx pltres]

--- Main ---

main :: IO ()
main = do
    print ("\nPotential regression RMSE:" :: String)
    print . sqrt . average $ square <$> regptndffs

    --- Load MNIST data
    jsnpth <- dataFilePath "mnist-compressed.json"
    mnistData <- decodeMNISTSize jsnpth
    imgs0 <- case mnistData of
        Just mnsts -> return mnsts
        Nothing -> error "Failed to decode MNISTSize JSON data."

    let imgs :: [S.Vector NN Int]
        imgs = take 1 $ fromJust . S.fromList <$> imgs0

    mapM_ print . breakEvery 28 . S.toList $ head imgs

    --- Initialization
    -- ltnt0 <- realize initializeGaussianBoltzmannHarmonium
    let ltnt0 = gbhrmhnd
    let ppc0 = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural ltnt0
        dns0 = observableDensities ltnt0 dnspltxys
    print $ sum dns0

    --- Training
    putStrLn $ "Initial Log-Likelihood: " ++ show (ppcLogLikelihood (chixyht, rhoxyht) imgs ppc0)

    let ny0s = ppcExpectationBiases imgs lkl
    kltnts <- iterateM nepchs (loggingEMStep ny0s imgs) (1, ltnt0)
    putStrLn ("\nTraining Complete" :: String)

    let ltnts = snd <$> kltnts
        lrndns = observableDensities (last ltnts) dnspltxys2
    print $ sum lrndns

    let jsonData =
            toJSON
                [ "tuning-curve-xys" .= ptnpltxys
                , "preferred-stimuli" .= S.map coordinates mus
                , "sum-of-tuning-curves" .= pltptns
                , "estimated-sum-of-tuning-curves" .= pltptnsht
                , "estimation-difference" .= pltptndffs
                , "regression-bounds" .= (regmn, regmx)
                , "density-xys" .= dnspltxys
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                ]

    rsltfl <- resultsFilePath "mnist-ppc.json"
    exportJSON rsltfl jsonData
    putStrLn "\nSimulation Complete\n"

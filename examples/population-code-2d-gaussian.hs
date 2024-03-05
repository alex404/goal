{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MonoLocalBinds #-}
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

--- Misc

import Data.Proxy (Proxy (..))

--- Globals ---

--- Types

type N = 28
type NN = N * N

type Neurons n = Replicated n Poisson

--- Linear Population

mumn, mumx, stp :: Double
mumn = -5
mumx = 5
stp = (mumx - mumn) / (fromIntegral (natValInt $ Proxy @N) - 1)

mu0s :: S.Vector N Double
mu0s = S.enumFromStepN mumn stp

mus :: S.Vector NN (Source # StandardNormal 2)
mus = S.concatMap (\x -> S.map (\y -> fromTuple (x, y)) mu0s) mu0s

var :: Double
var = 0.1

tcs :: S.Vector NN (Source # FullNormal 2)
tcs = S.map (`join` fromTuple (var, 0, var)) mus

gns :: Source # Neurons NN
gns = Point $ S.replicate 2

lkl :: Natural # Neurons NN <* FullNormal 2
lkl = joinPopulationCode (toNatural gns) (S.map toNatural tcs)

--- Regression

regmn, regmx :: Double
regmn = mumn + 1
regmx = mumx - 1
regres :: Int
regres = 200

regxys :: Sample (FullNormal 2)
regxys = [S.fromTuple (x, y) | x <- range regmn regmx regres, y <- range regmn regmx regres]

chixyht :: Double
rhoxyht :: Natural # FullNormal 2
(chixyht, rhoxyht) = conjugationParameterRegression regxys lkl

regptns, regptnsht, regptndffs :: [Double]
regptns = potential <$> lkl >$>* regxys
regptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> pltxys
regptndffs = zipWith (-) regptns regptnsht

--- Probabilistic Population Code

--- Gaussian
-- prrtru, prr0 :: Source # FullNormal 2
-- prrtru = fromTuple (2, 2, 1, -1, 2)
-- prr0 = fromTuple (0, 0, 1, 0, 1)

-- ppctru, ppc0 :: Natural # ProbabilisticPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) 3)
-- ppctru = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural prrtru
-- ppc0 = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural prr0

--- GMM
type K = 2
wghtstru :: Source # Categorical K
-- wghtstru = singleton 0.3
-- wghtstru = fromTuple (0.1, 0.2, 0.3)
wghtstru = fromTuple (0.2, 0.3)

nrmstru :: S.Vector (K + 1) (Source # FullNormal 2)
nrmstru =
    S.fromTuple
        ( fromTuple (2, 2, 1, 0, 2)
        , fromTuple (-2, 2, 2, 0, 1)
        , fromTuple (-2, -2, 1, 0, 2)
        )

-- nrmstru = S.fromTuple
--     ( fromTuple (2, 2, 1, 0, 2)
--     , fromTuple (1, -1, 2, 0, 1)
--     , fromTuple (2, 1, 2, 0, 1)
--     , fromTuple (1, 2, 1, 0, 2)
--     )

prrtru :: Natural # Mixture (FullNormal 2) K
prrtru = toNatural $ joinSourceMixture nrmstru wghtstru

ppctru :: Natural # ProbabilisticPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) K)
ppctru = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural prrtru

initializePPC ::
    (KnownPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) K)) =>
    Random (Natural # ProbabilisticPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) K))
initializePPC = do
    wghts <- uniformInitialize (-1, 1)
    cmpnts' <- S.replicateM $ uniformInitialize (-0.01, 0.01)
    let cmpnts = S.zipWith (+) cmpnts' (S.replicate standardNormal)
        prr = joinNaturalMixture cmpnts wghts
    return $ approximateJoinConjugatedHarmonium rhoxyht lkl prr

--- Training

eps :: Double
eps = 3e-3

nstps, ndatsmps, nalgsmps, nepchs :: Int
nstps = 10
ndatsmps = 1000
nalgsmps = 10
nepchs = 20

gp :: GradientPursuit
gp = defaultAdamPursuit

stochasticEMStep ::
    (KnownPopulationCode n (FullNormal 2) (Mixture (FullNormal 2) K), LegendreExponentialFamily (Mixture (FullNormal 2) K)) =>
    [S.Vector n Int] ->
    (Int, Natural # ProbabilisticPopulationCode n (FullNormal 2) (Mixture (FullNormal 2) K)) ->
    IO (Int, Natural # ProbabilisticPopulationCode n (FullNormal 2) (Mixture (FullNormal 2) K))
stochasticEMStep nns (k, ppc) = do
    let mxmdl = snd $ approximateSplitConjugatedHarmonium rhoxyht ppc
    let (cmpnts, wghts) = splitSourceMixture $ toSource mxmdl
    ppc' <- realize . iterateChain nstps $ ppcExpectationMaximization nns eps nalgsmps gp rhoxyht ppc
    putStrLn $
        concat
            [ "\nIteration: "
            , show k
            , "\nWeights:"
            , show wghts
            , "\nComponents:"
            , show cmpnts
            , " Log-Likelihood: "
            , show $ ppcLogLikelihood (chixyht, rhoxyht) nns ppc'
            , "\n"
            ]
    return (k + 1, ppc')

--- Plotting

pltsmps :: Int
pltsmps = 100

pltmn, pltmx :: Double
pltmn = mumn - 2
pltmx = mumx + 2

pltxys :: Sample (FullNormal 2)
pltxys = [S.fromTuple (x, y) | x <- range pltmn pltmx pltsmps, y <- range pltmn pltmx pltsmps]

pltptns, pltptnsht, pltptndffs :: [Double]
pltptns = potential <$> lkl >$>* pltxys
pltptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> pltxys
pltptndffs = zipWith (-) pltptns pltptnsht

--- Main ---

main :: IO ()
main = do
    print ("Potential regression RMSE:" :: String)
    print . sqrt . average $ square <$> regptndffs
    print ("Conjugation parameters:" :: String)
    print (chixyht, rhoxyht)

    --- Initialization
    ppc0 :: Natural # ProbabilisticPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) K) <-
        realize initializePPC
    xys <- realize $ samplePPC ndatsmps rhoxyht ppctru
    let prr0 = snd $ approximateSplitConjugatedHarmonium rhoxyht ppc0
        xs = fst <$> xys

    --- Training
    putStrLn $ "Initial Log-Likelihood: " ++ show (ppcLogLikelihood (chixyht, rhoxyht) xs ppc0)
    kppcs <- iterateM nepchs (stochasticEMStep xs) (1, ppc0)

    let ppcs = snd <$> kppcs
        prrs = snd . approximateSplitConjugatedHarmonium rhoxyht <$> ppcs

    --- Densities
    let trudns = observableDensities prrtru pltxys
        dns0 = observableDensities prr0 pltxys
        lrndns = observableDensities (last prrs) pltxys

    let jsonData =
            toJSON
                [ "xys" .= pltxys
                , "preferred-stimuli" .= S.map coordinates mus
                , "sum-of-tuning-curves" .= pltptns
                , "estimated-sum-of-tuning-curves" .= pltptnsht
                , "estimation-difference" .= pltptndffs
                , "regression-bounds" .= (regmn, regmx)
                , "true-density" .= trudns
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                ]

    rsltfl <- resultsFilePath "population-code-2d-gaussian.json"
    exportJSON rsltfl jsonData
    putStrLn "Simulation Complete"

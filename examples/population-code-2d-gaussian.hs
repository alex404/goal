{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MonoLocalBinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Proxy (Proxy (..))

--- Globals ---

--- Types

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
regptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> regxys
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
-- type K = 3
-- wghtstru :: Source # Categorical K
-- wghtstru = fromTuple (0.1, 0.2, 0.3)
--
-- nrmstru :: S.Vector (K + 1) (Source # FullNormal 2)
-- nrmstru =
--     S.fromTuple
--         ( fromTuple (2, 2, 1, 0, 1)
--         , fromTuple (2, -2, 1, 0, 1)
--         , fromTuple (-2, 2, 1, 0, 1)
--         , fromTuple (-2, -2, 1, 0, 1)
--         )
--
-- ltnttru :: Natural # Mixture (FullNormal 2) K
-- ltnttru = toNatural $ joinSourceMixture nrmstru wghtstru
--
-- ppctru :: Natural # ProbabilisticPopulationCode NN (FullNormal 2) (Mixture (FullNormal 2) K)
-- ppctru = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural ltnttru
--
-- initializeLatent :: Random (Natural # Mixture (FullNormal 2) K)
-- initializeLatent = do
--     -- wghts :: Natural # Categorical K <- uniformInitialize (-1, 1)
--     let wghts :: Natural # Categorical K
--         wghts = 0
--         cmpnt0 :: Source # StandardNormal 2
--         cmpnt0 = 0
--         cvr0 :: Source # CovarianceMatrix L.PositiveDefinite 2
--         cvr0 = snd $ split standardNormal
--     cmpnts' <- S.replicateM $ uniformInitialize (-1, 1)
--     let cmpnts = S.map (flip join cvr0 . (+ cmpnt0)) cmpnts'
--     return . toNatural . joinSourceMixture cmpnts $ toSource wghts

--- Gaussian Boltzmann Harmonium
bltztru :: Natural # Boltzmann BN
bltztru = join (fromTuple (1, 1, 1, 1, 0, 1, 1, 1, 1)) $ -10

mubnd :: Double
mubnd = 4

shftstru :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli)
shftstru =
    fromRows $
        S.fromTuple
            ( fromTuple (-mubnd, -mubnd, -mubnd, 0, 0, 0, mubnd, mubnd, mubnd)
            , fromTuple (-mubnd, 0, mubnd, -mubnd, 0, mubnd, -mubnd, 0, mubnd)
            )

smvntru :: Source # FullNormal OS
smvntru = standardNormal / 2

mvntru :: Natural # FullNormal OS
mvntru = toNatural smvntru

lmdltru :: Natural # BoltzmannLinearModel L.PositiveDefinite OS BN
lmdltru = join mvntru shftstru

ltnttru :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN
ltnttru = joinConjugatedHarmonium lmdltru bltztru

ppctru :: Natural # ProbabilisticPopulationCode NN (FullNormal 2) (GaussianBoltzmannHarmonium L.PositiveDefinite OS BN)
ppctru = approximateJoinConjugatedHarmonium rhoxyht lkl ltnttru

unibnd :: Double
unibnd = 4

initializeLatent ::
    Random (Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN)
initializeLatent = do
    bltz :: Natural # Boltzmann BN <- uniformInitialize (-unibnd, unibnd)
    shfts :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli) <-
        uniformInitialize (-unibnd / 2, unibnd / 2)
    let mvn = standardNormal
        lmdl = join mvn shfts
    return $ joinConjugatedHarmonium lmdl bltz

--- Training

eps :: Double
eps = 3e-3

nstps, ndatsmps, nepchs :: Int
nstps = 2000
ndatsmps = 1000
nepchs = 100

gp :: GradientPursuit
gp = defaultAdamPursuit

loggingEMStep ::
    (LinearSubspace y (FullNormal 2), LegendreExponentialFamily y) =>
    [Natural # FullNormal 2] ->
    Sample (Replicated NN Poisson) ->
    (Int, Natural # y) ->
    IO (Int, Natural # y)
loggingEMStep ny0s nss (k, nltnt) = do
    let nltnt' = ppcExpectationMaximizationAscent eps gp ny0s nltnt !! nstps
        ppc' = approximateJoinConjugatedHarmonium rhoxyht lkl nltnt'
    putStrLn
        . concat
        $ [ "\nIteration: "
          , show k
          , "\nLog-Likelihood: "
          , show $ ppcLogLikelihood (chixyht, rhoxyht) nss ppc'
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
dnspltxys3 :: Sample (FullNormal 2)
dnspltxys3 = [S.fromTuple (x, y) | x <- range dnspltmn dnspltmx pltres, y <- range dnspltmn dnspltmx pltres]

--- Main ---

main :: IO ()
main = do
    print ("\nPotential regression RMSE:" :: String)
    print . sqrt . average $ square <$> regptndffs
    let trudns = observableDensities ltnttru dnspltxys
    print $ sum trudns

    --- Initialization
    ltnt0 <- realize initializeLatent
    let ppc0 = approximateJoinConjugatedHarmonium rhoxyht lkl $ toNatural ltnt0
    xys <- realize $ samplePPC ndatsmps rhoxyht ppctru
    let xs = fst <$> xys
        dns0 = observableDensities ltnt0 dnspltxys2
    print $ sum dns0

    --- Training
    putStrLn $ "Initial Log-Likelihood: " ++ show (ppcLogLikelihood (chixyht, rhoxyht) xs ppc0)

    let ny0s = ppcExpectationBiases xs lkl
    kltnts <- iterateM nepchs (loggingEMStep ny0s xs) (1, ltnt0)
    putStrLn ("\nTraining Complete" :: String)

    let ltnts = snd <$> kltnts
        lrndns = observableDensities (last ltnts) dnspltxys3
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
                , "true-density" .= trudns
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                ]

    rsltfl <- resultsFilePath "population-code-2d-gaussian.json"
    exportJSON rsltfl jsonData
    putStrLn "\nSimulation Complete\n"

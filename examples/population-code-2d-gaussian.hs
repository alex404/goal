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
type BN = 6
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
gns = Point $ S.replicate 2

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

regys :: Sample (FullNormal 2)
regys = [S.fromTuple (y1, y2) | y1 <- range regmn regmx regres, y2 <- range regmn regmx regres]

chixyht :: Double
rhoxyht :: Natural # FullNormal 2
(chixyht, rhoxyht) = conjugationParameterRegression regys lkl

regptns, regptnsht, regptndffs :: [Double]
regptns = potential <$> lkl >$>* regys
regptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> regys
regptndffs = zipWith (-) regptns regptnsht

--- Probabilistic Population Code

unibnd :: Double
unibnd = 1

initializeLatent ::
    Random (Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN)
initializeLatent = do
    bltz :: Natural # Boltzmann BN <- uniformInitialize (-unibnd, unibnd)
    shfts :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli) <-
        uniformInitialize (-unibnd, unibnd)
    let mvn = standardNormal
        lmdl = join mvn shfts
    return $ joinConjugatedHarmonium lmdl bltz

--- Data generation

circle2d :: Double -> Double -> SamplePoint (FullNormal 2)
circle2d rds t = S.fromTuple (rds * cos t, rds * sin t)

nsmdl :: Source # FullNormal 2
nsmdl = fromTuple (0, 0, 0.2, 0, 0.2)

noisyCircle :: Double -> Int -> Random (Sample (FullNormal 2))
noisyCircle rds nsmps = do
    let intrvl = range 0 (2 * pi) nsmps
    mapM (noisyFunction nsmdl (circle2d rds)) intrvl

--- Training

eps :: Double
eps = 3e-2

nstps, nobs, nepchs :: Int
nstps = 1000
nobs = 1000
nepchs = 500

gp :: GradientPursuit
gp = defaultAdamPursuit

loggingEMStep ::
    (LinearSubspace y (FullNormal 2), LegendreExponentialFamily y) =>
    [Natural # FullNormal 2] ->
    Sample (Replicated NN Poisson) ->
    (Int, Double, Natural # y) ->
    IO (Int, Double, Natural # y)
loggingEMStep ny0s nss (k, ll, nltnt) = do
    putStrLn
        . concat
        $ [ "\nIteration: "
          , show k
          , "\nLog-Likelihood: "
          , show ll
          ]
    let nltnt1 = ppcExpectationMaximizationAscent eps gp ny0s nltnt !! nstps
        ll1 = ppcLogLikelihood (chixyht, rhoxyht) nss $ joinConjugatedHarmonium0 rhoxyht lkl nltnt1
    return (k + 1, ll1, nltnt1)

--- Plotting

pltres :: Int
pltres = 100

ptnpltmn, ptnpltmx :: Double
ptnpltmn = regmn - 2
ptnpltmx = regmx + 2

ptnpltys :: Sample (FullNormal 2)
ptnpltys = [S.fromTuple (y1, y2) | y1 <- range ptnpltmn ptnpltmx pltres, y2 <- range ptnpltmn ptnpltmx pltres]

pltptns, pltptnsht, pltptndffs :: [Double]
pltptns = potential <$> lkl >$>* ptnpltys
pltptnsht = map (+ chixyht) . dotMap rhoxyht $ sufficientStatistic <$> ptnpltys
pltptndffs = zipWith (-) pltptns pltptnsht

dnspltmn, dnspltmx :: Double
dnspltmn = regmn
dnspltmx = regmx

dnspltys :: Sample (FullNormal 2)
dnspltys = [S.fromTuple (y1, y2) | y1 <- range dnspltmn dnspltmx pltres, y2 <- range dnspltmn dnspltmx pltres]

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

    --- Initialization
    ltnt0 <- realize initializeLatent
    let ppc0 = joinConjugatedHarmonium0 rhoxyht lkl $ toNatural ltnt0
    let frc :: Double
        frc = 0.85
        nobs1 = round $ frc * fromIntegral nobs
        nobs2 = nobs - nobs1

    ys' <- realize $ noisyCircle 3 nobs1
    ys'' <- realize $ noisyCircle 0.1 nobs2
    ys <- realize . shuffleList $ ys' ++ ys''
    xs <- realize . mapM samplePoint $ lkl >$>* ys
    let dns0 = observableDensities ltnt0 dnspltxys2
    print $ sum dns0

    --- Training
    let ll0 = ppcLogLikelihood (chixyht, rhoxyht) xs ppc0
    putStrLn $ "Initial Log-Likelihood: " ++ show ll0

    let ny0s = ppcExpectationBiases xs lkl
    kltnts <- iterateM nepchs (loggingEMStep ny0s xs) (0, ll0, ltnt0)
    putStrLn ("\nTraining Complete" :: String)

    let (_, lls, ltnts) = unzip3 kltnts
        lrndns = observableDensities (last ltnts) dnspltxys3
    print $ sum lrndns

    let jsonData =
            toJSON
                [ "tuning-curve-ys" .= ptnpltys
                , "y-samples" .= ys
                , "x-samples" .= take 10 xs
                , "preferred-stimuli" .= S.map coordinates mus
                , "sum-of-tuning-curves" .= pltptns
                , "estimated-sum-of-tuning-curves" .= pltptnsht
                , "estimation-difference" .= pltptndffs
                , "regression-bounds" .= (regmn, regmx)
                , "density-ys" .= dnspltys
                , "initial-density" .= dns0
                , "learned-density" .= lrndns
                , "log-likelihoods" .= lls
                ]

    rsltfl <- resultsFilePath "population-code-2d-gaussian.json"
    exportJSON rsltfl jsonData
    putStrLn "\nSimulation Complete\n"

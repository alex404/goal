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
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Proxy (Proxy (..))

--- Globals ---

--- Types

type HierarchicalMixtureOfGaussians t n m k =
    AffineHarmonium
        L.Full
        (StandardNormal n)
        (StandardNormal m)
        (MultivariateNormal t n)
        (Mixture (FullNormal m) k)

type FullHMoG n m k = HierarchicalMixtureOfGaussians L.PositiveDefinite n m k
type DiagonalHMoG n m k = HierarchicalMixtureOfGaussians L.Diagonal n m k
type IsotropicHMoG n m k = HierarchicalMixtureOfGaussians L.Scale n m k

--- Instances

instance (KnownNat m, KnownNat k) => LinearSubspace (Mixture (FullNormal m) k) (StandardNormal m) where
    {-# INLINE (>+>) #-}
    (>+>) hrm nx0 =
        let (nx, nxz, nz) = splitHarmonium hrm
         in joinHarmonium (nx >+> nx0) nxz nz
    linearProjection hrm =
        let (nx, _, _) = splitHarmonium hrm
         in linearProjection nx

instance
    (KnownNat k, KnownNat m, KnownCovariance f n) =>
    ConjugatedLikelihood
        L.Full
        (StandardNormal n)
        (StandardNormal m)
        (MultivariateNormal f n)
        (Mixture (FullNormal m) k)
    where
    conjugationParameters lm =
        let (rho0, rprms) = conjugationParameters lm
         in (rho0, join (join rprms 0) 0)

--- HMoG

wdth, hght, sep :: Double
wdth = 1
hght = 4
sep = 6

trusz :: Source # Categorical 1
trusz = singleton 0.5

trusy0, trusy1 :: Source # FullNormal 1
trusy0 = fromTuple (-(sep / 2), 1)
trusy1 = fromTuple (sep / 2, 1)

trusyz :: Source # Mixture (FullNormal 1) 1
trusyz = joinSourceMixture (S.fromTuple (trusy0, trusy1)) trusz

trusx :: Source # DiagonalNormal 2
trusx = fromTuple (0, 0, wdth, hght)

trusfa :: Source # FactorAnalysis 2 1
trusfa = join trusx $ fromTuple (1, 0)

trunhmog :: Natural # DiagonalHMoG 2 1 1
trunhmog = joinConjugatedHarmonium (toNatural trusfa) $ toNatural trusyz

--- Training

eps :: Double
eps = 3e-3

nstps, smpn, nepchs :: Int
nstps = 2000
smpn = 1000
nepchs = 1000

gp :: GradientPursuit
gp = defaultAdamPursuit

loggingEMStep ::
    (KnownNat k, KnownNat m, KnownCovariance t n) =>
    [S.Vector n Double] ->
    (Int, Natural # HierarchicalMixtureOfGaussians t n m k) ->
    IO (Int, Natural # HierarchicalMixtureOfGaussians t n m k)
loggingEMStep xs (k, gbhrm) = do
    let gbhrms = expectationMaximizationAscent eps gp xs gbhrm
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show $ logLikelihood xs gbhrm
            ]
    return (k + 1, gbhrms !! nstps)

--- Plotting

pltres :: Int
pltres = 100

pltmnx1, pltmnx2, pltmxx1, pltmxx2 :: Double
pltmnx1 = -6
pltmnx2 = -6
pltmxx1 = 6
pltmxx2 = 6

pltx1s, pltx2s :: [Double]
pltx1s = range pltmnx1 pltmxx1 pltres
pltx2s = range pltmnx2 pltmxx2 pltres

--- Main ---

main :: IO ()
main = do
    xyzs <- realize $ sample smpn trunhmog
    let (xs, _) = unzip xyzs

    let trudnss = observableDensities trunhmog [S.fromTuple (x1, x2) | x2 <- pltx2s, x1 <- pltx1s]

    --- Export
    let json =
            toJSON
                [ "plot-range-x1" .= pltx1s
                , "plot-range-x2" .= pltx2s
                , "true-observable-density" .= trudnss
                , "observations" .= xs
                ]

    --- Process data
    flnm <- resultsFilePath "hmog.json"

    exportJSON flnm json

-- goalExport ldpth "data" $ S.toList <$> xs
-- goalExport ldpth "labels" $ Only . snd <$> yzs
-- writeHMoG ldpth nhmog
--
-- putStrLn "True Log-Likelihood:"
-- print $ logLikelihood xs nhmog

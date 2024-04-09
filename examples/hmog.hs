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

instance
    (KnownNat k, KnownNat m, KnownCovariance t n) =>
    Transition Natural Mean (HierarchicalMixtureOfGaussians t n m k)
    where
    transition hmog =
        let (nx, nxy, nhrmyz) = splitHarmonium hmog
            (ny, nyz, nz) = splitHarmonium nhrmyz
            hmog2 = joinHarmonium (joinHarmonium nx nxy ny) nyz nz
            mhmog2 = toMean hmog2
            (mhrmxy, myz, mz) = splitHarmonium mhmog2
            (mx, mxy, my) = splitHarmonium mhrmxy
         in joinHarmonium mx mxy $ joinHarmonium my myz mz

instance
    (KnownNat k, KnownNat m, KnownCovariance t n) =>
    Transition Mean Natural (HierarchicalMixtureOfGaussians t n m k)
    where
    transition mlgh =
        let (mx, msgmaxz, mz) = splitHarmonium mlgh
            mz' :: Mean # FullNormal m
            mz' = linearProjection mz
            lmdl = joinHarmonium mx msgmaxz mz'
            nlkl = fst . split $ toNatural lmdl
         in joinConjugatedHarmonium nlkl $ toNatural mz

--- HMoG

wdth, hght, sep :: Double
wdth = 0.1
hght = 8
sep = 5

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

truhmog :: Natural # DiagonalHMoG 2 1 1
truhmog = joinConjugatedHarmonium (toNatural trusfa) $ toNatural trusyz

--- Initialization

-- initializeHMoG :: Random (Natural # DiagonalHMoG 2 1 1)
-- initializeHMoG = do
--     nz <- uniformInitialize (-1, 1)
--     let sy1, sy2 :: Source # FullNormal 1
--         sy1 = fromTuple (-1, 10)
--         sy2 = fromTuple (1, 10)
--         ny1 = toNatural sy1
--         ny2 = toNatural sy2
--     let nyz = joinNaturalMixture (S.fromTuple (ny1, ny2)) nz
--         sx :: Source # DiagonalNormal 2
--         sx = fromTuple (0, 0, 10, 10)
--         nx = toNatural sx
--     nxy <- uniformInitialize (-1, 1)
--     let nfa = join nx nxy
--     return $ joinConjugatedHarmonium nfa nyz
--
initializeHMoG :: Random (Natural # DiagonalHMoG 2 1 1)
initializeHMoG = do
    nz <- uniformInitialize (-1, 1)
    ny1' <- uniformInitialize (-0.1, 0.1)
    ny2' <- uniformInitialize (-0.1, 0.1)
    let sy1, sy2 :: Source # FullNormal 1
        sy1 = fromTuple (-1.5, 1.5)
        sy2 = fromTuple (1.5, 1.5)
    let ny1 = toNatural sy1 + ny1'
        ny2 = toNatural sy2 + ny2'

    let nyz = joinNaturalMixture (S.fromTuple (ny1, ny2)) nz

    -- let mx = averageSufficientStatistic xs
    let sx :: Source # DiagonalNormal 2
        sx = fromTuple (0, 0, 2, 2)
    nx' <- uniformInitialize (-0.1, 0.1)
    let nx = toNatural sx + nx'
    nxy <- uniformInitialize (-1, 1)
    let nfa = join nx nxy
    return $ joinConjugatedHarmonium nfa nyz

--- Training

nobs, nepchs :: Int
nobs = 400
nepchs = 200

loggingEMStep ::
    (KnownNat k, KnownNat m, KnownCovariance t n) =>
    [S.Vector n Double] ->
    (Int, Double, Natural # HierarchicalMixtureOfGaussians t n m k) ->
    IO (Int, Double, Natural # HierarchicalMixtureOfGaussians t n m k)
loggingEMStep xs (k, ll, hmog) = do
    let hmog1 = expectationMaximization xs hmog
        ll1 = logLikelihood xs hmog1
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show ll
            ]
    return (k + 1, ll1, hmog1)

--- Plotting

pltres :: Int
pltres = 100

pltmnx1, pltmnx2, pltmxx1, pltmxx2 :: Double
pltmnx1 = -10
pltmnx2 = -10
pltmxx1 = 10
pltmxx2 = 10

pltx1s, pltx2s :: [Double]
pltx1s = range pltmnx1 pltmxx1 pltres
pltx2s = range pltmnx2 pltmxx2 pltres

pltmny, pltmxy :: Double
pltmny = -6
pltmxy = 6

pltys :: [Double]
pltys = range pltmny pltmxy pltres

--- Main ---

main :: IO ()
main = do
    xyzs <- realize $ sample nobs truhmog
    let (xs, _) = unzip xyzs

    hmog0 <- realize $ initializeHMoG
    let ll0 = logLikelihood xs hmog0

    -- khmogs1 <- iterateM nepchs (loggingEMAscentStep xs) (0, hmog0)

    kllhmogs <- iterateM nepchs (loggingEMStep xs) (0, ll0, hmog0)
    let (_, lls, hmogs) = unzip3 kllhmogs

    let hmog1 = last hmogs

    let truxdns = observableDensities truhmog [S.fromTuple (x1, x2) | x2 <- pltx2s, x1 <- pltx1s]
        hmog0xdns = observableDensities hmog0 [S.fromTuple (x1, x2) | x2 <- pltx2s, x1 <- pltx1s]
        hmog1xdns = observableDensities hmog1 [S.fromTuple (x1, x2) | x2 <- pltx2s, x1 <- pltx1s]

    let trumxmdl = snd $ splitConjugatedHarmonium truhmog
        mxmdl0 = snd $ splitConjugatedHarmonium hmog0
        mxmdl1 = snd $ splitConjugatedHarmonium hmog1

    let truydns = observableDensities trumxmdl $ S.singleton <$> pltys
        hmog0ydns = observableDensities mxmdl0 $ S.singleton <$> pltys
        hmog1ydns = observableDensities mxmdl1 $ S.singleton <$> pltys

    --- Export
    let json =
            toJSON
                [ "plot-range-x1" .= pltx1s
                , "plot-range-x2" .= pltx2s
                , "observations" .= xs
                , "log-likelihoods" .= lls
                , "observable-densities" .= [truxdns, hmog0xdns, hmog1xdns]
                , "plot-range-y" .= pltys
                , "mixture-densities" .= [truydns, hmog0ydns, hmog1ydns]
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

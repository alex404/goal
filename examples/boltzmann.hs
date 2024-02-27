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

-- Misc

import Control.Monad (replicateM)
import Data.Tuple.Only (Only (..))

--- Globals ---

--- Distribution

type Neurons = 4
type ObservationSpace = 1

inttru :: Natural # InteractionMatrix Neurons
inttru = -fromTuple (0.5, -0.5, 0.2, -0.5, 1, 0.5)

bsstru :: Natural # Replicated Neurons Bernoulli
bsstru = -fromTuple (0.1, -0.25, 0.25, -0.1)

bltztru :: Natural # Boltzmann Neurons
bltztru = join bsstru inttru

bltz0 :: Natural # Boltzmann Neurons
bltz0 = 0

--- Linear Model

smvn :: Source # FullNormal ObservationSpace
smvn = fromTuple (1, 2)
mvn = toNatural smvn

intltnt :: Natural # Tensor (StandardNormal ObservationSpace) (Replicated Neurons Bernoulli)
intltnt =
    fromRows . S.fromTuple $
        Only $
            fromTuple (1.5, 0.5, -1, -0.5)

lmdl :: Natural # BoltzmannLinearModel L.PositiveDefinite ObservationSpace Neurons
lmdl = join mvn intltnt

rho0 :: Double
rprms :: Natural # Boltzmann Neurons
(rho0, rprms) = linearBoltzmannHarmoniumConjugationParameters lmdl

hrmtru :: Natural # GaussianBoltzmannHarmonium L.PositiveDefinite ObservationSpace Neurons
hrmtru = joinConjugatedHarmonium lmdl bltztru

--- Simulation

smpsz :: Int
smpsz = 10000

ncycs :: Int
ncycs = 1000

eps :: Double
eps = 0.05

nepchs :: Int
nepchs = 500

--- Plotting

--- Functions

fitBoltzmann :: Sample (Boltzmann Neurons) -> [Natural # Boltzmann Neurons]
fitBoltzmann blss =
    vanillaGradientSequence (logLikelihoodDifferential blss) eps defaultAdamPursuit bltz0

bruteForceDensity ::
    S.Vector ObservationSpace Double ->
    Double
bruteForceDensity x =
    sum [density bltztru bls * density (lmdl >.>* bls) x | bls <- pointSampleSpace bltztru]

bruteForceDensity2 ::
    S.Vector ObservationSpace Double ->
    Double
bruteForceDensity2 x =
    sum [density hrmtru (x, bls) | bls <- pointSampleSpace bltztru]

-- bruteForceDensity3 ::
--     S.Vector ObservationSpace Double ->
--     Double
-- bruteForceDensity3 x = fst $ integrate 1e-4 (\y -> density hrm2 (x, S.singleton y)) (-5) 5
--
-- bruteForceDensity4 ::
--     S.Vector ObservationSpace Double ->
--     Double
-- bruteForceDensity4 x =
--     fst $ integrate 1e-4 (\y -> density prr2 (S.singleton y) * density (lmdl2 >.>* S.singleton y) x) (-5) 5
--
--- Integration testing

-- prr2 :: Natural # FullNormal 1
-- prr2 = standardNormal
--
-- int2 :: Natural # Tensor (StandardNormal 2) (StandardNormal 1)
-- int2 =
--     fromRows . S.fromTuple $
--         ( fromTuple $ Only 1.5
--         , fromTuple $ Only (-1.5)
--         )
--
-- lmdl2 :: Natural # FullLinearModel ObservationSpace 1
-- lmdl2 = join mvn int2
--
-- hrm2 :: Natural # FullGaussianHarmonium 2 1
-- hrm2 = joinConjugatedHarmonium lmdl2 prr2
--
--- Main ---

main :: IO ()
main = do
    let xs :: Sample (FullNormal ObservationSpace)
        -- xs = S.fromTuple <$> [(1, 0), (0, 1), (-1, 0), (0, -1)]
        xs = S.singleton <$> [-2 .. 2]

    -- cnjdns = exp <$> logConjugatedDensities rprms hrmtru xs
    -- brtdns = bruteForceDensity <$> xs
    -- brtdns2 = bruteForceDensity2 <$> xs

    print . sum . densities bltztru $ pointSampleSpace bltztru
    print . sum $ [(* density bltztru bls) . fst $ integrate 1e-4 (density (lmdl >.>* bls) . S.singleton) (-10) 10 | bls <- pointSampleSpace bltztru]
    print . fst $ integrate 1e-4 (\x -> exp . head $ logConjugatedDensities (rho0, rprms) hrmtru [S.singleton x]) (-10) 10
    print . sum $ [fst $ integrate 1e-4 (\x -> density hrmtru (S.singleton x, bls)) (-10) 10 | bls <- pointSampleSpace bltztru]
    print . fst $ integrate 1e-4 (\x -> sum [density hrmtru (S.singleton x, bls) | bls <- pointSampleSpace bltztru]) (-10) 10

-- print rprms
-- mapM_ print [(potential (lmdl >.>* bls), rprms <.> sufficientStatistic bls + rho0, bls) | bls <- pointSampleSpace bltztru]

-- smp <- replicateM smpsz . realize $ gibbsBoltzmann ncycs bltztru
--
-- let bltzs = take nepchs $ fitBoltzmann smp
--     dlls = logLikelihood smp <$> bltzs
--
-- print dlls
--
-- mapM_ print . zip (listCoordinates bltztru) . listCoordinates $ last bltzs

-- print cnjdns
-- print brtdns
-- print brtdns2

-- print $ exp <$> logObservableDensities hrm2 xs
-- print $ bruteForceDensity3 <$> xs
-- print $ bruteForceDensity4 <$> xs

--     drch = last drchs
--     dxrng = range dmn dmx pltsmps
--     dyrng = range dmn dmx pltsmps
--
-- let ddnss = do
--         y <- dyrng
--         x <- dxrng
--         return (density2d dtru (x, y), density2d drch (x, y))
--
-- let (dtdns, dldns) = unzip ddnss
--
-- let djson =
--         toJSON
--             [ "log-likelihood" .= dlls
--             , "xrange" .= dxrng
--             , "yrange" .= dyrng
--             , "true-density" .= dtdns
--             , "learned-density" .= dldns
--             , "samples" .= map S.toList dsmps
--             ]
--
--     --- Process data
--     mvnfl <- resultsFilePath "multivariate-normal.json"
--     drchfl <- resultsFilePath "multivariate-dirichlet.json"
--
--     exportJSON mvnfl njson
--     exportJSON drchfl djson

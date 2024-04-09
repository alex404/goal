{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S

--- Misc

import Data.List qualified as L
import Data.Proxy (Proxy (..))

--- Globals ---

-- Types --

type N = 8

type Neurons n = Replicated n Poisson

-- Affine Population

kp :: Double
kp = 1

mus :: S.Vector N Double
mus = S.generate (\k -> 2 * pi * fromIntegral k / fromIntegral (natVal (Proxy @N)))

tcs :: S.Vector N (Source # VonMises)
tcs = S.map (fromTuple . (,kp)) mus

zs :: [Double]
zs = range 0 (2 * pi) 1000

gnf :: Double -> Double
gnf z = 1 + 0.25 * sin z

gns :: Source # Neurons N
gns = Point $ S.map gnf mus

lkl :: Natural # Neurons N <* VonMises
lkl = joinPopulationCode (toNatural gns) (S.map toNatural tcs)

ptns :: [Double]
ptns = potential <$> lkl >$>* zs

chiht :: Double
rhoht :: Natural # VonMises
(chiht, rhoht) = conjugationParameterRegression zs lkl

ptnsht :: [Double]
ptnsht = conjugationCurve chiht rhoht zs

tcsmps :: [[Double]]
tcsmps = L.transpose $ listCoordinates . toMean <$> lkl >$>* zs

-- Prior

prr :: Source # VonMises
prr = fromTuple (2 * pi / 2, 1)

nprr :: Natural # VonMises
nprr = toNatural prr

prrdns :: [Double]
prrdns = densities prr zs

-- Posteriors

nobs :: Int
nobs = 2

observe :: Random [S.Vector N Int]
observe = do
    let zsmps = tail . init $ range 0 (2 * pi) (nobs + 2)
    mapM (samplePoint . (lkl >.>*)) zsmps

pstdnss :: [S.Vector N Int] -> [[Double]]
pstdnss obss =
    [densities pst zs | pst <- conjugatedBayesRule0 rhoht lkl nprr <$> obss]

-- Observable Covariance

nsmpcrl :: Int
nsmpcrl = 1000

estimateCorrelationMatrixRows :: Random [S.Vector N Double]
estimateCorrelationMatrixRows = do
    let zsmps = tail $ range 0 (2 * pi) nsmpcrl
    xsmps0 <- mapM (samplePoint . (lkl >.>*)) zsmps
    let xsmps = G.convert . G.map realToFrac <$> xsmps0
        dnss = [dns * 2 * pi / fromIntegral nsmpcrl | dns <- densities prr zs]
        mu = Point . weightedAverage $ zip dnss xsmps
        mcvr = Point . S.lowerTriangular . S.weightedAverageOuterProduct $ zip3 dnss xsmps xsmps
        mnrm :: Mean # FullNormal N
        mnrm = join mu mcvr
    return . S.toList . S.map coordinates . toRows . multivariateNormalCorrelations $ toSource mnrm

--- Main ---

main :: IO ()
main = do
    crlrws <- realize estimateCorrelationMatrixRows
    obss <- realize observe
    let jsonData =
            toJSON
                [ "plot-zs" .= zs
                , "preferred-stimuli" .= mus
                , "prior-density" .= prrdns
                , "sum-of-tuning-curves" .= ptns
                , "regression" .= ptnsht
                , "tuning-curves" .= tcsmps
                , "correlation-matrix" .= crlrws
                , "observations" .= obss
                , "posterior-densities" .= pstdnss obss
                ]

    rsltfl <- resultsFilePath "population-code-von-mises.json"
    exportJSON rsltfl jsonData

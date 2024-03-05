{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Misc

import Data.List qualified as L
import Data.Proxy (Proxy (..))

--- Globals ---

-- Types --

type N = 10

type Neurons n = Replicated n Poisson

-- Variables --

kp :: Double
kp = 1

mus :: S.Vector N Double
mus = S.generate (\k -> 2 * pi * fromIntegral k / fromIntegral (natVal (Proxy @N)))

tcs :: S.Vector N (Source # VonMises)
tcs = S.map (fromTuple . (,kp)) mus

zs :: [Double]
zs = range 0 (2 * pi) 1000

-- Linear Population

gns1 :: Source # Neurons N
gns1 = Point $ S.replicate 2

lkl1 :: Natural # Neurons 10 <* VonMises
lkl1 = joinPopulationCode (toNatural gns1) (S.map toNatural tcs)

ptns1 :: [Double]
ptns1 = potential <$> lkl1 >$>* zs

rh01ht :: Double
rprms1ht :: Natural # VonMises
(rh01ht, rprms1ht) = conjugationParameterRegression zs lkl1

ptnsht1 :: [Double]
ptnsht1 = map (+ rh01ht) . dotMap rprms1ht $ sufficientStatistic <$> zs

tcsmps1 :: [[Double]]
tcsmps1 = L.transpose $ listCoordinates . toMean <$> lkl1 >$>* zs

-- Affine Population

gnf :: Double -> Double
gnf z = 2 + 0.5 * sin z

gns2 :: Source # Neurons N
gns2 = Point $ S.map gnf mus

lkl2 :: Natural # Neurons 10 <* VonMises
lkl2 = joinPopulationCode (toNatural gns2) (S.map toNatural tcs)

ptns2 :: [Double]
ptns2 = potential <$> lkl2 >$>* zs

rh02ht :: Double
rprms2ht :: Natural # VonMises
(rh02ht, rprms2ht) = conjugationParameterRegression zs lkl2

ptnsht2 :: [Double]
ptnsht2 = map (+ rh02ht) . dotMap rprms2ht $ sufficientStatistic <$> zs

tcsmps2 :: [[Double]]
tcsmps2 = L.transpose $ listCoordinates . toMean <$> lkl2 >$>* zs

--- Main ---

main :: IO ()
main = do
    let jsonData =
            toJSON
                [ "zs" .= zs
                , "linear-sum-of-tuning-curves" .= ptns1
                , "affine-sum-of-tuning-curves" .= ptns2
                , "linear-regression" .= ptnsht1
                , "affine-regression" .= ptnsht2
                , "linear-tuning-curves" .= tcsmps1
                , "affine-tuning-curves" .= tcsmps2
                ]

    rsltfl <- resultsFilePath "population-code-von-mises.json"
    exportJSON rsltfl jsonData

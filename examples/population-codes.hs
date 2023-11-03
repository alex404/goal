{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TupleSections #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Misc

import Data.List qualified as L
import Data.Proxy (Proxy (..))

--- Globals ---

-- Types --

type N = 10

type Neurons n = Replicated n Poisson

-- Functions --

joinVonMisesIndependent ::
    (KnownNat n) =>
    -- | Gains
    Natural # Neurons n ->
    -- | Von Mises Curves
    S.Vector n (Natural # VonMises) ->
    -- | Population Likelihood
    Natural # Neurons n <* VonMises
joinVonMisesIndependent nz0 nps =
    let mtx = fromRows nps
        nz = nz0 - Point (S.map potential nps)
     in join nz mtx

-- Variables --

kp :: Double
kp = 1

mus :: S.Vector N Double
mus = S.generate (\k -> 2 * pi * fromIntegral k / fromIntegral (natVal (Proxy @N)))

tcs :: S.Vector N (Source # VonMises)
tcs = S.map (fromTuple . (,kp)) mus

xs :: [Double]
xs = range 0 (2 * pi) 1000

sx :: Double -> S.Vector 3 Double
sx x = S.fromTuple (1, cos x, sin x)

-- Linear Population

gns1 :: Source # Neurons N
gns1 = Point $ S.replicate 2

fzx1 :: Natural # Neurons 10 <* VonMises
fzx1 = joinVonMisesIndependent (toNatural gns1) (S.map toNatural tcs)

ys1 :: [Double]
ys1 = potential <$> fzx1 >$>* xs

bts1 :: S.Vector 3 Double
bts1 = S.linearLeastSquares (sx <$> xs) ys1

yhts1 :: [Double]
yhts1 = S.dotMap bts1 $ sx <$> xs

fzxss1 :: [[Double]]
fzxss1 = L.transpose $ listCoordinates . toMean <$> fzx1 >$>* xs

-- Affine Population

gnf :: Double -> Double
gnf x = 2 + 0.5 * sin x

gns2 :: Source # Neurons N
gns2 = Point $ S.map gnf mus

fzx2 :: Natural # Neurons 10 <* VonMises
fzx2 = joinVonMisesIndependent (toNatural gns2) (S.map toNatural tcs)

ys2 :: [Double]
ys2 = potential <$> fzx2 >$>* xs

bts2 :: S.Vector 3 Double
bts2 = S.linearLeastSquares (sx <$> xs) ys2

yhts2 :: [Double]
yhts2 = S.dotMap bts2 $ sx <$> xs

fzxss2 :: [[Double]]
fzxss2 = L.transpose $ listCoordinates . toMean <$> fzx2 >$>* xs

--- Main ---

main :: IO ()
main = do
    let jsonData =
            toJSON
                [ "xs" .= xs
                , "linear-sum-of-tuning-curves" .= ys1
                , "affine-sum-of-tuning-curves" .= ys2
                , "linear-regression" .= yhts1
                , "affine-regression" .= yhts2
                , "linear-tuning-curves" .= fzxss1
                , "affine-tuning-curves" .= fzxss2
                ]

    rsltfl <- resultsFilePath "population-codes.json"
    exportJSON rsltfl jsonData

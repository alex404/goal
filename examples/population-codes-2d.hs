{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeOperators #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Misc

import Data.List qualified as L

--- Globals ---

-- Types --

type N = 25

type Neurons n = Replicated n Poisson

-- Functions --

joinPopulationCode ::
    (KnownNat n, LegendreExponentialFamily x) =>
    -- | Gains
    Natural # Neurons n ->
    -- | Von Mises Curves
    S.Vector n (Natural # x) ->
    -- | Population Likelihood
    Natural # Neurons n <* x
joinPopulationCode nz0 nps =
    let mtx = fromRows nps
        nz = nz0 - Point (S.map potential nps)
     in join nz mtx

-- Variables --

mu0s :: S.Vector 5 Double
mu0s = S.fromTuple (-2, -1, 0, 1, 2)

mus :: S.Vector N (Source # StandardNormal 2)
mus = S.concatMap (\x -> S.map (\y -> fromTuple (x, y)) mu0s) mu0s

tcs :: S.Vector N (Source # DiagonalNormal 2)
tcs = S.map (`join` 0.1) mus

mn, mx :: Double
mn = -3
mx = 3

xys :: Sample (DiagonalNormal 2)
xys = [S.fromTuple (x, y) | x <- range mn mx 200, y <- range mn mx 200]

-- Linear Population

gns1 :: Source # Neurons N
gns1 = Point $ S.replicate 2

fzx1 :: Natural # Neurons N <* DiagonalNormal 2
fzx1 = joinPopulationCode (toNatural gns1) (S.map toNatural tcs)

ys1 :: [Double]
ys1 = potential <$> fzx1 >$>* xys

fzxss1 :: [[Double]]
fzxss1 = L.transpose $ listCoordinates . toMean <$> fzx1 >$>* xys

--- Main ---

main :: IO ()
main = do
    let jsonData =
            toJSON
                [ "xys" .= xys
                , "sum-of-tuning-curves" .= ys1
                , "tuning-curves" .= fzxss1
                ]

    rsltfl <- resultsFilePath "population-codes-2d.json"
    exportJSON rsltfl jsonData

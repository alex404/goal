{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Globals ---

pltsmps :: Int
pltsmps = 100

--- Multivariate Normal

type TrueNormal = FullNormal 2
type FitNormal = FullNormal 2

nnsmps :: Int
nnsmps = 20

mux, muy, vrx, cvr, vry :: Double
mux = 1
muy = -1
vrx = 3
cvr = 1
vry = 2

ntru :: Source # TrueNormal
ntru = fromTuple (mux, muy, vrx, cvr, vry)

nxmn, nxmx, nymn, nymx :: Double
nxmn = -6
nxmx = 10
nymn = -6
nymx = 8

--- Dirichlet

-- Globals
alphs :: S.Vector 3 Double
alphs = S.fromTuple (3, 7, 5)

dtru :: Natural # Dirichlet 3
dtru = Point alphs

dmn, dmx :: Double
dmn = 0
dmx = 1

-- Training
eps :: Double
eps = 0.05

dnsmps :: Int
dnsmps = 20

nepchs :: Int
nepchs = 500

drch0 :: Natural # Dirichlet 3
drch0 = fromTuple (1, 1, 1)

-- Functions
fitDirichlet :: Sample (Dirichlet 3) -> [Natural # Dirichlet 3]
fitDirichlet xyzs =
    vanillaGradientSequence (logLikelihoodDifferential xyzs) eps defaultAdamPursuit drch0

density2d :: Natural # Dirichlet 3 -> (Double, Double) -> Double
density2d drch (x, y) =
    let z = 1 - x - y
     in if x + y < 0.995
            then density drch $ S.fromTuple (x, y, z)
            else 0

--- Main ---

main :: IO ()
main = do
    --- Dirichlet simulation

    dsmps <- realize $ sample dnsmps dtru

    let drchs = take nepchs $ fitDirichlet dsmps
        dlls = logLikelihood dsmps <$> drchs
        drch = last drchs
        dxrng = range dmn dmx pltsmps
        dyrng = range dmn dmx pltsmps

    let ddnss = do
            y <- dyrng
            x <- dxrng
            return (density2d dtru (x, y), density2d drch (x, y))

    let (dtdns, dldns) = unzip ddnss

    let djson =
            toJSON
                [ "log-likelihood" .= dlls
                , "xrange" .= dxrng
                , "yrange" .= dyrng
                , "true-density" .= dtdns
                , "learned-density" .= dldns
                , "samples" .= map S.toList dsmps
                ]

    --- Multivariate Normal simulation
    nsmps <- realize $ sample nnsmps ntru

    let mmvn :: Mean # FitNormal
        mmvn = averageSufficientStatistic nsmps
        smvn = toSource mmvn
        nmvn = toNatural mmvn
        nxrng = range nxmn nxmx pltsmps
        nyrng = range nymn nymx pltsmps

    let ndnss = do
            y <- nyrng
            x <- nxrng
            return (density ntru $ S.doubleton x y, density smvn $ S.doubleton x y, density nmvn $ S.doubleton x y)

    let (ntdns, nsdns, nndns) = unzip3 ndnss

    let njson =
            toJSON
                [ "xrange" .= nxrng
                , "yrange" .= nyrng
                , "true-density" .= ntdns
                , "source-fit-density" .= nsdns
                , "natural-fit-density" .= nndns
                , "samples" .= map S.toList nsmps
                ]

    --- Process data
    mvnfl <- resultsFilePath "multivariate-normal.json"
    drchfl <- resultsFilePath "multivariate-dirichlet.json"

    exportJSON mvnfl njson
    exportJSON drchfl djson

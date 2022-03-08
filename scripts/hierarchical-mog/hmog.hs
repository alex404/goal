#! stack runghc

{-# LANGUAGE ScopedTypeVariables,TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds,KindSignatures,TypeFamilies,NoStarIsType,UndecidableInstances
   #-}

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G


import qualified Data.List as L

--- Globals ---


-- Types --

type HMOG n m k = AffineMixture (MultivariateNormal m) (LinearGaussianHarmonium n m) k
type HMOG2 n m k = AffineHarmonium Tensor (MVNMean n) (MVNMean m) (MultivariateNormal n) (Mixture (MultivariateNormal m) k)

type N = 10
type M = 2
type K = 2

type HMOG' = HMOG N M K
type HMOG2' = HMOG2 N M K

-- Data Processing --

-- Reproduced tutorial here: https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/factor-analysis/A-simple-example-of-FA/index.html
-- Data file pulled from:
-- https://userpage.fu-berlin.de/soga/300/30100_data_sets/food-texture.csv
ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/kernels_selected"


main :: IO ()
main = do

    print $ natValInt (Proxy @(Dimension HMOG'))
    print $ natValInt (Proxy @(Dimension HMOG2'))
    print $ natValInt (Proxy @(Dimension (LinearGaussianHarmonium N M)))
    print $ natValInt (Proxy @(Dimension (Mixture (MultivariateNormal M) K)))
    print $ natValInt (Proxy @(Dimension (MultivariateNormal M)))

    lrxss <- goalImport csvpth
    let xss :: [[Double]]
        xss = case lrxss of
                Left err -> error err
                Right csvs -> csvs

    print $ length xss
    print . L.nub $ length <$> xss

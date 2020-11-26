#! stack runghc

{-# LANGUAGE
    DataKinds,
    ScopedTypeVariables,
    DeriveGeneric,
    TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Data.List as L


--- Program ---


-- Globals --

nsmps :: Int
nsmps = 1000

mu,nu :: Double
mu = 20
nu = 0.3

tru :: Source # CoMPoisson
tru = fromTuple (mu,nu)

-- Plot

mx :: Int
mx = 50

xs :: Sample CoMPoisson
xs = [0..mx]

-- Main --

main :: IO ()
main = do

    let [s1,s2] = listCoordinates $ toMean tru

    smps <- realize $ sample nsmps tru

    let s1',s2' :: Double
        s1' = average $ fromIntegral <$> smps
        s2' = average $ logFactorial <$> smps

    putStrLn $ "\nTrue Mean:" ++ take 5 (show s1)
    putStrLn $ "True Average Log-Factorial:" ++ take 5 (show s2)

    putStrLn $ "\nSample Mean:"  ++ take 5 (show s1')
    putStrLn $ "Sample Average Log-Factorial:" ++ take 5 (show s2')

    let hsts0 = fromIntegral . subtract 1 . L.length <$> L.group (L.sort $ smps ++ xs)
        hsts = (/(sum hsts0 :: Double)) <$> hsts0

    let ldpth = "data"

    goalExport ldpth "histogram" $ zip xs hsts
    goalExport ldpth "densities" . zip xs $ density (toNatural tru) <$> xs

    runGnuplot ldpth "com-poisson"

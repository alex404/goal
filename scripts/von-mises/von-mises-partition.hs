#! stack runghc

{-# LANGUAGE
    DataKinds,
    TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified Criterion.Main as C


--- Program ---


-- Globals --

nsmps :: Int
nsmps = 100000

mu,kap :: Double
mu = 2
kap = 2

tru :: Source # VonMises
tru = Point $ S.doubleton mu kap

-- Plot


potential' :: Natural # VonMises -> Double
potential' nvm =
    logIntegralExp 1e-6  (head . unnormalizedLogDensities nvm . (:[])) 0 (2*pi) (range 0 (2*pi) 100)

vonMisesFromPrecision :: Double -> Natural # VonMises
vonMisesFromPrecision prcs =
    let svm  :: Source # VonMises
        svm = Point $ S.doubleton pi prcs
     in toNatural svm

testAtPrecision :: Double -> IO ()
testAtPrecision prcs = do
    let nvm = vonMisesFromPrecision prcs
    putStrLn $ "\nPrecision: " ++ show prcs
    putStrLn $ "GSL Log-Partition Function: " ++ show (potential nvm)
    putStrLn $ "LIE Log-Partition Function: " ++ show (potential' nvm)

-- Main --

main :: IO ()
main = do

    --mapM_ testAtPrecision [0,1,1e2,1e-1,1e-2,1e-3,1e-4,1e-5]
    let prcsns = [0,1,1e2,1e-1,1e-2,1e-3,1e-4,1e-5]
        nvms = vonMisesFromPrecision <$> prcsns

    C.defaultMain
       [ C.bench "GSL" $ C.nf (sum . map potential) nvms
       , C.bench "LIE" $ C.nf (sum . map potential') nvms ]

#! stack runghc

{-# LANGUAGE ScopedTypeVariables,TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds
   #-}


--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S


--- Globals ---


-- Data file pulled from:
-- https://userpage.fu-berlin.de/soga/300/30100_data_sets/food-texture.csv
ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/food-texture.dat"

-- Initialization --

mvx :: Natural # Replicated 5 Normal
mvx = joinReplicated $ S.zipWith (curry fromTuple) 0 (-0.1)

fa0 :: Natural # FactorAnalysis 5 2
fa0 = join mvx (-0.1)

-- Training --

nepchs :: Int
nepchs = 20

-- Functions --

parseCSV :: String -> Sample (MultivariateNormal 5)
parseCSV csvstr = do
    csv <- tail $ lines csvstr
    let lst = read $ '[' : drop 5 csv ++ "]"
    return . fromJust $ S.fromList lst

getCorrelations
    :: Source # MultivariateNormal 5
    -> [Double]
getCorrelations mvn = do
    let crrs = multivariateNormalCorrelations mvn
    (idx,crw) <- zip [1 :: Int ..] . S.toList $ S.toRows crrs
    drop idx $ S.toList crw


--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let smps = parseCSV csvstr

    let emfas = take nepchs $ iterate (factorAnalysisExpectationMaximization smps) fa0

    let lls = do
            mvnz <- factorAnalysisObservableDistribution <$> emfas
            return $ logLikelihood smps mvnz

    putStrLn "LL Ascent:"
    mapM_ print lls


    let mvn :: Source # MultivariateNormal 5
        mvn = mle smps
    let fa = last emfas
        famvn = factorAnalysisObservableDistribution fa
        crrs = getCorrelations mvn
        facrrs = getCorrelations $ transition famvn

    goalExport ldpth "correlations" $ zip crrs facrrs
    runGnuplot ldpth "correlations-scatter"

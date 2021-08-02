#! stack runghc

{-# LANGUAGE ScopedTypeVariables,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds
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
csvpth :: FilePath
csvpth = "data/food-texture.csv"


-- Initialization --

mvx :: Natural # MultivariateNormal 5
mvx = joinNaturalMultivariateNormal 0 (-1)

mvz :: Natural # MultivariateNormal 2
mvz = joinNaturalMultivariateNormal 0 (-1)

fa0 :: Natural # FactorAnalysis 5 2
fa0 = fst . split $ joinHarmonium mvx (-1) mvz

-- Training --

nepchs :: Int
nepchs = 20

-- Functions --

parseCSV :: String -> Sample (MultivariateNormal 5)
parseCSV csvstr = do
    csv <- tail $ lines csvstr
    let lst = read $ '[' : drop 5 csv ++ "]"
    return . fromJust $ S.fromList lst


--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let smps = parseCSV csvstr

    let emfas = take nepchs $ iterate (factorAnalysisExpectationMaximization smps) fa0

    let lls = do
            mvnz <- factorAnalysisObservableDistribution <$> emfas
            return $ logLikelihood smps mvnz

    print lls

    let fa = last emfas
        mnvz = factorAnalysisObservableDistribution fa
        (muz,sgmaz) = splitMultivariateNormal $ transition mnvz
    print muz
    print sgmaz

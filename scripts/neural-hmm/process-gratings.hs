#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE ScopedTypeVariables #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Graphical

import qualified Data.Vector.Storable as S
--import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L
import qualified Data.Map as M


listEquals :: Eq x => [x] -> Bool
listEquals [] = True
listEquals (x:xs) = all (== x) xs

--- Main ---

nstps :: Int
nstps = 4

toResponse :: [Double] -> S.Vector Int
toResponse xs =
    S.fromList $ round <$> xs

exportData :: Int -> IO ()
exportData n = do

    let flnm = "data/gratings-session-" ++ show n

    putStrLn $ "\nProcessing File: " ++ flnm

    estr <- goalImport flnm
    csvs :: [[Double]]
        <- case estr of
          Left csv -> error csv
          Right csv -> return csv

    let cntrsts0:coris0:soris0:stps:rspnss0 = L.transpose csvs
        rspnss = breakEvery nstps $ toResponse <$> L.transpose rspnss0
        cntrsts = head <$> breakEvery nstps cntrsts0
        coris = head <$> breakEvery nstps coris0
        soris = head <$> breakEvery nstps soris0

    putStrLn "\nSane Subsequences?"
    print . all (== [1..fromIntegral nstps]) $ breakEvery nstps stps
    print . all listEquals $ breakEvery nstps cntrsts0
    print . all listEquals $ breakEvery nstps coris0
    print . all listEquals $ breakEvery nstps soris0

    let dmp = conditionalDataMap $ zip rspnss (zip3 cntrsts coris soris)
        (kys,elms) = unzip $ M.toAscList dmp
        kyprs = zip [0 :: Int ..] kys

    putStrLn "\nStimulus Indices:"
    mapM_ print kyprs

    let dmp' = M.fromAscList $ zip (fst <$> kyprs) elms
        kmp = M.fromAscList kyprs

    writeFile ("data/gratings-session-" ++ show n ++ ".map") $ show dmp'
    writeFile ("data/gratings-session-" ++ show n ++ ".keys") $ show kmp

main :: IO ()
main = do

    let ns = [1,2,3]

    mapM_ exportData ns


{-# LANGUAGE ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core

import Data.List
import qualified Data.Map as M


--- Globals ---


--- Main ---


poolData :: KohnExperiment nn -> [String] -> IO ()
poolData kxp exps = do

    let kpdr = kohnProjectPath kxp

    let poolPairs ex (x,nmp) = (x, M.mapKeys (toPooledNeuronID ex) nmp)

    dts <- sequence $ do

        ex <- exps

        let expdr = "kohn-data/" ++ protocol kxp ++ "/" ++ ex

        return $ do

            (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile expdr "stmstrm0"
            (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile expdr "stmstrm1"
            (stmttls0 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile expdr "stmttls0"
            (stmttls1 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile expdr "stmttls1"

            return ( poolPairs ex <$> stmstrm0, poolPairs ex <$> stmstrm1
                   , poolPairs ex <$> stmttls0, poolPairs ex <$> stmttls1 )

    let (stmstrm0s,stmstrm1s,stmttls0s,stmttls1s) = unzip4 dts

    let plstmstrm0 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm0s
    let plstmstrm1 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm1s
        plstmttls0 = foldl1 (M.unionWith (\(k,mp1) (_,mp2) -> (k,M.union mp1 mp2))) stmttls0s
        plstmttls1 = foldl1 (M.unionWith (\(k,mp1) (_,mp2) -> (k,M.union mp1 mp2))) stmttls1s

    let nrns = M.keys . snd . snd . head . M.toList $ plstmttls0

    putStrLn $ "\nPROTOCOL: " ++ protocol kxp ++ " | EXPERIMENT: " ++ experiment kxp ++ "\n"

    putStr "Number of Neurons: "
    print $ length nrns
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length plstmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length plstmstrm1

    goalWriteFile kpdr "stmstrm0" $ show plstmstrm0
    goalWriteFile kpdr "stmstrm1" $ show plstmstrm1
    goalWriteFile kpdr "stmttls0" $ show plstmttls0
    goalWriteFile kpdr "stmttls1" $ show plstmttls1

main :: IO ()
main =

    poolData small40Pooled ["112l44", "112l45", "112r35", "112r36" ]
    --poolData big40Pooled ["105r62", "107l114", "112l16", "112r32" ]


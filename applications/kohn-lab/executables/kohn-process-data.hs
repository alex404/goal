{-# LANGUAGE ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core
import qualified Goal.Core.Vector.Boxed as B

import Data.List
import qualified Data.Map as M


--- Globals ---


--- Main ---

kpp :: String
kpp = "patterson-2013"

poolData :: forall nn . KnownNat nn => KohnExperiment nn -> [String] -> IO ()
poolData kxp exps = do

    let poolPairs ex (x,nmp) = (x, M.mapKeys (toPooledNeuronID ex) nmp)

    dts <- sequence $ do

        ex <- exps

        let expdr = "patterson-2013/" ++ protocol kxp ++ "/" ++ ex

        return $ do

            (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile expdr "stmstrm0"
            (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile expdr "stmstrm1"
            (stmttls0 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile expdr "stmttls0"
            (stmttls1 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile expdr "stmttls1"

            return ( poolPairs ex <$> stmstrm0, poolPairs ex <$> stmstrm1
                   , poolPairs ex <$> stmttls0, poolPairs ex <$> stmttls1 )

    let (stmstrm0s,stmstrm1s,stmttls0s,_) = unzip4 dts

    let plstmstrm0 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm0s
    let plstmstrm1 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm1s
        plstmttls0 = foldl1 (M.unionWith (\(k,mp1) (_,mp2) -> (k,M.union mp1 mp2))) stmttls0s
        --plstmttls1 = foldl1 (M.unionWith (\(k,mp1) (_,mp2) -> (k,M.union mp1 mp2))) stmttls1s

    let nrns = M.keys . snd . snd . head . M.toList $ plstmttls0

    putStrLn $ "\nPROTOCOL: " ++ protocol kxp ++ " | EXPERIMENT: " ++ experiment kxp ++ "\n"

    putStr "Number of Neurons: "
    print $ length nrns
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length plstmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length plstmstrm1

    --goalWriteFile kpdr "stmstrm0" $ show plstmstrm0
    --goalWriteFile kpdr "stmstrm1" $ show plstmstrm1
    --goalWriteFile kpdr "stmttls0" $ show plstmttls0
    --goalWriteFile kpdr "stmttls1" $ show plstmttls1

    let predts = Dataset $ experiment kxp ++ "-pre-adapt"
        pstdts = Dataset $ experiment kxp ++ "-post-adapt"

    let zxs0,zxs1 :: [(B.Vector nn Int, Stimulus)]
        zxs0 = converter <$> plstmstrm0
        zxs1 = converter <$> plstmstrm1

    goalWriteDataset kpp predts $ show zxs0
    goalWriteDataset kpp pstdts $ show zxs1


processData :: forall nn . KnownNat nn => KohnExperiment nn -> IO ()
processData kxp = do

    let kpdr = kohnProjectPath kxp

    putStrLn $ "\nPROTOCOL: " ++ protocol kxp ++ " | EXPERIMENT: " ++ experiment kxp ++ "\n"

    ecss <- getSpikes kxp
    bids <- getBIDs kxp
    chns <- getChannels kxp
    adpt <- getAdaptor kxp
    --let chns = Nothing

    let strm0 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        strm0 = blockStream chns bids ecss

    let prtclttl,prtcln1 :: Int
        prtclttl = length $ dropWhile ((/= 0) .  fst . fst) $ reverse strm0
        prtcln1 = length $ takeWhile ((/= 0) .  fst . fst) strm0

    let bidstrm,bidstrm0,bidstrm1 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        bidstrm = take prtclttl strm0
        bidstrm0 = take prtcln1 bidstrm
        bidstrm1 = drop (prtcln1 + 1) bidstrm

    let stmstrm0,stmstrm1 :: [(Stimulus, M.Map NeuronID [SpikeTime])]
        stmstrm0 = blockToStimulusStream adpt bidstrm0
        stmstrm1 = blockToStimulusStream adpt bidstrm1

    let zxs0,zxs1 :: [(B.Vector nn Int, Stimulus)]
        zxs0 = converter <$> stmstrm0
        zxs1 = converter <$> stmstrm1

    let bidttls0,bidttls1 :: M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
        bidttls0 = blockIDTotals bids bidstrm0
        bidttls1 = blockIDTotals bids bidstrm1

    let stmttls0,stmttls1 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])
        stmttls0 = blockIDToStimulusTotals adpt bidttls0
        stmttls1 = blockIDToStimulusTotals adpt bidttls1

    let nrns = M.keys . snd . head $ bidstrm

    putStr "Number of Neurons: "
    print $ length nrns
    putStr "Adaptor: "
    print adpt
    putStr "Adaptor (Radians): "
    print $ adaptorToRads adpt
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length stmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length stmstrm1
    putStrLn "Block ID Trial Counts: "
    print . map length . group . sort $ fst . fst <$> bidstrm
    putStrLn "Pre Block ID Trial Counts: "
    print . map length . group . sort $ fst . fst <$> bidstrm0
    putStrLn "Post Block ID Trial Counts: "
    print . map length . group . sort $ fst . fst <$> bidstrm1

    --goalWriteFile kpdr "bidstrm" $ show bidstrm
    --goalWriteFile kpdr "bidstrm0" $ show bidstrm0
    --goalWriteFile kpdr "bidstrm1" $ show bidstrm1
    --goalWriteFile kpdr "bidttls0" $ show bidttls0
    --goalWriteFile kpdr "bidttls1" $ show bidttls1
    goalWriteFile kpdr "stmstrm0" $ show stmstrm0
    goalWriteFile kpdr "stmstrm1" $ show stmstrm1
    goalWriteFile kpdr "stmttls0" $ show stmttls0
    goalWriteFile kpdr "stmttls1" $ show stmttls1
    let predts = Dataset $ experiment kxp ++ "-pre-adapt"
        pstdts = Dataset $ experiment kxp ++ "-post-adapt"
    goalWriteDataset kpp predts $ show zxs0
    goalWriteDataset kpp pstdts $ show zxs1

main :: IO ()
main = do

    processData experiment112l44
    processData experiment112l45
    processData experiment112r35
    processData experiment112r36
    processData experiment105r62
    processData experiment107l114
    processData experiment112l16
    processData experiment112r32
    poolData small40Pooled ["112l44", "112l45", "112r35", "112r36" ]

    let dst = Dataset <$>
            [ experiment experiment112l44
            , experiment experiment112l45
            , experiment experiment112r35
            , experiment experiment112r36
            , experiment experiment105r62
            , experiment experiment107l114
            , experiment experiment112l16
            , experiment experiment112r32
            , experiment small40Pooled ]

    goalWriteCSV kpp "datasets" dst
    putStrLn "\n"

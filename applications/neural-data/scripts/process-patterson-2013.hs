{-# LANGUAGE ScopedTypeVariables,DataKinds #-}


--- Imports ---


-- Goal --

import Goal.Core

import NeuralData

-- Other --

import qualified Data.Map as M
import qualified Data.List as L


--- Globals ---


--- Experiments ---


-- Type Synonyms --

type BlockID = Int
type Stimulus = Double

-- Types --

data NeuronID =
    NeuronID (Int,Int) | PooledNeuronID (String, Int, Int)
    deriving (Eq, Read, Show, Ord)

-- Patterson Lab Record --

data PattersonExperiment = PattersonExperiment
    { protocol :: String
    , experiment :: String }

-- Experiments --

experiment112l44 :: PattersonExperiment
experiment112l44 = PattersonExperiment "small40" "112l44"

experiment112l45 :: PattersonExperiment
experiment112l45 = PattersonExperiment "small40" "112l45"

experiment112r35 :: PattersonExperiment
experiment112r35 = PattersonExperiment "small40" "112r35"

experiment112r36 :: PattersonExperiment
experiment112r36 = PattersonExperiment "small40" "112r36"

experiment105r62 :: PattersonExperiment
experiment105r62 = PattersonExperiment "big40" "105r62"

experiment107l114 :: PattersonExperiment
experiment107l114 = PattersonExperiment "big40" "107l114"

experiment112l16 :: PattersonExperiment
experiment112l16 = PattersonExperiment "big40" "112l16"

experiment112r32 :: PattersonExperiment
experiment112r32 = PattersonExperiment "big40" "112r32"

big40Pooled :: PattersonExperiment
big40Pooled = PattersonExperiment "big40" "big40-pooled"

small40Pooled :: PattersonExperiment
small40Pooled = PattersonExperiment "small40" "small40-pooled"


--- Functions ---


converter :: (Stimulus,M.Map NeuronID Int) -> ([Int],Stimulus)
converter (s,mp) = (M.elems mp,s)


-- Sample Stream Builder --

nullNeuronMap :: [(Int,Int,Double)] -> M.Map NeuronID Int
nullNeuronMap ecss =
    let ecs = [ NeuronID (e,c) | (e,c,_) <- ecss, e /= 2000 ]
     in M.fromList . zip (L.nub ecs) $ repeat 0

-- | Filter out weird channels, and cut off the initial (meaningless?) part of the data stream.
preFilterSpikes :: Maybe [Int] -> [(Int,Int,Double)] -> [(Int,Int,Double)]
preFilterSpikes mchns ecss =
    let ecss' = dropWhile (\(e,_,_) -> e /= 2000) ecss
     in filter (\(e,c,_) -> e == 2000 || (maybe True (elem e) mchns && 0 < c && c < 10)) ecss'

-- | Breaks the data stream into a stream of pairs of blockIDs + times, and the stream of spikes
-- that happened in that block.
breakSpikeBlocks :: [BlockID] -> [(Int,Int,Double)] -> [(BlockID,[NeuronID])]
breakSpikeBlocks [] _ = []
breakSpikeBlocks (bid:bids) (_:ecss) =
    let (ecss1,ecss') = span (\(e,_,_) -> e /= 2000) ecss
     in ( bid, [ NeuronID (e,c) | (e,c,_) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: Maybe [Int] -> [Int] -> [(Int,Int,Double)] -> [(BlockID,M.Map NeuronID Int)]
blockStream mchns bids ecss =
    let ecss' = preFilterSpikes mchns ecss
        mp0 = nullNeuronMap ecss'
        spkblks = breakSpikeBlocks bids ecss'
     in [ (bid, flip M.union mp0 . M.fromListWith (+) $ zip nidspks $ repeat 1)
            | (bid,nidspks) <- spkblks ]

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockToStimulusStream
    :: Maybe Double
    -> [(BlockID,M.Map NeuronID Int)]
    -> [(Stimulus,M.Map NeuronID Int)]
blockToStimulusStream madpt = mapMaybe patternMatch
    where patternMatch (bid,mp)
            | bid == 0 || bid == 1 = Nothing
            | bid <= 9 = Just (bidToStimulus madpt bid,mp)
            | bid <= 17 = patternMatch (bid-8,mp)
            | otherwise = Nothing

adaptorToRads :: Double -> Double
adaptorToRads adpt =
    let adptpi0 = 2*pi*adpt/360
     in 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

bidToStimulus :: Maybe Double -> BlockID -> Stimulus
bidToStimulus madpt bid =
    let adpt' = maybe 0 maybeer madpt
        bid' = circulate $ bid - 2
     in ((pi/4)*) . fromIntegral . circulate $ bid' - adpt' + 4
    where maybeer adpt = circulate . round $ adpt / 22.5
          circulate k
            | k < 0 = circulate $ k + 8
            | k > 7 = circulate $ k - 8
            | otherwise = k

--bidToCentredStimulus :: Double -> BlockID -> Stimulus
--bidToCentredStimulus adpt bid =
--    let adpt' = circulate . round $ adpt / 22.5
--        bid' = circulate $ bid - 2
--     in ((pi/4)*) . fromIntegral . circulate $ bid' - adpt' + 4
--    where circulate k
--            | k < 0 = circulate $ k + 8
--            | k > 7 = circulate $ k - 8
--            | otherwise = k
--


--- IO ---


pattersonPath :: FilePath
pattersonPath = "patterson-2013"

expmnt :: Experiment
expmnt = Experiment prjnm pattersonPath

importPattersonData :: FromRecord r => PattersonExperiment -> String -> IO [r]
importPattersonData pxp flnm = do
    let flpth = concat [experiment pxp, "/", flnm]
    goalImport expmnt flpth

getBIDs :: PattersonExperiment -> IO [Int]
getBIDs pxp = map head <$> importPattersonData pxp "blockIDs"

getSpikes :: PattersonExperiment -> IO [(Int,Int,Double)]
getSpikes pxp = importPattersonData pxp "spikes"

--getChannels :: PattersonExperiment -> IO (Maybe [Int])
--getChannels pxp = do
--
--    csvpth <- ( ++ "/channels.csv") <$> pattersonRawDataPath pxp
--    bl <- doesFileExist csvpth
--
--    if bl
--       then do
--           chnstr <- readFile csvpth
--           return . Just . map read $ lines chnstr
--       else return Nothing

getAdaptor :: PattersonExperiment -> IO Double
getAdaptor pxp = head . head <$> importPattersonData pxp "adaptor"

--- Main ---

poolData
    :: PattersonExperiment
    -> [([(Stimulus,M.Map NeuronID Int)],[(Stimulus,M.Map NeuronID Int)])]
    -> IO ()
poolData pxp stmstrms = do

    let (stmstrm0s,stmstrm1s) = unzip stmstrms
        stmstrm0s' = L.sortOn fst <$> stmstrm0s
        stmstrm1s' = L.sortOn fst <$> stmstrm1s
        stm0s = fst <$> head stmstrm0s'
        stm1s = fst <$> head stmstrm1s'
        stmsstrm0s'' = map snd <$> stmstrm0s'
        stmsstrm1s'' = map snd <$> stmstrm1s'

    putStrLn "\nNumber of Neurons?"
    print . length . L.nub . concat $ concatMap (M.keys . snd) <$> stmstrm0s
    print . sum $ length . M.keys . snd . head <$> stmstrm0s

    let zxs0 = flip zip stm0s $ concatMap M.elems <$> L.transpose stmsstrm0s''
    let zxs1 = flip zip stm1s $ concatMap M.elems <$> L.transpose stmsstrm1s''
    let k = length . fst $ head zxs0

    putStrLn "\nPooled:"
    putStr "Number of Pooled Neurons: "
    print k
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length zxs0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length zxs1

    let predts = experiment pxp ++ "-pre-adapt"
        pstdts = experiment pxp ++ "-post-adapt"

    goalWriteDataset expmnt predts $ show (k,zxs0)
    goalWriteDataset expmnt pstdts $ show (k,zxs1)

processData
    :: PattersonExperiment
    -> IO ([(Stimulus,M.Map NeuronID Int)],[(Stimulus,M.Map NeuronID Int)])
processData pxp = do

    putStrLn $ "\nPROTOCOL: " ++ protocol pxp ++ " | EXPERIMENT: " ++ experiment pxp ++ "\n"

    ecss <- getSpikes pxp
    bids <- getBIDs pxp
    --chns <- getChannels pxp
    adpt <- getAdaptor pxp
    --let chns = Nothing

    let strm0 :: [(BlockID, M.Map NeuronID Int)]
        strm0 = blockStream Nothing bids ecss

    let prtclttl,prtcln1 :: Int
        prtclttl = length $ dropWhile ((/= 0) . fst) $ reverse strm0
        prtcln1 = length $ takeWhile ((/= 0) . fst) strm0

    let bidstrm,bidstrm0,bidstrm1 :: [(BlockID, M.Map NeuronID Int)]
        bidstrm = take prtclttl strm0
        bidstrm0 = take prtcln1 bidstrm
        bidstrm1 = drop (prtcln1 + 1) bidstrm

    let stmstrm0,stmstrm1 :: [(Stimulus, M.Map NeuronID Int)]
        stmstrm0 = blockToStimulusStream (Just adpt) bidstrm0
        stmstrm1 = blockToStimulusStream (Just adpt) bidstrm1

    let zxs0,zxs1 :: [([Int],Stimulus)]
        zxs0 = converter <$> stmstrm0
        zxs1 = converter <$> stmstrm1

    let k = length . M.keys . snd . head $ bidstrm

    putStr "Number of Neurons: "
    print k
    putStr "Adaptor: "
    print adpt
    putStr "Adaptor (Radians): "
    print $ adaptorToRads adpt
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length stmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length stmstrm1
    putStrLn "Block ID Trial Counts: "
    print . map length . L.group . L.sort $ fst <$> bidstrm
    putStrLn "Pre Block ID Trial Counts: "
    print . map length . L.group . L.sort $ fst <$> bidstrm0
    putStrLn "Post Block ID Trial Counts: "
    print . map length . L.group . L.sort $ fst <$> bidstrm1

    let predts = experiment pxp ++ "-pre-adapt"
        pstdts = experiment pxp ++ "-post-adapt"

    goalWriteDataset expmnt predts $ show (k,zxs0)
    goalWriteDataset expmnt pstdts $ show (k,zxs1)

    return (stmstrm0, stmstrm1)

main :: IO ()
main = do

    stmstrms <- mapM processData
        [ experiment112l44, experiment112l45, experiment112r35, experiment112r36
        , experiment105r62, experiment107l114, experiment112l16, experiment112r32 ]
    poolData small40Pooled $ take 4 stmstrms



    let dsts = concat $ do
            xp <- [ experiment experiment112l44, experiment experiment112l45
                  , experiment experiment112r35, experiment experiment112r36
                  , experiment experiment105r62, experiment experiment107l114
                  , experiment experiment112l16, experiment experiment112r32
                  , experiment small40Pooled ]
            return [xp ++ "-pre-adapt", xp ++ "-post-adapt"]

    goalWriteDatasetsCSV expmnt dsts
    putStrLn "\n"

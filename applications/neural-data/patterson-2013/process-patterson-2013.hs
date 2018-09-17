{-# LANGUAGE ScopedTypeVariables,DataKinds #-}


--- Imports ---


-- Goal --

import Goal.Core

-- Other --

import qualified Data.Map as M
import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V
import Data.List
import qualified Data.Csv as CSV


--- Globals ---


--- Experiments ---


-- Type Synonyms --

type BlockID = Int
type BlockEvent = (BlockID,SpikeTime)
type Stimulus = Double
type SpikeTime = Double

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

--big40Pooled :: PattersonExperiment
--big40Pooled = PattersonExperiment "big40" "big40-pooled"

small40Pooled :: PattersonExperiment
small40Pooled = PattersonExperiment "small40" "small40-pooled"



--- Functions ---

converter :: (Stimulus,M.Map NeuronID [SpikeTime]) -> ([Int],Stimulus)
converter (s,mp) =
    (length <$> M.elems mp,s)



-- Sample Stream Builder --

nullNeuronMap :: Monoid a => [(Int,Int,Double)] -> M.Map NeuronID a
nullNeuronMap ecss =
    let ecs = [ NeuronID (e,c) | (e,c,_) <- ecss, e /= 2000 ]
     in M.fromList . zip (nub ecs) $ repeat mempty

-- | Filter out weird channels, and cut off the initial (meaningless?) part of the data stream.
preFilterSpikes :: Maybe [Int] -> [(Int,Int,Double)] -> [(Int,Int,Double)]
preFilterSpikes mchns ecss =
    let ecss' = dropWhile (\(e,_,_) -> e /= 2000) ecss
     in filter (\(e,c,_) -> e == 2000 || (maybe True (elem e) mchns && 0 < c && c < 10)) ecss'

-- | Breaks the data stream into a stream of pairs of blockIDs + times, and the stream of spikes
-- that happened in that block.
breakSpikeBlocks :: [BlockID] -> [(Int,Int,Double)] -> [(BlockEvent,[(NeuronID,[SpikeTime])])]
breakSpikeBlocks [] _ = []
breakSpikeBlocks (bid:bids) ((_,_,s0):ecss) =
    let (ecss1,ecss') = span (\(e,_,_) -> e /= 2000) ecss
     in ( (bid,s0), [ (NeuronID (e,c),[s]) | (e,c,s) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: Maybe [Int] -> [Int] -> [(Int,Int,Double)] -> [(BlockEvent,M.Map NeuronID [SpikeTime])]
blockStream mchns bids ecss =
    let ecss' = preFilterSpikes mchns ecss
        mp0 = nullNeuronMap ecss'
        bnspks = breakSpikeBlocks bids ecss'
     in [ (bspk, flip M.union mp0 $ M.fromListWith (++) nspks) | (bspk,nspks) <- bnspks ]

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockToStimulusStream :: Double -> [(BlockEvent,M.Map NeuronID [SpikeTime])] -> [(Stimulus,M.Map NeuronID [SpikeTime])]
blockToStimulusStream adpt = mapMaybe patternMatch
    where patternMatch ((k,tm),mp)
            | k == 0 || k == 1 = Nothing
            | k <= 9 = Just (bidToCentredStimulus adpt k,mp)
            | k <= 17 = patternMatch ((k-8,tm),mp)
            | otherwise = Nothing


--- Plots ---


-- Tuning Curves --

-- Fit --

adaptorToRads :: Double -> Double
adaptorToRads adpt =
    let adptpi0 = 2*pi*adpt/360
     in 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

bidToCentredStimulus :: Double -> BlockID -> Stimulus
bidToCentredStimulus adpt bid =
    let adpt' = circulate . round $ adpt / 22.5
        bid' = circulate $ bid - 2
     in ((pi/4)*) . fromIntegral . circulate $ bid' - adpt' + 4
    where circulate k
            | k < 0 = circulate $ k + 8
            | k > 7 = circulate $ k - 8
            | otherwise = k


--- IO ---

pattersonRawDataPath :: PattersonExperiment -> FilePath
pattersonRawDataPath kd = "raw-data/" ++ protocol kd ++ "/" ++ experiment kd

getBIDs :: PattersonExperiment -> IO [Int]
getBIDs pxp = do
    let csvpth = pattersonRawDataPath pxp ++ "/blockIDs.csv"
    bidstr <- readFile csvpth
    return $ read <$> lines bidstr

getSpikes :: PattersonExperiment -> IO [(Int,Int,Double)]
getSpikes pxp = do

    let csvpth = pattersonRawDataPath pxp ++ "/spikes.csv"
    ecsstr <- BS.readFile csvpth
    let (Right ecssV) = C.decode C.NoHeader ecsstr
    return $ V.toList ecssV

getChannels :: PattersonExperiment -> IO (Maybe [Int])
getChannels pxp = do

    let csvpth = pattersonRawDataPath pxp ++ "/channels.csv"
    bl <- doesFileExist csvpth

    if bl
       then do
           chnstr <- readFile csvpth
           return . Just . map read $ lines chnstr
       else return Nothing

getAdaptor :: PattersonExperiment -> IO Double
getAdaptor pxp = do
    let csvpth = pattersonRawDataPath pxp ++ "/adaptor.csv"
    adpstr <- readFile csvpth
    return . head $ read <$> lines adpstr

--- Main ---

poolData
    :: PattersonExperiment
    -> [([(Stimulus,M.Map NeuronID [SpikeTime])],[(Stimulus,M.Map NeuronID [SpikeTime])])]
    -> IO ()
poolData pxp stmstrms = do

    let (stmstrm0s,stmstrm1s) = unzip stmstrms

    let plstmstrm0 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm0s
    let plstmstrm1 = foldl1 (zipWith (\(stm,mp1) (_,mp2) -> (stm,M.union mp1 mp2)))
            $ sortOn fst <$> stmstrm1s

    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length plstmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length plstmstrm1

    let predts = experiment pxp ++ "-pre-adapt.dat"
        pstdts = experiment pxp ++ "-post-adapt.dat"

    let zxs0 = converter <$> plstmstrm0
        zxs1 = converter <$> plstmstrm1

    writeFile ("data/" ++ predts) $ show zxs0
    writeFile ("data/" ++ pstdts) $ show zxs1


processData
    :: PattersonExperiment
    -> IO ([(Stimulus,M.Map NeuronID [SpikeTime])],[(Stimulus,M.Map NeuronID [SpikeTime])])
processData pxp = do

    putStrLn $ "\nPROTOCOL: " ++ protocol pxp ++ " | EXPERIMENT: " ++ experiment pxp ++ "\n"

    ecss <- getSpikes pxp
    bids <- getBIDs pxp
    chns <- getChannels pxp
    adpt <- getAdaptor pxp
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

    let zxs0 = converter <$> stmstrm0
        zxs1 = converter <$> stmstrm1

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

    let predts = experiment pxp ++ "-pre-adapt.dat"
        pstdts = experiment pxp ++ "-post-adapt.dat"

    writeFile ("data/" ++ predts) $ show zxs0
    writeFile ("data/" ++ pstdts) $ show zxs1

    return (stmstrm0, stmstrm1)

main :: IO ()
main = do

    stmstrms <- mapM processData
        [ experiment112l44, experiment112l45, experiment112r35, experiment112r36
        , experiment105r62, experiment107l114, experiment112l16, experiment112r32 ]
    poolData small40Pooled stmstrms

    let dsts = fmap Dataset . concat $  do
            xp <- [ experiment experiment112l44 , experiment experiment112l45 , experiment experiment112r35
                  , experiment experiment112r36 , experiment experiment105r62 , experiment experiment107l114
                  , experiment experiment112l16 , experiment experiment112r32 , experiment small40Pooled ]
            return [xp ++ "-pre-adapt", xp ++ "-post-adapt"]

    BS.writeFile "datasets.csv" $ CSV.encodeDefaultOrderedByName dsts
    putStrLn "\n"

module Goal.Datasets.KohnLab
    ( -- * Functions
--      timeToSampleStream
--    , blockToTimeStream
--    , streamAverage
      blockStream
    , averageBlockIDs
    , averageBlockIDsToStimuli
    -- * IO
    , getBIDs
    , getSpikes
    , getChannels
    , getAdaptor
    -- * Types
    , BlockID
    , BlockEvent
    , NeuronID
    , SpikeTime
    , SpikeInterval
    , Stimulus
    ) where

import qualified Data.Map as M
import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V


import Goal.Core
import Data.List


--- Functions ---


-- New --

type BlockID = Int
type BlockEvent = (BlockID,SpikeTime)
type NeuronID = (Int,Int)
type Stimulus = Double
type SpikeTime = Double
type SpikeInterval = Double


--- Sample Stream Builder ---

nullNeuronMap :: Monoid a => [(Int,Int,Double)] -> M.Map NeuronID a
nullNeuronMap ecss =
    let ecs = [ (e,c) | (e,c,_) <- ecss, e /= 2000 ]
     in M.fromList . zip (nub ecs) $ repeat mempty

nullBlockIDMap :: Monoid a => [BlockID] -> M.Map BlockID a
nullBlockIDMap bids =
     M.fromList . zip bids $ repeat mempty

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
     in ( (bid,s0), [ ((e,c),[s]) | (e,c,s) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: Maybe [Int] -> [Int] -> [(Int,Int,Double)] -> [(BlockEvent,M.Map NeuronID [SpikeTime])]
blockStream mchns bids ecss =
    let ecss' = preFilterSpikes mchns ecss
        mp0 = nullNeuronMap ecss'
        bnspks = breakSpikeBlocks bids ecss'
     in [(bspk, flip M.union mp0 $ M.fromListWith (++) nspks) | (bspk,nspks) <- bnspks]

averageBlockIDs :: [BlockID] -> [(BlockEvent,M.Map NeuronID [SpikeTime])] -> M.Map BlockID (M.Map NeuronID [SpikeTime])
averageBlockIDs bids bstrm =
    let bstrm' = [(bid,nmp) | ((bid,_),nmp) <- bstrm]
     in flip M.union (nullBlockIDMap bids) $ M.fromListWith (M.unionWith (++)) bstrm'

averageBlockIDsToStimuli :: M.Map BlockID (M.Map NeuronID [SpikeTime]) -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
averageBlockIDsToStimuli bidmp =
    let bidstms = zip [2..9] $ range 0 (2*pi) 8
     in foldr foldfun mempty bidstms
    where foldfun (bid,stm) = M.insert stm (M.unionWith (++) (bidmp M.! bid) (bidmp M.! (bid + 8)))

---- | Converts the stream of blockIDs+times into blockIDs+duration, and shifts the neural spike times
---- to be relative to the blockid time.
--spikeBlockIntervalStream
--    :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
--    -> [((BlockID,SpikeInterval),[(NeuronID,SpikeInterval)])]
--spikeBlockIntervalStream [] = []
--spikeBlockIntervalStream [((bid,s0),nidss)] = [((bid, subtract s0 . snd $ last nidss),neuronIDRelativeSpikes s0 nidss)]
--spikeBlockIntervalStream (((bid,s0),nidss):bnss) =
--    ((bid, subtract s0 . snd . fst $ head bnss),neuronIDRelativeSpikes s0 nidss) : spikeBlockIntervalStream bnss
--
--spikeBlockIntervalStreamToSamples
--    :: [((BlockID,SpikeInterval),[(NeuronID,SpikeInterval)])]
--    -> [(BlockID,[Int])]
--spikeBlockIntervalStreamToSamples strm =
--    let cntmp0 = M.fromList [ (nid,0) | nid <- nub . map fst . concat $ snd <$> strm ]
--     in concat [ zip (repeat bid) $ spikeBlockSlicer cntmp0 (round int) nidss | ((bid,int),nidss) <- strm ]
--
--blockToTimeStream :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])] -> [(BlockID,[Int])]
--blockToTimeStream = spikeBlockIntervalStreamToSamples . spikeBlockIntervalStream
--
--timeToSampleStream :: Double -> [(BlockID,[Int])] -> [(Double,[Int])]
--timeToSampleStream adpt tstrm =
--    let blockIDMapper (bid,spks) = (\alph -> (alph,spks)) <$> blockIDMap adpt bid
--     in mapMaybe blockIDMapper tstrm
--
--streamAverage :: Ord x => [(x,[Int])] -> [(x,[Double])]
--streamAverage tstrm =
--    let bgrps = groupBy (\(bid1,_) (bid2,_) -> bid1 == bid2) (sortBy (comparing fst) tstrm)
--        nrms = genericLength <$> bgrps
--        bids = fst . head <$> bgrps
--        rts = [ map ((/ nrm) . fromIntegral) . foldl1 (zipWith (+)) $ snd <$> bgrp | (nrm,bgrp) <- zip nrms bgrps ]
--     in zip bids rts
--
---- | Converts a block of neuron IDs and spike times into a vector of spike counts.
--spikeBlockToSpikeCounts :: M.Map NeuronID Int -> [(NeuronID,SpikeTime)] -> [Int]
--spikeBlockToSpikeCounts cntmp0 nidss =
--    map snd . M.toAscList $ foldr (\(nid,_) cntmp' -> M.adjust (+1) nid cntmp') cntmp0 nidss
--
--spikeBlockSlicer :: M.Map NeuronID Int -> Int -> [(NeuronID,SpikeInterval)] -> [[Int]]
--spikeBlockSlicer cntmp0 1 nidss = [spikeBlockToSpikeCounts cntmp0 nidss]
--spikeBlockSlicer cntmp0 n nidss =
--    let (nids0,nidss') = span (\(_,s) -> s < 1) nidss
--     in spikeBlockToSpikeCounts cntmp0 nids0 : spikeBlockSlicer cntmp0 (n-1) (neuronIDRelativeSpikes 1 nidss')
--
----- Block Stream Builder --
--
--
---- | Shifts the spike times of the neurons back by the given double.
--neuronIDRelativeSpikes :: Double -> [(NeuronID,SpikeTime)] -> [(NeuronID,SpikeTime)]
--neuronIDRelativeSpikes s0 nidss = [(nid,s - s0) | (nid,s) <- nidss]
--
--
----- BlockID Mapper ---
--
--
--blockIDUnit :: Double
--blockIDUnit = 22.5
--
--blockIDMap0 :: Double -> BlockID -> Maybe Double
--blockIDMap0 adpt 0 = Just adpt
--blockIDMap0 _ 1 = Nothing
--blockIDMap0 _ n
--  | n < 19 = Just $ fromIntegral (n - 2) * blockIDUnit
--  | otherwise = error "BlockID Out of Bounds"
--
----blockIDMap = blockIDMap0
--blockIDMap :: Double -> BlockID -> Maybe Double
--blockIDMap adpt n =
--    mapper <$> blockIDMap0 adpt n
--    where mapper x =
--              let xpi = 2*pi*x/360
--                  k :: Int
--                  k = floor $ xpi / pi
--               in roundSD 5 $ (xpi - fromIntegral k*pi)*2
--
--
--- IO ---


getBIDs :: String -> String -> IO [Int]
getBIDs dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    bidstr <- readFile $ csvdr ++ "blockIDs.csv"
    return $ read <$> lines bidstr

getSpikes :: String -> String -> IO [(Int,Int,Double)]
getSpikes dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    ecsstr <- BS.readFile $ csvdr ++ "spikes.csv"
    let (Right ecssV) = C.decode C.NoHeader ecsstr
    return $ V.toList ecssV

getChannels :: String -> String -> IO (Maybe [Int])
getChannels dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    bl <- doesFileExist $ csvdr ++ "channels.csv"

    if bl
       then do
           chnstr <- readFile $ csvdr ++ "channels.csv"
           return . Just . map read $ lines chnstr
       else return Nothing

getAdaptor :: String -> String -> IO Double
getAdaptor dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    adpstr <- readFile $ csvdr ++ "adaptor.csv"
    return . head $ read <$> lines adpstr




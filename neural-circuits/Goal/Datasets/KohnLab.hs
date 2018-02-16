module Goal.Datasets.KohnLab
    ( -- * Functions
      timeToSampleStream
    , blockToTimeStream
    , streamAverage
    , blockStream
    -- * IO
    , getBIDs
    , getSpikes
    , getChannels
    , getAdaptor
    -- * Types
    , BlockID
    , NeuronID
    , SpikeTime
    , SpikeInterval
    ) where

import qualified Data.Map.Strict as M
import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V


import Goal.Core
import Data.List


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

getChannels :: String -> String -> IO [Int]
getChannels dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    chnstr <- readFile $ csvdr ++ "channels.csv"
    return $ read <$> lines chnstr

getAdaptor :: String -> String -> IO Double
getAdaptor dr flnm = do

    csvdr <- goalDatasetLocation dr flnm

    adpstr <- readFile $ csvdr ++ "adaptor.csv"
    return . head $ read <$> lines adpstr


--- Functions ---


-- New --

type BlockID = Int
type NeuronID = (Int,Int)
type SpikeTime = Double
type SpikeInterval = Double


--- Sample Stream Builder ---

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: [Int] -> [Int] -> [(Int,Int,Double)] -> [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
blockStream chns bids ecss = breakSpikeBlocks bids $ preFilterSpikes chns ecss

blockToTimeStream :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])] -> [(BlockID,[Int])]
blockToTimeStream = spikeBlockIntervalStreamToSamples . spikeBlockIntervalStream

timeToSampleStream :: Double -> [(BlockID,[Int])] -> [(Double,[Int])]
timeToSampleStream adpt tstrm =
    let blockIDMapper (bid,spks) = (\alph -> (alph,spks)) <$> blockIDMap adpt bid
     in mapMaybe blockIDMapper tstrm

streamAverage :: Ord x => [(x,[Int])] -> [(x,[Double])]
streamAverage tstrm =
    let bgrps = groupBy (\(bid1,_) (bid2,_) -> bid1 == bid2) (sortBy (comparing fst) tstrm)
        nrms = genericLength <$> bgrps
        bids = fst . head <$> bgrps
        rts = [ map ((/ nrm) . fromIntegral) . foldl1 (zipWith (+)) $ snd <$> bgrp | (nrm,bgrp) <- zip nrms bgrps ]
     in zip bids rts

-- | Converts the stream of blockIDs+times into blockIDs+duration, and shifts the neural spike times
-- to be relative to the blockid time.
spikeBlockIntervalStream
    :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
    -> [((BlockID,SpikeInterval),[(NeuronID,SpikeInterval)])]
spikeBlockIntervalStream [] = []
spikeBlockIntervalStream [((bid,s0),nidss)] = [((bid, subtract s0 . snd $ last nidss),neuronIDRelativeSpikes s0 nidss)]
spikeBlockIntervalStream (((bid,s0),nidss):bnss) =
    ((bid, subtract s0 . snd . fst $ head bnss),neuronIDRelativeSpikes s0 nidss) : spikeBlockIntervalStream bnss

-- | Converts a block of neuron IDs and spike times into a vector of spike counts.
spikeBlockToSpikeCounts :: M.Map NeuronID Int -> [(NeuronID,SpikeTime)] -> [Int]
spikeBlockToSpikeCounts cntmp0 nidss =
    map snd . M.toAscList $ foldr (\(nid,_) cntmp' -> M.adjust (+1) nid cntmp') cntmp0 nidss

spikeBlockSlicer :: M.Map NeuronID Int -> Int -> [(NeuronID,SpikeInterval)] -> [[Int]]
spikeBlockSlicer cntmp0 1 nidss = [spikeBlockToSpikeCounts cntmp0 nidss]
spikeBlockSlicer cntmp0 n nidss =
    let (nids0,nidss') = span (\(_,s) -> s < 1) nidss
     in spikeBlockToSpikeCounts cntmp0 nids0 : spikeBlockSlicer cntmp0 (n-1) (neuronIDRelativeSpikes 1 nidss')

spikeBlockIntervalStreamToSamples
    :: [((BlockID,SpikeInterval),[(NeuronID,SpikeInterval)])]
    -> [(BlockID,[Int])]
spikeBlockIntervalStreamToSamples strm =
    let cntmp0 = M.fromList [ (nid,0) | nid <- nub . map fst . concat $ snd <$> strm ]
     in concat [ zip (repeat bid) $ spikeBlockSlicer cntmp0 (round int) nidss | ((bid,int),nidss) <- strm ]


--- Block Stream Builder --


-- | Filter out weird channels, and cut off the initial (meaningless?) part of the data stream.
preFilterSpikes :: [Int] -> [(Int,Int,Double)] -> [(Int,Int,Double)]
preFilterSpikes chns ecss =
    let ecss' = dropWhile (\(e,_,_) -> e /= 2000) ecss
     in filter (\(e,c,_) -> e == 2000 || (elem e chns && 0 < c && c < 10)) ecss'

-- | Breaks the data stream into a stream of pairs of blockIDs + times, and the stream of spikes
-- that happened in that block.
breakSpikeBlocks :: [BlockID] -> [(Int,Int,Double)] -> [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
breakSpikeBlocks [] _ = []
breakSpikeBlocks (bid:bids) ((_,_,s0):ecss) =
    let (ecss1,ecss') = span (\(e,_,_) -> e /= 2000) ecss
     in ( (bid,s0), [ ((e,c),s) | (e,c,s) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Shifts the spike times of the neurons back by the given double.
neuronIDRelativeSpikes :: Double -> [(NeuronID,SpikeTime)] -> [(NeuronID,SpikeTime)]
neuronIDRelativeSpikes s0 nidss = [(nid,s - s0) | (nid,s) <- nidss]


--- BlockID Mapper ---


blockIDUnit :: Double
blockIDUnit = 22.5

blockIDMap0 :: Double -> BlockID -> Maybe Double
blockIDMap0 adpt 0 = Just adpt
blockIDMap0 _ 1 = Nothing
blockIDMap0 _ n
  | n < 19 = Just $ fromIntegral (n - 2) * blockIDUnit
  | otherwise = error "BlockID Out of Bounds"

--blockIDMap = blockIDMap0
blockIDMap :: Double -> BlockID -> Maybe Double
blockIDMap adpt n =
    mapper <$> blockIDMap0 adpt n
    where mapper x =
              let xpi = 2*pi*x/360
                  k :: Int
                  k = floor $ xpi / pi
               in roundSD 5 $ (xpi - fromIntegral k*pi)*2

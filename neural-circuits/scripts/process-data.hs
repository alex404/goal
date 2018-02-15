{-# LANGUAGE DataKinds #-}

import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V

import Goal.Core
import Goal.NeuralCircuits
import Goal.Datasets.KohnLab.Old
import Data.List


--- Globals ---

dr,flnm :: String
dr = "adaptation/small40"
flnm = "112l44/"

--- Functions ---


--- Main ---

main :: IO ()
main = do

    csvdr <- goalDatasetLocation dr flnm

    ecsstr <- BS.readFile $ csvdr ++ "spikes.csv"
    let (Right ecssV) = C.decode C.NoHeader ecsstr
        ecss0 :: [(Int,Int,Double)]
        ecss0 = V.toList ecssV

    bidstr <- readFile $ csvdr ++ "blockIDs.csv"
    let bids00 :: [Int]
        bids00 = read <$> lines bidstr
        bids0 = bids00

    chnstr <- readFile $ csvdr ++ "channels.csv"
    let chns :: [Int]
        chns = read <$> lines chnstr

    let bidspks0 = zip bids0 $ filter (\(e,_,_) -> e == 2000) ecss0
        (bidspks,(_,(_,_,spk)):_) = span (\(bid,_) -> bid /= 0) bidspks0
        bids = fst <$> bidspks
        ecss = takeWhile (\(_,_,spk') -> spk' < spk) ecss0

    let hst1,hst2,hst3 :: [[Int]]
        hst1 = sampleStreamHistogram . sampleStream $ blockStream chns bids ecss
        hst2 = rawDataHistogram $ rawDataMap chns bids ecss
        hst3 = spikeCountHistogram $ tuningMap chns bids ecss

    let hstvls = zip [(0 :: Int)..] . transpose $ foldl1 (zipWith (+)) <$> [hst1,hst2,hst3]

    let (bidspks',(_,(_,_,spk')):_) = span (\(bid,_) -> bid /= 0) $ reverse bidspks0
        bids' = fst <$> bidspks'
        ecss' = takeWhile (\(_,_,spk) -> spk' < spk) $ reverse ecss0

    print spk'
    print $ last [ spk | ((_,spk),_) <- blockStream chns bids0 ecss0 ]

    let blkstrm' = dropWhile (\((_,spk),_) -> spk < spk') $ blockStream chns bids0 ecss0

    print $ length blkstrm'

    let hst1',hst2',hst3' :: [[Int]]
        hst1' = sampleStreamHistogram $ sampleStream blkstrm'
        hst2' = rawDataHistogram $ rawDataMap chns bids' ecss'
        hst3' = spikeCountHistogram $ tuningMap chns bids' ecss'

    let hstvls' = zip [(0 :: Int)..] . transpose $ foldl1 (zipWith (+)) <$> [hst1',hst2',hst3']

    let rnbl1 = toRenderable . execEC $ do

            layout_title .= "All Spike Histogram"
            plot . fmap plotBars . liftEC $ plot_bars_values .= drop 2 hstvls

    let rnbl2 = toRenderable . execEC $ do

            layout_title .= "All Spike Histogram"
            plot . fmap plotBars . liftEC $ plot_bars_values .= drop 2 hstvls'

    let rnbl3 k = toRenderable . execEC $ do

            let hstvls'' = zip [(0 :: Int)..] . transpose $ (!! k) <$> [hst1,hst1']

            layout_title .= "Comparison"
            plot . fmap plotBars . liftEC $ plot_bars_values .= drop 2 hstvls''

    --goalRenderableToSVG ("neural-circuits/" ++ flnm) "pre-histogram" 1600 800 rnbl1
    --goalRenderableToSVG ("neural-circuits/" ++ flnm) "post-histogram" 1600 800 rnbl2
    sequence_ [ goalRenderableToSVG ("neural-circuits/" ++ flnm) ("comparison-histogram-" ++ show k) 1600 800
        $ rnbl3 k | k <- [0..54]]

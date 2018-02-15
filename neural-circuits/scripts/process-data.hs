{-# LANGUAGE DataKinds #-}

import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V

import Goal.Core
import Goal.NeuralCircuits
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
        ecss :: [(Int,Int,Double)]
        ecss = V.toList ecssV

    bidstr <- readFile $ csvdr ++ "blockIDs.csv"
    let bids :: [Int]
        bids = read <$> lines bidstr

    chnstr <- readFile $ csvdr ++ "channels.csv"
    let chns :: [Int]
        chns = read <$> lines chnstr

    let allbstrm,prebstrm,pstbstrm :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
        allbstrm = blockStream chns bids ecss
        prebstrm = takeWhile (\((bid,_),_) -> bid /= 0) allbstrm
        pstbstrm = reverse . takeWhile (\((bid,_),_) -> bid /= 0) $ reverse allbstrm

    let alltstrm,pretstrm,psttstrm :: [(BlockID,[Int])]
        alltstrm = blockToTimeStream allbstrm
        pretstrm = blockToTimeStream prebstrm
        psttstrm = blockToTimeStream pstbstrm

    let alltavg,pretavg,psttavg :: [(BlockID,[Double])]
        alltavg = timeStreamAverage alltstrm
        pretavg = timeStreamAverage pretstrm
        psttavg = timeStreamAverage psttstrm

    let ttlrnbl = toRenderable . execEC $ do

            let avgs :: [[Double]]
                avgs = map (average . snd) <$> [alltavg,pretavg,psttavg]
                ttlavgs :: [(BlockID,[Double])]
                ttlavgs = zip (fst <$> alltavg) $ transpose avgs
            plot . fmap plotBars . liftEC $ do
                plot_bars_titles .= ["all","pre","post"]
                plot_bars_values .= ttlavgs

    goalRenderableToSVG ("neural-circuits/" ++ flnm) "average-rate-histogram" 1600 800 ttlrnbl
{-
    let rnbl2 = toRenderable . execEC $ do

            layout_title .= "All Spike Histogram"
            plot . fmap plotBars . liftEC $ plot_bars_values .= drop 2 hstvls'

    let rnbl3 k = toRenderable . execEC $ do

            let hstvls'' = zip [(0 :: Int)..] . transpose $ (!! k) <$> [hst1,hst1']

            layout_title .= "Comparison"
            plot . fmap plotBars . liftEC $ plot_bars_values .= drop 2 hstvls''

    goalRenderableToSVG ("neural-circuits/" ++ flnm) "post-histogram" 1600 800 rnbl2
    sequence_ [ goalRenderableToSVG ("neural-circuits/" ++ flnm) ("comparison-histogram-" ++ show k) 1600 800
        $ rnbl3 k | k <- [0..54]]
        -}

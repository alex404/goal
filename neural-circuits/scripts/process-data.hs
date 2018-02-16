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

    adpstr <- readFile $ csvdr ++ "adaptor.csv"
    let adpts :: [Double]
        adpts = read <$> lines adpstr
        [adpt] = adpts
        adrd = 2*pi*adpt/360

    let allbstrm,prebstrm,pstbstrm :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
        allbstrm = blockStream chns bids ecss
        prebstrm = takeWhile (\((bid,_),_) -> bid /= 0) allbstrm
        pstbstrm = reverse . takeWhile (\((bid,_),_) -> bid /= 0) $ reverse allbstrm

    let alltstrm,pretstrm,psttstrm :: [(BlockID,[Int])]
        alltstrm = blockToTimeStream allbstrm
        pretstrm = blockToTimeStream prebstrm
        psttstrm = blockToTimeStream pstbstrm

    let alltavg,pretavg,psttavg :: [(BlockID,[Double])]
        alltavg = streamAverage alltstrm
        pretavg = streamAverage pretstrm
        psttavg = streamAverage psttstrm

    let allsstrm,presstrm,pstsstrm :: [(Double,[Int])]
        allsstrm = timeToSampleStream adpt alltstrm
        presstrm = timeToSampleStream adpt pretstrm
        pstsstrm = timeToSampleStream adpt psttstrm

    let allsavg,presavg,pstsavg :: [(Double,[Double])]
        allsavg = streamAverage allsstrm
        presavg = streamAverage presstrm
        pstsavg = streamAverage pstsstrm

    let bids' = fst <$> alltavg

    let tttlrnbl = toRenderable . execEC $ do

            let ttlavgs = zip bids' . transpose $ map (average . snd) <$> [pretavg,alltavg,psttavg]

            plot . fmap plotBars . liftEC $ do
                plot_bars_titles .= ["pre","all","post"]
                plot_bars_values .= ttlavgs
                plot_bars_alignment .= BarsLeft

    let tnrnrnbl k = toRenderable . execEC $ do

            let ttls = zip bids' . transpose $ map ((!! k) . snd) <$> [pretavg,alltavg,psttavg]
            plot . fmap plotBars . liftEC $ do
                plot_bars_titles .= ["pre","all","post"]
                plot_bars_values .= ttls
                plot_bars_alignment .= BarsLeft

    let sttlrnbl = toRenderable . execEC $ do

            plot . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [[ (stm, sum spks) | (stm,spks) <- presavg ]]
                plot_lines_style .= solidLine 4 (opaque blue)

            plot . liftEC $ do
                plot_lines_title .= "all"
                plot_lines_values .= [[ (stm, sum spks) | (stm,spks) <- allsavg ]]
                plot_lines_style .= solidLine 4 (opaque red)

            plot . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [[ (stm, sum spks) | (stm,spks) <- pstsavg ]]
                plot_lines_style .= solidLine 4 (opaque green)

            plot . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adrd

    let snrnrnbl k = toRenderable . execEC $ do

            plot . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [[ (stm, spks !! k) | (stm,spks) <- presavg ]]
                plot_lines_style .= solidLine 4 (opaque blue)

            plot . liftEC $ do
                plot_lines_title .= "all"
                plot_lines_values .= [[ (stm, spks !! k) | (stm,spks) <- allsavg ]]
                plot_lines_style .= solidLine 4 (opaque red)

            plot . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [[ (stm, spks !! k) | (stm,spks) <- pstsavg ]]
                plot_lines_style .= solidLine 4 (opaque green)

            plot . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adrd


    goalRenderableToSVG ("neural-circuits/" ++ flnm) "average-rate-histogram" 1600 800 tttlrnbl
    goalRenderableToSVG ("neural-circuits/" ++ flnm) "average-rate-lines" 1600 800 sttlrnbl

    sequence_ $ do
        k <- [0..54]
        return . goalRenderableToSVG
            ("neural-circuits/" ++ flnm) ("comparison-histogram-" ++ show k) 1600 800 $ tnrnrnbl k

    sequence_ $ do
        k <- [0..54]
        return . goalRenderableToSVG
            ("neural-circuits/" ++ flnm) ("comparison-lines-" ++ show k) 1600 800 $ snrnrnbl k

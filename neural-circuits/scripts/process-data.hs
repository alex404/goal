{-# LANGUAGE DataKinds #-}

import Goal.Core
import Goal.NeuralCircuits
import Data.List


--- Globals ---

dr,flnm :: String
dr = "adaptation/small40"
flnm = "112l44/"


--- Main ---


main :: IO ()
main = do

    ecss <- getSpikes dr flnm
    bids <- getBIDs dr flnm
    chns <- getChannels dr flnm
    adpt <- getAdaptor dr flnm

    let adrd = 2*pi*adpt/360

    let allbstrm,prebstrm,pstbstrm :: [((BlockID,SpikeTime),[(NeuronID,SpikeTime)])]
        allbstrm = blockStream chns bids ecss
        prebstrm = takeWhile (\((bid,_),_) -> bid /= 0) allbstrm
        pstbstrm = takeWhile (\((bid,_),_) -> bid == 0) $ dropWhile (\((bid,_),_) -> bid /= 0) allbstrm

    print $ length pstbstrm

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

    let nnrns = length . snd $ head allsavg

    let tttlrnbl = toRenderable . execEC $ do

            goalLayout

            let ttlgrps = groupBy (\(x,_) (y,_) -> x == y) . sort $ concat [pretavg,alltavg,psttavg]
                bids' = fst . head <$> ttlgrps
                ttlavgs = zip bids' $ map (sum . snd) <$> ttlgrps

            plot . fmap plotBars . liftEC $ do
                plot_bars_titles .= ["pre","all","post"]
                plot_bars_values .= ttlavgs
                plot_bars_alignment .= BarsLeft

    let tnrnrnbl k = toRenderable . execEC $ do

            goalLayout

            let ttlgrps = groupBy (\(x,_) (y,_) -> x == y) . sort $ concat [pretavg,alltavg,psttavg]
                bids' = fst . head <$> ttlgrps
                ttls = zip bids' $ map ((!! k) . snd) <$> ttlgrps

            plot . fmap plotBars . liftEC $ do
                plot_bars_titles .= ["pre","all","post"]
                plot_bars_values .= ttls
                plot_bars_alignment .= BarsLeft

    let sttlrnbl = toRenderable . execEC $ do

            goalLayout
            radiansAbscissa

            plot . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, sum spks) | (stm,spks) <- presavg ] ]
                plot_lines_style .= solidLine 4 (opaque blue)

            plot . liftEC $ do
                plot_lines_title .= "all"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, sum spks) | (stm,spks) <- allsavg ] ]
                plot_lines_style .= solidLine 4 (opaque red)

            plot . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, sum spks) | (stm,spks) <- pstsavg ] ]
                plot_lines_style .= solidLine 4 (opaque green)

            plot . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adrd

    let snrnrnbl k = toRenderable . execEC $ do

            goalLayout
            radiansAbscissa
            let mxy = maximum $ (!! k) . snd <$> (presavg ++ pstsavg)
            layout_y_axis . laxis_generate .= scaledAxis def (0,1.5*mxy)

            plot . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, spks !! k) | (stm,spks) <- presavg ] ]
                plot_lines_style .= solidLine 4 (opaque blue)
{-

            plot . liftEC $ do
                plot_lines_title .= "all"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, spks !! k) | (stm,spks) <- allsavg ] ]
                plot_lines_style .= solidLine 4 (opaque red)
                -}

            plot . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [ loopRadiansPlotData [ (stm, spks !! k) | (stm,spks) <- pstsavg ] ]
                plot_lines_style .= solidLine 4 (opaque green)

            plot . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adrd

    goalRenderableToSVG ("neural-circuits/" ++ flnm) "average-rate-histogram" 1600 800 tttlrnbl
    goalRenderableToSVG ("neural-circuits/" ++ flnm) "average-rate-lines" 1600 800 sttlrnbl

    sequence_ $ do
        k <- take nnrns [0..]
        return . goalRenderableToSVG
            ("neural-circuits/" ++ flnm) ("comparison-histogram-" ++ show k) 1600 800 $ tnrnrnbl k

    sequence_ $ do
        k <- take nnrns [0..]
        return . goalRenderableToSVG
            ("neural-circuits/" ++ flnm) ("comparison-lines-" ++ show k) 1600 800 $ snrnrnbl k

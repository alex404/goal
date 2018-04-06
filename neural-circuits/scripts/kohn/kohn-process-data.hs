{-# LANGUAGE DataKinds #-}

import Goal.Core
import Goal.NeuralCircuits

import qualified Data.Map as M

--- Globals ---

dr,flnm,sbdr :: String
dr = "adaptation/small40"
flnm = "112l45/"
sbdr = "neural-circuits/kohn-data/" ++ flnm


--- Main ---


main :: IO ()
main = do

    ecss <- getSpikes dr flnm
    bids <- getBIDs dr flnm
    chns <- getChannels dr flnm
    adpt <- getAdaptor dr flnm

    let adrd = (+2) . round $ adpt/22.5
        adrd' = if adrd >= 10 then adrd - 8 else adrd + 8
        adptpi0 = 2*pi*adpt/360
        adptpi =
            2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

    let allbstrm,prebstrm,pstbstrm :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        allbstrm = take 1531 $ blockStream chns bids ecss
        prebstrm = take 850 allbstrm
        pstbstrm = drop 851 allbstrm

    let prebids,pstbids :: M.Map BlockID (M.Map NeuronID [SpikeTime])
        prebids = averageBlockIDs bids prebstrm
        pstbids = averageBlockIDs bids pstbstrm

    let prestms,pststms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])
        prestms = averageBlockIDsToStimuli prebids
        pststms = averageBlockIDsToStimuli pstbids

    let bidrnbl preln pstln = toRenderable . execEC $ do

            goalLayoutLR

            let dffln = zipWith (\(x1,y1) (_,y2) -> (x1,y1 - y2)) preln pstln

            plotRight . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adrd

            plotRight . return $ vlinePlot "adaptor" (solidLine 4 $ withOpacity black 0.5) adrd'

            plotRight . liftEC $ do
                plot_lines_values .= [[(0,0),(1,0)]]
                plot_lines_style .= solidLine 4 (opaque white)

            plotLeft . liftEC $ do
                plot_lines_values .= [[(0,0),(1,0)]]
                plot_lines_style .= solidLine 4 (opaque white)

            plotRight . liftEC $ do
                plot_lines_title .= "diff"
                plot_lines_values .= [dffln]
                plot_lines_style .= dashedLine 4 [4,4] (opaque blue)

            plotLeft . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [preln]
                plot_lines_style .= solidLine 4 (opaque black)

            plotLeft . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [pstln]
                plot_lines_style .= solidLine 4 (opaque red)

    let stmrnbl preln pstln = toRenderable . execEC $ do

            goalLayoutLR
            radiansAbscissaLR

            let dffln = zipWith (\(x1,y1) (_,y2) -> (x1,y1 - y2)) preln pstln

            plotRight . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adptpi

            plotRight . liftEC $ do
                plot_lines_values .= [[(0,0),(1,0)]]
                plot_lines_style .= solidLine 4 (opaque white)

            plotLeft . liftEC $ do
                plot_lines_values .= [[(0,0),(1,0)]]
                plot_lines_style .= solidLine 4 (opaque white)

            plotRight . liftEC $ do
                plot_lines_title .= "diff"
                plot_lines_values .= [dffln]
                plot_lines_style .= dashedLine 4 [4,4] (opaque blue)

            plotLeft . liftEC $ do
                plot_lines_title .= "pre"
                plot_lines_values .= [preln]
                plot_lines_style .= solidLine 4 (opaque black)

            plotLeft . liftEC $ do
                plot_lines_title .= "post"
                plot_lines_values .= [pstln]
                plot_lines_style .= solidLine 4 (opaque red)

    let nrns = M.keys . snd . head $ allbstrm
        shave = reverse . drop 1 . reverse . drop 2

        bidttlpreln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList prebids ]
        bidttlpstln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList pstbids ]
        bidnrnpreln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList prebids ]
        bidnrnpstln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList pstbids ]

        bidttlrnbl = bidrnbl bidttlpreln bidttlpstln
        bidnrnrnbl nrn = bidrnbl (bidnrnpreln nrn) (bidnrnpstln nrn)

        stmttlpreln = [ (stm, sum . M.elems $ length <$> nmp) | (stm,nmp) <- M.toList prestms ]
        stmttlpstln = [ (stm, sum . M.elems $ length <$> nmp) | (stm,nmp) <- M.toList pststms ]
        stmnrnpreln nrn = [ (stm, length $ nmp M.! nrn) | (stm,nmp) <- M.toList prestms ]
        stmnrnpstln nrn = [ (stm, length $ nmp M.! nrn) | (stm,nmp) <- M.toList pststms ]

        stmttlrnbl = stmrnbl stmttlpreln stmttlpstln
        stmnrnrnbl nrn = stmrnbl (stmnrnpreln nrn) (stmnrnpstln nrn)

    putStr "Number of Neurons: "
    print $ length nrns

    goalRenderableToSVG sbdr "bid-total-spikes" 1200 800 bidttlrnbl
    goalRenderableToSVG sbdr "stimulus-total-spikes" 1200 800 stmttlrnbl

    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            sbdr ("bid-spikes-neuron-" ++ show nrn) 1200 800 $ bidnrnrnbl nrn

    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            sbdr ("stimulus-spikes-neuron-" ++ show nrn) 1200 800 $ stmnrnrnbl nrn

    goalWriteFile sbdr "prestms" $ show prestms
    goalWriteFile sbdr "pststms" $ show pststms

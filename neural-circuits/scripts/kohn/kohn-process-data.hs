{-# LANGUAGE DataKinds #-}

import Goal.Core
import Goal.NeuralCircuits

import Data.List
import qualified Data.Map as M


--- Globals ---


--- Main ---


processData :: (KnownNat nn, KnownNat t1, KnownNat t2) => KohnExperiment nn t1 t2 -> IO ()
processData kxp = do

    putStrLn $ "\nPROTOCOL: " ++ protocol kxp ++ " | EXPERIMENT: " ++ experiment kxp ++ "\n"

    let sbdr = kohnProjectPath kxp

    ecss <- getSpikes kxp
    bids <- getBIDs kxp
    chns <- getChannels kxp
    adpt <- getAdaptor kxp

    let bstrm :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        bstrm = blockStream chns bids ecss

    let prtclttl,prtcln1 :: Int
        prtclttl = length $ dropWhile ((/= 0) .  fst . fst) $ reverse bstrm
        prtcln1 = length $ takeWhile ((/= 0) .  fst . fst) bstrm

    let allbstrm,prebstrm,pstbstrm :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        allbstrm = take prtclttl bstrm
        prebstrm = take prtcln1 allbstrm
        pstbstrm = drop (prtcln1 + 1) allbstrm

    let prebids,pstbids :: M.Map BlockID (M.Map NeuronID [SpikeTime])
        prebids = averageBlockIDs bids prebstrm
        pstbids = averageBlockIDs bids pstbstrm

    let prestms,pststms :: M.Map Stimulus (M.Map NeuronID [SpikeTime])
        prestms = averageBlockIDsToStimuli prebids
        pststms = averageBlockIDsToStimuli pstbids

    let nrns = M.keys . snd . head $ allbstrm
        shave = reverse . drop 1 . reverse . drop 2

        bidttlpreln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList prebids ]
        bidttlpstln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList pstbids ]
        bidnrnpreln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList prebids ]
        bidnrnpstln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList pstbids ]

        bidttlrnbl = toRenderable $ blockIDTuningCurves adpt bidttlpreln bidttlpstln
        bidnrnrnbl nrn = toRenderable $ blockIDTuningCurves adpt (bidnrnpreln nrn) (bidnrnpstln nrn)

        stmttlpreln = [ (stm, sum . M.elems $ length <$> nmp) | (stm,nmp) <- M.toList prestms ]
        stmttlpstln = [ (stm, sum . M.elems $ length <$> nmp) | (stm,nmp) <- M.toList pststms ]
        stmnrnpreln nrn = [ (stm, length $ nmp M.! nrn) | (stm,nmp) <- M.toList prestms ]
        stmnrnpstln nrn = [ (stm, length $ nmp M.! nrn) | (stm,nmp) <- M.toList pststms ]

        stmttlrnbl = toRenderable $ stimulusTuningCurves adpt stmttlpreln stmttlpstln
        stmnrnrnbl nrn = toRenderable $ stimulusTuningCurves adpt (stmnrnpreln nrn) (stmnrnpstln nrn)

    let prestrm = blockToStimulusStream prebstrm
        pststrm = blockToStimulusStream pstbstrm

    putStr "Number of Neurons: "
    print $ length nrns
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length prestrm
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length pststrm
    putStrLn "Block ID Trial Counts: "
    print . map length . group . sort $ fst . fst <$> allbstrm

    goalRenderableToSVG sbdr "bid-total-spikes" 1200 800 bidttlrnbl
    goalRenderableToSVG sbdr "stimulus-total-spikes" 1200 800 stmttlrnbl

    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            (sbdr ++ "/individual") ("bid-spikes-neuron-" ++ show nrn) 1200 800 $ bidnrnrnbl nrn

    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            (sbdr ++ "/individual") ("stimulus-spikes-neuron-" ++ show nrn) 1200 800 $ stmnrnrnbl nrn

    goalWriteFile sbdr "prestms" $ show prestms
    goalWriteFile sbdr "pststms" $ show pststms
    goalWriteFile sbdr "prestrm" $ show prestrm
    goalWriteFile sbdr "pststrm" $ show pststrm

main :: IO ()
main = do

    processData experiment112l44
    processData experiment112l45
    processData experiment112r35
    processData experiment112r36
    processData experiment105r62
    processData experiment107l114
    processData experiment112l16
    processData experiment112r29
    processData experiment112r32
    putStrLn "\n"

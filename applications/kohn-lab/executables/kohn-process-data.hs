{-# LANGUAGE DataKinds #-}

import KohnLab

import Goal.Core

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

    let bidstrm,bidstrm0,bidstrm1 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]
        bidstrm = take prtclttl bstrm
        bidstrm0 = take prtcln1 bidstrm
        bidstrm1 = drop (prtcln1 + 1) bidstrm

    let stmstrm0,stmstrm1 :: [(Stimulus, M.Map NeuronID [SpikeTime])]
        stmstrm0 = blockToStimulusStream bidstrm0
        stmstrm1 = blockToStimulusStream bidstrm1

    let bidttls0,bidttls1 :: M.Map BlockID (M.Map NeuronID [SpikeTime])
        bidttls0 = blockIDTotals bids bidstrm0
        bidttls1 = blockIDTotals bids bidstrm1

    let stmttls0,stmttls1 :: M.Map Stimulus (M.Map NeuronID [SpikeTime])
        stmttls0 = blockIDToStimulusTotals bidttls0
        stmttls1 = blockIDToStimulusTotals bidttls1

    let nrns = M.keys . snd . head $ bidstrm

    putStr "Number of Neurons: "
    print $ length nrns
    putStr "Adaptor: "
    print adpt
    putStr "Number of Filtered Pre-Adaptation Trials: "
    print $ length stmstrm0
    putStr "Number of Filtered Post-Adaptation Trials: "
    print $ length stmstrm1
    putStrLn "Block ID Trial Counts: "
    print . map length . group . sort $ fst . fst <$> bidstrm

    goalWriteFile sbdr "bidstrm0" $ show bidstrm0
    goalWriteFile sbdr "bidstrm1" $ show bidstrm1
    goalWriteFile sbdr "bidttls0" $ show bidttls0
    goalWriteFile sbdr "bidttls1" $ show bidttls1
    goalWriteFile sbdr "stmstrm0" $ show stmstrm0
    goalWriteFile sbdr "stmstrm1" $ show stmstrm1
    goalWriteFile sbdr "stmttls0" $ show stmttls0
    goalWriteFile sbdr "stmttls1" $ show stmttls1

main :: IO ()
main = do

    processData experiment112l44
    processData experiment112l45
    processData experiment112r35
    processData experiment112r36
    processData experiment105r62
    processData experiment107l114
    processData experiment112l16
    processData experiment112r32
    putStrLn "\n"

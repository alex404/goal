{-# LANGUAGE ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core

import qualified Data.Map as M


--- Globals ---


--- Main ---


analyzeBIDData :: (KnownNat nn, KnownNat t1, KnownNat t2) => KohnExperiment nn t1 t2 -> IO ()
analyzeBIDData kxp = do

    let kpp = kohnProjectPath kxp
        biddr = kpp ++ "/bid-analysis"

    adpt <- getAdaptor kxp

    --(bidstrm0 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpp "bidstrm0"
    --(bidstrm1 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpp "bidstrm1"
    (bidttls0 :: M.Map BlockID (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpp "bidttls0"
    (bidttls1 :: M.Map BlockID (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpp "bidttls1"

    let nrns = M.keys . (!! 2) $ M.elems bidttls0
        shave = reverse . drop 1 . reverse . drop 2

        bidttlpreln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList bidttls0 ]
        bidttlpstln = shave [ (bid, sum . M.elems $ length <$> nmp) | (bid,nmp) <- M.toList bidttls1 ]
        bidnrnpreln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList bidttls0 ]
        bidnrnpstln nrn = shave [ (bid, length $ nmp M.! nrn) | (bid,nmp) <- M.toList bidttls1 ]

        bidttlrnbl = toRenderable $ blockIDTuningCurves adpt bidttlpreln bidttlpstln
        bidnrnrnbl nrn = toRenderable $ blockIDTuningCurves adpt (bidnrnpreln nrn) (bidnrnpstln nrn)

    goalRenderableToSVG biddr "bid-total-spikes" 1200 800 bidttlrnbl

    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            (biddr ++ "/individual") ("bid-spikes-neuron-" ++ show nrn) 1200 800 $ bidnrnrnbl nrn

main :: IO ()
main = do

    analyzeBIDData experiment112l44
    analyzeBIDData experiment112l45
    analyzeBIDData experiment112r35
    analyzeBIDData experiment112r36
    analyzeBIDData experiment105r62
    analyzeBIDData experiment107l114
    analyzeBIDData experiment112l16
    analyzeBIDData experiment112r32
    putStrLn "\n"

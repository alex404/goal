{-# LANGUAGE ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core

import qualified Data.Map as M
import Data.List


--- Globals ---


--- Plots ---


shave :: [a] -> [a]
shave = reverse . drop 1 . reverse . drop 2

blockIDLiner
    :: Maybe NeuronID
    -> M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
    -> [(BlockID,Double)]
blockIDLiner (Just nrn) bidttls = shave $ do
    (bid,(k,nmp)) <- M.toList bidttls
    return (bid, (/ fromIntegral k) . genericLength $ nmp M.! nrn)
blockIDLiner Nothing bidttls = shave $ do
    (bid,(k,nmp)) <- M.toList bidttls
    return (bid, (/ fromIntegral k) . sum . M.elems $ genericLength <$> nmp)

blockIDTuningCurves
    :: Stimulus
    -> M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
    -> M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
    -> Maybe NeuronID
    -> LayoutLR Int Double Double
blockIDTuningCurves adpt bidttls0 bidttls1 mnrn = execEC $ do

    goalLayoutLR

    let preln = blockIDLiner mnrn bidttls0
        pstln = blockIDLiner mnrn bidttls1

    let adrd = (+2) . round $ adpt/22.5
        adrd' = if adrd >= 10 then adrd - 8 else adrd + 8

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



--- Main ---


analyzeBIDData :: (KnownNat nn, KnownNat t1, KnownNat t2) => KohnExperiment nn t1 t2 -> IO ()
analyzeBIDData kxp = do

    let kpdr = kohnProjectPath kxp

    adpt <- getAdaptor kxp

    --(bidstrm0 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "bidstrm0"
    --(bidstrm1 :: [(BlockEvent, M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "bidstrm1"
    (bidttls0 :: M.Map BlockID (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "bidttls0"
    (bidttls1 :: M.Map BlockID (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "bidttls1"

    let nrns = M.keys . snd . (!! 2) $ M.elems bidttls0

    renderNeuralLayouts kpdr "bid-spike-average" nrns $ blockIDTuningCurves adpt bidttls0 bidttls1

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

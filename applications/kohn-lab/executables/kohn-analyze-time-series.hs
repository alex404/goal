{-# LANGUAGE ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core

import qualified Data.Map as M


--- Globals ---


--- Plots ---

spikeTimeSeriesLayout
    :: Int
    -> [(Stimulus, M.Map NeuronID [SpikeTime])]
    -> Maybe NeuronID
    -> Layout Int Int
spikeTimeSeriesLayout adptn stmstrm mnrn = execEC $ do

    goalLayout
    layout_x_axis . laxis_title .= "Trial"
    layout_y_axis . laxis_title .= "Number of Spikes"

    let spks = snd <$> streamToSpikeCounts mnrn stmstrm

    plot . liftEC $ do

        plot_lines_style .= solidLine 2 (opaque blue)
        plot_lines_values .= [ zip [(0 :: Int)..] spks ]

    plot . return $ vlinePlot "" (solidLine 4 $ opaque black) adptn


--- Main ---


analyzeTimeSeries :: (KnownNat nn, KnownNat t1, KnownNat t2) => KohnExperiment nn t1 t2 -> IO ()
analyzeTimeSeries kxp = do

    let kpdr = kohnProjectPath kxp

    (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm0"
    (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm1"

    let nrns = M.keys . snd . head $ stmstrm0
        adptn = length stmstrm0
        stmstrm = stmstrm0 ++ stmstrm1

    renderNeuralLayouts kpdr "spike-count-time-series" nrns $ spikeTimeSeriesLayout adptn stmstrm

main :: IO ()
main = do

    analyzeTimeSeries experiment112l44
    analyzeTimeSeries experiment112l45
    analyzeTimeSeries experiment112r35
    analyzeTimeSeries experiment112r36
    analyzeTimeSeries experiment105r62
    analyzeTimeSeries experiment107l114
    analyzeTimeSeries experiment112l16
    analyzeTimeSeries experiment112r32

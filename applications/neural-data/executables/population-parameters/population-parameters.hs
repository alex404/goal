{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData
import NeuralData.VonMises

import Paths_neural_data

import Goal.Core
import Goal.Probability


--- Globals ---


ananm :: String
ananm = "population-parameters"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000


--- CLI ---


data AnalysisOpts = AnalysisOpts String String

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument ( help "Which data collection to analyze" )
    <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg) = do

    dsts <- if null dstarg
               then fromJust <$> goalReadDatasetsCSV (Experiment prjnm expnm)
               else return [dstarg]

    tcgpi <- getDataFileName "population-parameters/tuning-curves.gpi"
    ppgpi <- getDataFileName "population-parameters/population-parameter-histogram.gpi"

    let expmnt = Experiment prjnm expnm


    forM_ dsts $ \dst -> do

        let msbexpt = Just $ SubExperiment "tuning-curves" dst
            msbexph = Just $ SubExperiment "histograms" dst

        (k,(zxs0 :: [([Int], Double)])) <- getNeuralData expnm dst

        (tcss,hstcsv:hstcsvs) <- realize $ case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do
                let zxs :: [(Response k, Double)]
                    zxs = strengthenNeuralData zxs0
                lkl <- fitIPLikelihood zxs
                return (analyzeTuningCurves xsmps lkl,populationParameterHistogram 10 lkl)

        goalWriteAnalysis True expmnt msbexpt tcss

        runGnuplot expmnt msbexpt defaultGnuplotOptions tcgpi

        goalWriteNamedAnalysis True expmnt msbexph hstcsv
        mapM_ (goalWriteNamedAnalysis False expmnt msbexph) hstcsvs

        runGnuplot expmnt msbexph defaultGnuplotOptions ppgpi



--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts

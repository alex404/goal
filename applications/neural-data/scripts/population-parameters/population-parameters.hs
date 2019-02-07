{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData
import NeuralData.VonMises

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Globals ---


ananm :: String
ananm = "population-parameters"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000


--- CLI ---


runOpts :: ExperimentOpts -> IO ()
runOpts expopts@(ExperimentOpts expnm _) = do

    dsts <- readDatasets expopts

    let tcgpi = "tuning-curves.gpi"
        ppgpi = "population-parameter-histogram.gpi"

    let expmnt = Experiment prjnm expnm

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let msbexpt = Just $ SubExperiment "tuning-curves" dst
            msbexph = Just $ SubExperiment "histograms" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        (tcss,(hstcsv:hstcsvs,ppds),(sgns,smus,sprcs)) <- realize $ case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do
                let zxs :: [(Response k, Double)]
                    zxs = strengthenNeuralData zxs0
                lkl <- fitIPLikelihood zxs
                let (hstcsvppds,prmss) = populationParameters 20 lkl
                return (analyzeTuningCurves xsmps lkl,unzip hstcsvppds,prmss)

        putStrLn "Gains Log-Normal Parameters:"
        print $ listCoordinates sgns
        putStrLn "Preferred Stimuli Von Mises Parameters:"
        print $ listCoordinates smus
        putStrLn "Precisions Log-Normal Parameters:"
        print $ listCoordinates sprcs

        goalExport True expmnt msbexpt tcss

        runGnuplot expmnt msbexpt defaultGnuplotOptions tcgpi

        goalExportNamed True expmnt msbexph hstcsv
        mapM_ (goalExportNamed False expmnt msbexph) hstcsvs
        mapM_ (goalExportNamed False expmnt msbexph) ppds

        runGnuplot expmnt msbexph defaultGnuplotOptions ppgpi



--- Main ---


main :: IO ()
main = do

    let prgstr =
            "Analyze the basic parameters of an independent Poisson likelihood \
            \model. Produces csvs and plots of the tuning-curves of the model, as \
            \well as histograms of the parameters."
        hdrstr = "Analyze and plot statistics of the parameters of a given neural population."
        opts = info (experimentOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

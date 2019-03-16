{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    ScopedTypeVariables,
    DataKinds
    #-}


--- Imports ---


import NeuralData
import NeuralData.Mixture

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---


nbns :: Int
nbns = 10

xsmps :: [Double]
xsmps = init $ range mnx mxx 101


--- CLI ---


data ValidationOpts = ValidationOpts Int NatNumber Double Double Double Int Int Double Double

validationOpts :: Parser ValidationOpts
validationOpts = ValidationOpts
    <$> option auto
        ( short 'n'
        <> long "n-population"
        <> help "Number of shotgun populations to generate."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 'm'
        <> long "n-components"
        <> help "Number of components."
        <> showDefault
        <> value 7 )
    <*> option auto
        ( short 'M'
        <> long "concentration"
        <> help "Concetration of mixture weights."
        <> showDefault
        <> value 2 )
    <*> option auto
        ( short 'l'
        <> long "learning-rate"
        <> help "The learning rate."
        <> showDefault
        <> value (-0.05) )
    <*> option auto
        ( short 'w'
        <> long "weight-decay"
        <> help "Weight decay rate."
        <> showDefault
        <> value 0.001 )
    <*> option auto
        ( short 'b'
        <> long "n-batch"
        <> help "Batch size."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 'e'
        <> long "n-epochs"
        <> help "Number of batches to run the learning over."
        <> showDefault
        <> value 5000 )
    <*> option auto
        ( short 'p'
        <> long "log-mu-precision"
        <> help "The mu parameter of the initial precision log-normal."
        <> showDefault
        <> value (-1) )
    <*> option auto
        ( short 'P'
        <> long "log-sd-precision"
        <> help "The sd parameter of the initial precision log-normal."
        <> showDefault
        <> value 0.5 )

data AllOpts = AllOpts ExperimentOpts ValidationOpts

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> validationOpts

runOpts :: AllOpts -> IO ()
runOpts ( AllOpts expopts@(ExperimentOpts expnm _)
    (ValidationOpts npop nmx cnc eps dcy nbtch nepchs pmu psd) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

        lgnrm :: Natural # LogNormal
        lgnrm = toNatural . Point @ Source $ S.doubleton pmu psd

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let shtanl = Just $ Analysis "shotgun-training" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        putStrLn "\nNumber of Samples:"
        print $ length zxs0

        case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do

                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs2 <- realize $ shuffleList zxs1

                let zxs = zxs2

                case someNatVal (nmx-1)
                    of SomeNat (Proxy :: Proxy m) -> do

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                        (sgdnrms, nnans, mlkls) <- realize
                            $ shotgunFitMixtureLikelihood npop eps 0 nbtch nepchs drch lgnrm zxs

                        (dcysgdnrms, dcynnans, dcymlkls) <- realize
                            $ shotgunFitMixtureLikelihood npop eps dcy nbtch nepchs drch lgnrm zxs

                        let mlkl = last mlkls
                            dcymlkl = last dcymlkls

                        putStrLn $ concat ["\nNumber of NaNs: ", show nnans , " / ", show npop]
                        putStrLn $ concat ["\nNumber of NaNs (Decayed): ", show dcynnans , " / ", show npop]

                        goalExportNamed True expmnt shtanl sgdnrms
                        goalExportNamed False expmnt shtanl dcysgdnrms

                        let rltv = "../population-parameters/"
                            ttl = "shotgun-training"
                            dcyttl = "decayed-shotgun-training"

                        runPopulationParameterAnalyses expmnt dst xsmps nbns rltv ttl Nothing Nothing mlkl
                        runPopulationParameterAnalyses
                            expmnt dst xsmps nbns rltv dcyttl Nothing Nothing dcymlkl

                runGnuplot expmnt shtanl defaultGnuplotOptions "cross-entropy-descent.gpi"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

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
import NeuralData.VonMises
import NeuralData.Mixture

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L


--- Globals ---



--- CLI ---


data ValidationOpts = ValidationOpts Int Double [Double] Double Int Int Double Double

validationOpts :: Parser ValidationOpts
validationOpts = ValidationOpts
    <$> option auto
        ( short 'n'
        <> long "n-samples"
        <> help "Number of sample populations to generate."
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 'v'
        <> long "validation-percent"
        <> help "Percent of samples to withhold for validation."
        <> showDefault
        <> value 0.2 )
    <*> many (option auto
        ( short 'm'
        <> long "dirichlet"
        <> help "Dirichlet parameters (and consequently number of mixers)") )
    <*> option auto
        ( short 'l'
        <> long "learning-rate"
        <> help "The learning rate."
        <> showDefault
        <> value (-0.05) )
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
        <> value 500 )
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
    (ValidationOpts npop vprcnt mxs eps nbtch nepchs pmu psd) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

        lgnrm :: Natural # LogNormal
        lgnrm = toNatural . Point @ Source $ S.doubleton pmu psd

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let ceexp = Just $ SubExperiment "cross-validation" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        (sgdnrms, nnans, _) <- realize $ case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do
                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs <- shuffleList zxs1
                let (tzxs,vzxs) = splitAt (round . (*vprcnt) . fromIntegral $ length zxs) zxs
                case someNatVal (L.genericLength mxs - 1)
                    of SomeNat (Proxy :: Proxy m) -> do
                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point . fromJust $ S.fromList mxs
                        multiFitMixtureLikelihood npop eps nbtch nepchs drch lgnrm tzxs vzxs

        goalExportNamed True expmnt ceexp sgdnrms

        putStrLn $ concat ["Number of NaNs: ", show nnans , " / ", show npop]

        runGnuplot expmnt ceexp defaultGnuplotOptions "cross-entropy-descent.gpi"



--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

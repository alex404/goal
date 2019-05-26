{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
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


data ValidationOpts = ValidationOpts Int NatNumber Double Double Double Int Int

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
        <> help "Concentration of mixture weights."
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

data AllOpts = AllOpts ExperimentOpts ValidationOpts

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> validationOpts

runOpts :: AllOpts -> IO ()
runOpts ( AllOpts expopts@(ExperimentOpts expnm _)
    (ValidationOpts npop nmx cnc eps _ nbtch nepchs) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        putStrLn "\nNumber of Samples:"
        print $ length zxs0

        case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do

                let zs1 :: [Response k]
                    zs1 = fst <$> strengthenNeuralData zxs0
                zs <- realize $ shuffleList zs1

                --let skp = round $ fromIntegral (length zs) / (fromIntegral nbtch :: Double)

                case someNatVal (nmx-1)
                    of SomeNat (Proxy :: Proxy m) -> do

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                        (sgdnrms, nnans, mmdl) <- realize
                            $ shotgunFitMixture npop eps 0 nbtch nepchs drch zs

                        --(dcysgdnrms, dcynnans, dcymmdls) <- realize
                        --    $ shotgunFitMixture skp npop eps dcy nbtch nepchs drch zs

                        (emsgdnrms, emnnans, emmmdl) <- realize
                            $ shotgunEMFitMixture npop 100 drch zs

                        putStrLn $ concat ["\nNumber of NaNs: ", show nnans , " / ", show npop]
                        putStrLn $ concat ["\nNumber of NaNs (EM): ", show emnnans , " / ", show npop]

                        let stmianl = Just $ Analysis "stimulus-independent" dst

                        goalExportNamed True expmnt stmianl sgdnrms
                        goalExportNamed False expmnt stmianl emsgdnrms

                        runGnuplot expmnt stmianl defaultGnuplotOptions "cross-entropy-descent.gpi"

                        let mdlcrs = mixturePopulationNoiseCorrelations mmdl
                            emmdlcrs = mixturePopulationNoiseCorrelations emmmdl
                            mvncrs = estimateCorrelations zs
                            mlticrs = S.combineTriangles (S.replicate 1) mdlcrs mvncrs
                            emmlticrs = S.combineTriangles (S.replicate 1) emmdlcrs mvncrs

                            rws = S.toList <$> S.toList (S.toRows mlticrs)
                            emrws = S.toList <$> S.toList (S.toRows emmlticrs)

                        goalExport False expmnt stmianl rws
                        goalExport False expmnt stmianl emrws

                        runGnuplot expmnt stmianl defaultGnuplotOptions "noise-correlations.gpi"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

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

nstms :: Int
nstms = 8

stms :: [Double]
stms = tail $ range mnx mxx (nstms + 1)

xsmps :: [Double]
xsmps = init $ range mnx mxx 101


--- CLI ---


data ValidationOpts = ValidationOpts Int Int Int Double NatNumber Double Int Int Double Double

validationOpts :: Parser ValidationOpts
validationOpts = ValidationOpts
    <$> option auto
        ( short 'n'
        <> long "n-population"
        <> help "Number of sample populations to generate."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 'f'
        <> long "k-fold-validation"
        <> help "Number of (k-)folds."
        <> showDefault
        <> value 5 )
    <*> option auto
        ( short 'm'
        <> long "dirichlet"
        <> help "Number of mixture model counts to test."
        <> showDefault
        <> value 8 )
    <*> option auto
        ( short 'M'
        <> long "concentration"
        <> help "Concetration of mixture weights."
        <> showDefault
        <> value 2 )
    <*> option auto
        ( short 's'
        <> long "mixture-step"
        <> help "Number of mixture counts to step each iteration."
        <> showDefault
        <> value 1 )
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
    (ValidationOpts npop kfld nmx cnc nstp eps nbtch nepchs pmu psd) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

        lgnrm :: Natural # LogNormal
        lgnrm = toNatural . Point @ Source $ S.doubleton pmu psd

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let cesbexp = Just $ SubExperiment "cross-validation" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        putStrLn "\nNumber of Mixers:"

        case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do

                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs <- realize $ shuffleList zxs1

                let idxs = take nmx [0,nstp..]

                cvls <- forM idxs $ \m -> case someNatVal m

                    of SomeNat (Proxy :: Proxy m) -> do

                        print $ m+1

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                        (sgdnrms, nnans, mlkl) <- realize
                            $ shotgunFitMixtureLikelihood npop eps nbtch nepchs drch lgnrm zxs

                        goalExportNamed (m==0) expmnt cesbexp sgdnrms

                        putStrLn $ concat ["\nNumber of NaNs: ", show nnans , " / ", show npop]

                        let hstnm = "population-parameters-" ++ show (m+1) ++ "-mixers"

                            hstsbexp = Just $ SubExperiment hstnm dst
                            (tcs,gps,ccrvs,ctcrvs) = analyzePopulationCurves xsmps mlkl

                        goalExportNamed True expmnt hstsbexp tcs
                        goalExportNamed False expmnt hstsbexp gps
                        goalExportNamed False expmnt hstsbexp ccrvs
                        goalExportNamed False expmnt hstsbexp ctcrvs

                        let tcgpi = "../population-parameters/tuning-curves.gpi"
                            gpgpi = "../population-parameters/gain-profiles.gpi"
                            ccgpi = "../population-parameters/conjugacy-curves.gpi"
                            ctgpi = "../population-parameters/category-dependence.gpi"

                        mapM_ (runGnuplot expmnt hstsbexp defaultGnuplotOptions) [tcgpi,gpgpi,ccgpi,ctgpi]

                        let (pfshst,pfsft,_) = preferredStimulusHistogram nbns mlkl Nothing

                        goalExportNamed False expmnt hstsbexp pfshst
                        goalExportNamed False expmnt hstsbexp pfsft

                        let (prcshst,prcsft,_) = precisionsHistogram nbns mlkl Nothing

                        goalExportNamed False expmnt hstsbexp prcshst
                        goalExportNamed False expmnt hstsbexp prcsft

                        let (gnhsts,gnfts,_) = gainsHistograms nbns mlkl Nothing

                        mapM_ (goalExportNamed False expmnt hstsbexp) gnhsts
                        mapM_ (goalExportNamed False expmnt hstsbexp) gnfts

                        let phgpi = "../population-parameters/population-histogram.gpi"
                        runGnuplot expmnt hstsbexp defaultGnuplotOptions phgpi

                        let crsnm = "noise-correlations-" ++ show (m+1) ++ "-mixers"

                        let crsbexp = Just $ SubExperiment crsnm dst

                        let (mtxln:mtxlns) = do
                                let smlkl = sortVonMisesMixturePopulationEncoder mlkl
                                mtx <- mixturePopulationNoiseCorrelations smlkl <$> xsmps
                                return $ S.toList <$> S.toList (S.toRows mtx)

                        goalExport True expmnt crsbexp mtxln
                        mapM_ (goalExport False expmnt crsbexp) mtxlns

                        let aniopts = defaultGnuplotOptions { whetherPNG = False, whetherAnimate = True }
                            crgpi = "../population-parameters/noise-correlations.gpi"

                        runGnuplot expmnt crsbexp aniopts crgpi

                        realize $ crossValidateMixtureLikelihood kfld npop eps nbtch nepchs drch lgnrm zxs

                goalExportNamed False expmnt cesbexp cvls

        runGnuplot expmnt cesbexp defaultGnuplotOptions "cross-entropy-descent.gpi"
        runGnuplot expmnt cesbexp defaultGnuplotOptions "cross-validation.gpi"




--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

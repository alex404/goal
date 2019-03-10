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


data ValidationOpts = ValidationOpts Int Int Int Int Double Double Int Int Double Double

validationOpts :: Parser ValidationOpts
validationOpts = ValidationOpts
    <$> option auto
        ( short 'n'
        <> long "n-population"
        <> help "Number of shotgun populations to generate per training set size."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 'm'
        <> long "min-traing-set-size"
        <> help "Minimum size of the training set."
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 's'
        <> long "steps"
        <> help "Number of different sample sizes to test."
        <> showDefault
        <> value 8 )
    <*> option auto
        ( short 'r'
        <> long "repeats"
        <> help "Number of repeats per sample size."
        <> showDefault
        <> value 8 )
    <*> option auto
        ( short 'c'
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
    (ValidationOpts nshtgn mnsz nstps nsts cnc eps nbtch nepchs pmu psd) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

        lgnrm :: Natural # LogNormal
        lgnrm = toNatural . Point @ Source $ S.doubleton pmu psd

    let mxsz = nbtch * nepchs
        bas = fromIntegral (mxsz - mnsz + 1) ** (recip (fromIntegral $ nstps - 1) :: Double)
        szs = (+(mnsz-1)) . round . (bas **) . fromIntegral <$> [0..nstps-1]

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let ceanl = Just $ Analysis "upper-bound-descent" dst

        (k,m,cs) <- getFittedMixtureLikelihood expnm dst

        putStrLn "\nSize of Dataset:"

        case someNatVal k of

            SomeNat (Proxy :: Proxy k) -> case someNatVal m of

                SomeNat (Proxy :: Proxy m) ->

                    forM_ szs $ \sz -> do

                        print sz

                        let plkl :: Mean #> Natural # MixtureGLM m (Neurons k) VonMises
                            plkl = strengthenMixtureLikelihood cs

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                        (ddsts,qlkl) <- realize $ mixtureLikelihoodDataDependence
                            nsts nshtgn eps nbtch nepchs drch lgnrm plkl xsmps sz

                        goalExportNamed (sz==head szs) expmnt ceanl ddsts

                        let rltv = "../population-parameters/"
                            ttl = show sz ++ "-training-data"

                        runPopulationParameterAnalyses expmnt dst xsmps nbns rltv ttl Nothing Nothing qlkl

        --runGnuplot expmnt ceanl defaultGnuplotOptions "cross-entropy-descent.gpi"
        --runGnuplot expmnt ceanl defaultGnuplotOptions "cross-validation.gpi"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

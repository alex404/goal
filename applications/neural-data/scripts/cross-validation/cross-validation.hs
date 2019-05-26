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

nstms :: Int
nstms = 8

stms :: [Double]
stms = tail $ range mnx mxx (nstms + 1)

xsmps :: [Double]
xsmps = init $ range mnx mxx 101


--- CLI ---


data ValidationOpts = ValidationOpts Int Int Int Double NatNumber Double Int Int

validationOpts :: Parser ValidationOpts
validationOpts = ValidationOpts
    <$> option auto
        ( short 'n'
        <> long "n-population"
        <> help "Number of shotgun populations to generate."
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
        <> value 2 )
    <*> option auto
        ( short 'l'
        <> long "learning-rate"
        <> help "The learning rate."
        <> showDefault
        <> value (-0.002) )
    <*> option auto
        ( short 'b'
        <> long "n-batch"
        <> help "Batch size."
        <> showDefault
        <> value 500 )
    <*> option auto
        ( short 'e'
        <> long "n-epochs"
        <> help "Number of batches to run the learning over."
        <> showDefault
        <> value 100 )

data AllOpts = AllOpts ExperimentOpts ValidationOpts

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> validationOpts

runOpts :: AllOpts -> IO ()
runOpts ( AllOpts expopts@(ExperimentOpts expnm _)
    (ValidationOpts nshtgn kfld nmx cnc nstp _ _ nepchs) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let ceanl = Just $ Analysis "cross-validation" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do

                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs <- realize $ shuffleList zxs1

                let vtzxss = kFoldConditionalDataset kfld zxs

                let idxs = take nmx [0,nstp..]

                mcvs <- forM idxs $ \m -> case someNatVal m

                    of SomeNat (Proxy :: Proxy m) -> do

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                        let inteps = -0.05
                            intnepchs = 500
                            intnbtch = 100

                            nbtch = 50
                            eps = 0.005
                            dcy = 0.000

                        let intl = DataInitialization drch inteps intnepchs intnbtch (-0.1,0.1) zxs
                            optm = Hybrid (-eps) dcy nbtch

                        realize $ crossValidateMixtureLikelihood intl optm nshtgn nepchs vtzxss

                let cvs = fromJust <$> mcvs
                    mxcv = maximum cvs
                    igns = negate . subtract mxcv <$> cvs
                    cvstts = zipWith3 CrossValidationStats ((+1) . fromIntegral <$> idxs) cvs igns

                goalExportNamed True expmnt ceanl cvstts

        runGnuplot expmnt ceanl defaultGnuplotOptions "cross-validation.gpi"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

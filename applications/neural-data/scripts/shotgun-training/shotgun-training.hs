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
import NeuralData.Conditional.VonMises.Training

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
        <> value 0.002 )
    <*> option auto
        ( short 'b'
        <> long "n-batch"
        <> help "Batch size (for sgd)."
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
    (ValidationOpts nsht nmx cnc _ _ _ nepchs) ) = do

    dsts <- readDatasets expopts

    let expmnt = Experiment prjnm expnm

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let shtanl = Just $ Analysis "shotgun-training" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        putStrLn "\nNumber of Samples:"
        print $ length zxs0

        putStrLn "\nNumber of Neurons:"
        print . length . fst $ head zxs0

        case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do

                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs <- realize $ shuffleList zxs1

                let inteps = -0.05
                    intnepchs = 500
                    intnbtch = 100

                let lkl = last $ fitIPLikelihood inteps intnbtch intnepchs zxs

                let (zs,xs) = unzip zxs

                let ubnd = stochasticConditionalCrossEntropy xs zs lkl

                mtrulkl0 <- getMixtureLikelihood expnm dst

                if null mtrulkl0

                   then do

                       let lklbnd = LikelihoodBounds ubnd Nothing

                       goalExportNamed True expmnt shtanl $ replicate nepchs lklbnd

                   else do

                       let (_,m',cs) = fromJust mtrulkl0

                       case someNatVal m'

                           of SomeNat (Proxy :: Proxy m') -> do

                               let trulkl :: Natural #> ConditionalMixture (Neurons k) m' VonMises
                                   trulkl = strengthenMixtureLikelihood cs

                               let lklbnd = LikelihoodBounds ubnd . Just
                                       $ mixtureStochasticConditionalCrossEntropy xs zs trulkl

                               goalExportNamed True expmnt shtanl $ replicate nepchs lklbnd

                case someNatVal (nmx-1)
                    of SomeNat (Proxy :: Proxy m) -> do

                        let drch :: Natural # Dirichlet (m+1)
                            drch = Point $ S.replicate cnc

                            intts = zip
                                [ DataInitialization drch inteps intnepchs intnbtch (-0.1,0.1) zxs ]
                                [ "standard" ]

                            nbtch = 50
                            nstps = div (length zxs0) nbtch
                            eps = 0.005
                            dcy = 0.000

                            optms = zip
                                [ StochasticGradientDescent (-eps) dcy nbtch
                                , ExpectationMaximization eps nbtch nstps
                                , Hybrid (-eps) dcy nbtch
                                , Hybrid2 eps nbtch nstps ]
                                [ "sgd-e0.005-b50", "gd-em", "hybrid-em", "hybrid2-em" ]

                        sequence_ $ do

                            (intt,inttnm) <- intts
                            (optm,optmnm) <- optms

                            let algnm = inttnm ++ '-':optmnm

                            return $ do

                                putStrLn "\nAlgorithm: "
                                print algnm

                                (dscnts, nnans, mxmlkl) <- realize
                                    $ shotgunFitMixtureLikelihood intt optm nsht nepchs zxs zxs

                                putStrLn $ concat ["\nNumber of NaNs: ", show nnans , " / ", show nsht]

                                goalExportNamed False expmnt shtanl dscnts

                                let rltv = "../population-parameters/"

                                runPopulationParameterAnalyses
                                    expmnt dst xsmps nbns rltv algnm Nothing Nothing zxs mxmlkl

                runGnuplot expmnt shtanl defaultGnuplotOptions "cross-entropy-descent.gpi"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

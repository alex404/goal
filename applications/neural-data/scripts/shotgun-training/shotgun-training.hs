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
import NeuralData.Conditional.VonMises
import NeuralData.Conditional.VonMises.Training
import NeuralData.Conditional.VonMises.Analysis

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
runOpts ( AllOpts (ExperimentOpts expmnt dst) (ValidationOpts nsht nmx cnc _ _ _ nepchs) ) = do

    let ldpth = loadPath expmnt dst
        tdst = dst ++ "/true"

    (k,zxs0) <- readDataset expmnt dst

    putStrLn "\nNumber of Samples:"
    print $ length zxs0

    putStrLn "\nNumber of Neurons:"
    print k

    case someNatVal k of

        SomeNat (Proxy :: Proxy k) -> do

            let zxs1 :: [(Response k, Double)]
                zxs1 = strengthenNeuralData zxs0

            zxs <- realize $ shuffleList zxs1

            let lkl = fitIPLikelihood 0.05 500 zxs

            let ubnd = conditionalLogLikelihood zxs lkl

            flbl <- doesFileExist $ parametersPath expmnt tdst

            if flbl

               then do

                   (_,m',cs) <- readMixtureLikelihood expmnt tdst

                   case someNatVal m'

                       of SomeNat (Proxy :: Proxy m') -> do

                           let trulkl :: Natural #> ConditionalMixture (Neurons k) m' Tensor VonMises
                               trulkl = strengthenMixtureLikelihood cs

                           let lklbnd = LikelihoodBounds ubnd . Just
                                   $ conditionalLogLikelihood zxs trulkl

                           goalExportNamed ldpth "bounds" $ replicate nepchs lklbnd

               else do

                   let lklbnd = LikelihoodBounds ubnd Nothing
                   goalExportNamed ldpth "bounds" $ replicate nepchs lklbnd

            case someNatVal (nmx-1)
                of SomeNat (Proxy :: Proxy m) -> do

                    let drch :: Natural # Dirichlet (m+1)
                        drch = Point $ S.replicate cnc

                        intts = zip
                            [ DataInitialization drch 0.05 500 100 (-0.1,0.1) zxs ]
                            [ "standard" ]

                        nbtch = 50
                        eps = 0.005
                        nstps = 10

                        optms = zip
                            [ StochasticGradientDescent eps nbtch
                            , ExpectationMaximization eps nbtch nstps
                            , Hybrid eps nbtch ]
                            [ "gd-ml", "gd-em", "hybrid-em" ]

                    sequence_ $ do

                        (intt,inttnm) <- intts
                        (optm,optmnm) <- optms

                        let algnm = inttnm ++ '-':optmnm

                        return $ do

                            putStrLn "\nAlgorithm: "
                            print algnm

                            let dst' = concat [dst,"/",algnm]
                                ldpth' = loadPath expmnt dst'

                            (ascnts, nnans, mxmlkl) <- realize
                                $ shotgunFitMixtureLikelihood intt optm nsht nepchs zxs zxs

                            putStrLn $ concat ["\nNumber of NaNs: ", show nnans , " / ", show nsht]

                            goalExportNamed ldpth' "log-likelihood-ascent" ascnts
                            writeMixtureLikelihood expmnt dst' mxmlkl

                            let rltv = "../population-parameters/"

                            runPopulationParameterAnalyses
                                expmnt dst' xsmps nbns rltv Nothing Nothing zxs mxmlkl

            runGnuplot ldpth "log-likelihood-ascent"


--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

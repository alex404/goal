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


ananm :: String
ananm = "fitting-validation"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000


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
    <*> option auto
        ( short 'm'
        <> long "dirichlet"
        <> help "Dirichlet parameters (and consequently number of mixers)"
        <> showDefault
        <> value [] )
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

    let cegpi = "cross-entropy-descent.gpi"

        expmnt = Experiment prjnm expnm

        lgnrm :: Natural # LogNormal
        lgnrm = toNatural . Point @ Source $ S.doubleton pmu psd

    forM_ dsts $ \dst -> do

        putStrLn "\nDataset:"
        putStrLn dst

        let ceexp = Just $ SubExperiment ananm dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        cstss <- realize $ case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do
                let zxs1 :: [(Response k, Double)]
                    zxs1 = strengthenNeuralData zxs0
                zxs <- shuffleList zxs1
                let (tzxs,vzxs) = splitAt (round . (*vprcnt) . fromIntegral $ length zxs) zxs
                    (vzs,vxs) = unzip vzxs
                if null mxs
                   then do
                       lklss <- replicateM npop $ fitIPLikelihood eps nbtch nepchs lgnrm tzxs
                       let cost = stochasticConditionalCrossEntropy vxs vzs
                       return . L.transpose $ map cost <$> lklss
                   else case someNatVal (L.genericLength mxs - 1) of
                          SomeNat (Proxy :: Proxy m) -> do
                              let drch :: Natural # Dirichlet (m+1)
                                  drch = Point . fromJust $ S.fromList mxs
                              mlklss <- replicateM npop
                                  $ fitMixtureLikelihood eps nbtch nepchs drch lgnrm tzxs
                              let cost = mixtureStochasticConditionalCrossEntropy vxs vzs
                                  tracer mlkls =
                                      let nanbl = any isNaN . listCoordinates $ last mlkls
                                          weighter = listCoordinates . toSource . snd
                                              . splitMixtureModel . fst . splitBottomSubLinear
                                          iwghts = weighter $ head mlkls
                                          lwghts = weighter $ last mlkls
                                          trcstr = concat [ "\nAny NaNs?\n"
                                                          , show nanbl
                                                          , "\nInitial Mixture Weights:\n"
                                                          , show iwghts
                                                          , "\nFinal Mixture Weights:\n"
                                                          , show lwghts ]
                                       in trace trcstr mlkls
                              return . L.transpose $ map cost . tracer <$> mlklss

        goalExport True expmnt ceexp cstss

        runGnuplot expmnt ceexp defaultGnuplotOptions cegpi



--- Main ---


main :: IO ()
main = do

    let prgstr = "Stress test the fitting of likelihoods."
        hdrstr = "Stress test the fitting of likelihoods."
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

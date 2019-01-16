{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    BangPatterns,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData
import NeuralData.VonMises

import Goal.Core
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B


--- Analysis ---


ananm :: String
ananm = "informations"

printInformations :: forall k . KnownNat k => Proxy k -> PopulationCodeInformations -> IO ()
printInformations prxk
    (PopulationCodeInformations mimu misd lnmu lnsd lnrtomu lnrtosd
        affmu affsd affrtomu affrtosd mdcdmu mdcdsd mdcdrtomu mdcdrtosd) = do
    let k = 1 + natVal prxk
        sg = 3
    putStrLn . concat $ ["\nSubpopulation size: ", show k]
    putStrLn . concat $
        [ "Mutual Information Mean/SD: ", show $ roundSD sg mimu, ", ", show $ roundSD sg misd ]
    let printFun ttl mu sd rtomu rtosd =
            putStrLn . concat $
                [ ttl, " Mean: ", show $ roundSD sg mu
                , "; SD: ", show $ roundSD sg sd
                , "; Ratio Mean: ", show $ roundSD sg rtomu
                , "; Ratio SD: ", show $ roundSD sg rtosd ]
    printFun "Linear Posterior Divergence" lnmu lnsd lnrtomu lnrtosd
    printFun "Affine Posterior Divergence" affmu affsd affrtomu affrtosd
    when (not $ null mdcdmu) $ printFun "Affine Posterior Divergence"
        (fromJust mdcdmu) (fromJust mdcdsd) (fromJust mdcdrtomu) (fromJust mdcdrtosd)

runInformationAnalysis
    :: forall k . KnownNat k
    => Int
    -> Int
    -> Int
    -> Maybe Int
    -> Int
    -> [([Int], Double)]
    -> Proxy k
    -> IO [PopulationCodeInformations]
runInformationAnalysis nrct ncntr nmcmc mndcd nsub zxs0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0

    lkl <- realize $ fitIPLikelihood zxs

    (alldvgs0 :: B.Vector k PopulationCodeInformations) <- B.generatePM $ \prxk -> do
        pci <- realize (analyzeInformations nrct ncntr nmcmc mndcd nsub lkl prxk)
        printInformations prxk pci
        return pci

    return $ B.toList alldvgs0


--- CLI ---


data InformationOpts = InformationOpts Int Int Int Int Int

data AllOpts = AllOpts ExperimentOpts InformationOpts

informationOpts :: Parser InformationOpts
informationOpts = InformationOpts
    <$> option auto
        ( short 'r'
        <> long "n-rectification"
        <> help "Number of samples for estimating rectification parameters"
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 'c'
        <> long "n-centering"
        <> help "Number of samples for centering partition function integrals."
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 'm'
        <> long "n-monte-carlo"
        <> help "Number of samples to use to estimate general model expectations."
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 'l'
        <> long "n-decoder"
        <> help "Number of samples used to fit the linear decoder. Decoder \
            \divergences will not be computed if this number is non-positive."
        <> showDefault
        <> value 1000 )
    <*> option auto
        ( short 's'
        <> long "n-subpopulation"
        <> help "Number of subpopulations to average for each population size."
        <> showDefault
        <> value 100 )

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> informationOpts

runOpts :: AllOpts -> IO ()
runOpts (AllOpts (expopts@(ExperimentOpts expnm _)) (InformationOpts nrct ncntr nmcmc ndcd nsub)) = do

    let expmnt = Experiment prjnm expnm

    dsts <- readDatasets expopts

    let infgpi = "informations.gpi"

        mndcd = if ndcd < 1 then Nothing else Just ndcd

    forM_ dsts $ \dst -> do

        putStrLn $ "\nDataset: " ++ dst

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs = case someNatVal k of
                    SomeNat prxk -> runInformationAnalysis nrct ncntr nmcmc mndcd nsub zxs prxk

        infs <- rinfs

        let msbexp = (Just $ SubExperiment ananm dst)

        goalWriteNamedAnalysis True expmnt msbexp infs

        runGnuplot expmnt msbexp defaultGnuplotOptions infgpi


--- Main ---


main :: IO ()
main = do

    let opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr =
            "Estimate the average divergence of various approximate posteriors \
            \from the true posteriors of a model neural population.  Analyses are \
            \performed on sub populations of all sizes, to elucidate the effect \
            \of population size on divergence measures."
    runOpts =<< execParser opts

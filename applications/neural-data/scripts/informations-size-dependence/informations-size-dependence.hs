{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData
import NeuralData.VonMises

import Goal.Core
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B


--- Analysis ---


ananm :: String
ananm = "informations-size-dependence"

printInformations :: forall k . KnownNat k => Proxy k -> LogNormalInformations -> IO ()
printInformations prxk (LogNormalInformations ln _ _ aff _ _ mdcd _ _) = do
    let k = 1 + natVal prxk
    putStrLn . concat $ ["\nSubpopulation size: ", show k]
    putStrLn $ "Homogeneous Divergence Ratio: " ++ show ln
    putStrLn $ "Conjugate Divergence Ratio: " ++ show aff
    unless (null mdcd) . putStrLn $ "Decoder Posterior Divergence: " ++ show (fromJust mdcd)

runSubsamplingAnalysis
    :: forall k . KnownNat k
    => Int
    -> Int
    -> Int
    -> Maybe Int
    -> Int
    -> [([Int], Double)]
    -> Proxy k
    -> IO [LogNormalInformations]
runSubsamplingAnalysis nrct ncntr nmcmc mndcd nsub zxs0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0

    lkl <- realize $ fitIPLikelihood zxs

    (alldvgs0 :: B.Vector k LogNormalInformations) <- B.generatePM $ \prxk -> do
        pci <- realize ( logNormalInformationStatistics
                       <$> informationSubsamplingAnalysis nrct ncntr nmcmc mndcd nsub lkl prxk )
        printInformations prxk pci
        return pci

    return $ B.toList alldvgs0

runResamplingAnalysis
    :: forall k n . (KnownNat k, KnownNat n)
    => Int
    -> Int
    -> Int
    -> Maybe Int
    -> Int
    -> [([Int], Double)]
    -> Proxy k
    -> Proxy n
    -> IO [LogNormalInformations]
runResamplingAnalysis nrct ncntr nmcmc mndcd nsub zxs0 _ prxn = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0

    lkl <- realize $ fitIPLikelihood zxs

    (alldvgs0 :: B.Vector n LogNormalInformations) <- B.generatePM $ \prxn -> do
        pci <- realize ( logNormalInformationStatistics
                       <$> informationResamplingAnalysis nrct ncntr nmcmc mndcd nsub lkl prxn)
        printInformations prxn pci
        return pci

    return $ B.toList alldvgs0


--- CLI ---


data InformationOpts = InformationOpts Int Int Int Int Int NatNumber Bool

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
        <> help "Number of populations to average for each population size."
        <> showDefault
        <> value 100 )
    <*> option auto
        ( short 'R'
        <> long "n-resample"
        <> help "Max population size for resample analysis."
        <> showDefault
        <> value 200 )
    <*> switch
        ( short 'S'
        <> long "subsample"
        <> help "Run subsample analysis (default is resample analysis)." )


allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> informationOpts

runOpts :: AllOpts -> IO ()
runOpts (AllOpts expopts@(ExperimentOpts expnm _) (InformationOpts nrct ncntr nmcmc ndcd nsub npop sbl)) = do

    let expmnt = Experiment prjnm expnm

    dsts <- readDatasets expopts

    let infgpi = "informations-size-dependence.gpi"

        mndcd = if ndcd < 1 then Nothing else Just ndcd

    forM_ dsts $ \dst -> do

        putStrLn $ "\nDataset: " ++ dst

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs =
                case someNatVal k of
                       SomeNat prxk ->
                           if sbl
                              then runSubsamplingAnalysis nrct ncntr nmcmc mndcd nsub zxs prxk
                              else case someNatVal npop of
                                     SomeNat prxn ->
                                         runResamplingAnalysis nrct ncntr nmcmc mndcd nsub zxs prxk prxn

        infs <- rinfs

        let msbexp = Just $ SubExperiment ananm dst

        goalExportNamed True expmnt msbexp infs

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

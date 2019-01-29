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
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Analysis ---


ananm :: String
ananm = "resample-informations"


--- CLI ---


data InformationOpts = InformationOpts Int Int Int Int Int Int

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
        <> value 200 )
    <*> option auto
        ( short 'l'
        <> long "n-decoder"
        <> help "Number of samples used to fit the linear decoder. Decoder \
            \divergences will not be computed if this number is non-positive."
        <> showDefault
        <> value 0 )
    <*> option auto
        ( short 's'
        <> long "n-resamples"
        <> help "Number of populations to sample."
        <> showDefault
        <> value 2000 )
    <*> option auto
        ( short 'b'
        <> long "n-bins"
        <> help "Number of bins in the histogram."
        <> showDefault
        <> value 20 )

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> informationOpts

runInformationAnalysis
    :: forall k r . KnownNat k
    => Int
    -> Int
    -> Int
    -> Maybe Int
    -> Int
    -> Int
    -> [([Int], Double)]
    -> Proxy k
    -> Random r (Informations,[InformationCounts])
runInformationAnalysis nrct ncntr nmcmc mndcd npop nbns zxs0 prxk = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0

    lkl <- fitIPLikelihood zxs

    inf <- estimateInformations nrct ncntr nmcmc mndcd lkl

    let (_,_,pdnsprms) = unzip3 $ populationParameters nbns lkl
    let (ParameterDensityParameters gmu gsd) = head pdnsprms
        (ParameterDensityParameters pmu psd) = pdnsprms !! 2
    let rgns = Point $ S.fromTuple (gmu, gsd)
        rprcs = Point $ S.fromTuple (pmu, psd)

    infs <- informationResamplingAnalysis nrct ncntr nmcmc mndcd npop rgns rprcs prxk
    return (inf,histogramInformationStatistics nbns infs)

runOpts :: AllOpts -> IO ()
runOpts (AllOpts expopts@(ExperimentOpts expnm _) (InformationOpts nrct ncntr nmcmc ndcd nsmp nbns)) = do

    let expmnt = Experiment prjnm expnm

    dsts <- readDatasets expopts

    let infgpi = "resample-informations.gpi"

        mndcd = if ndcd < 1 then Nothing else Just ndcd

    forM_ dsts $ \dst -> do

        putStrLn $ "\nDataset: " ++ dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs = case someNatVal k of
                    SomeNat prxk -> realize
                        $ runInformationAnalysis nrct ncntr nmcmc mndcd nsmp nbns zxs0 prxk

        (inf,infs) <- rinfs

        let msbexp = Just $ SubExperiment ananm dst

        goalExportNamed True expmnt msbexp infs
        goalExportNamed False expmnt msbexp [informationsToRatios inf]


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

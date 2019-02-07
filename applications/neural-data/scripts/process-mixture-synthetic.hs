{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeApplications,
    TypeOperators
    #-}


--- Imports ---


-- Goal --

import NeuralData
import NeuralData.VonMises

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
xsmps = tail $ range mnx mxx 100

ananm :: String
ananm = "true-tuning-curves"


-- Convolutional --


cnvmus :: KnownNat k => S.Vector k Double
cnvmus = S.init $ S.range mnx mxx

cnvps :: KnownNat k => Double -> S.Vector k (Natural # VonMises)
cnvps pmu = S.map (toNatural . Point @ Source . flip S.doubleton pmu) cnvmus

cnvlkl :: KnownNat k => Double -> Double -> Mean #> Natural # Neurons k <* VonMises
cnvlkl gmu pmu = vonMisesPopulationEncoder True (Left (log gmu)) (cnvps pmu)

-- Random --

rlklr
    :: KnownNat k
    => Double
    -> Double
    -> Double
    -> Double
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
rlklr lgmu lgvr lpmu lpvr = do
    let rgns = Point $ S.doubleton lgmu lgvr
        rprcs = Point $ S.doubleton lpmu lpvr
    randomLikelihood rgns zero rprcs

mcnvlklr
    :: KnownNat k
    => Double
    -> Double
    -> Double
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
mcnvlklr lgmu lgvr pmu = do
    gns <- randomGains . Point $ S.doubleton lgmu lgvr
    return . vonMisesPopulationEncoder True (Right $ transition gns) $ cnvps pmu

normalizeLikelihood
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Mean #> Natural # Neurons k <* VonMises
normalizeLikelihood lkl0 =
    let (nz,nzx) = splitAffine lkl0
        bnd = 0.0001
        eps = -0.005
        cauchify = last . take 10000 . cauchySequence euclideanDistance bnd
        rho0 = average $ potential <$> lkl0 >$>* xsmps
        diff = populationCodeRectificationDifferential rho0 zero xsmps nzx
        nz' = cauchify $ vanillaGradientSequence diff eps defaultAdamPursuit nz
     in joinAffine nz' nzx

combineStimuli :: [[Response k]] -> [([Int],Double)]
combineStimuli zss =
    concat $ zipWith (\zs x -> zip (toList <$> zs) $ repeat x) zss stms

-- IO --

synthesizeData :: forall k . KnownNat k
               => String -> Proxy k -> NatNumber -> Double -> Double -> Double -> Double -> IO ()
synthesizeData expnm prxk nsmps0 lgmu lgsd lpmu lpsd = do

    let lgvr = square lgsd
        gmu = exp $ lgmu + lgvr/2
        lpvr = square lpsd
        pmu = exp $ lpmu + lgvr/2

    putStrLn "\nMean Gain:"
    print gmu
    putStrLn "\nMean Precision:"
    print pmu

    let nsmps = round (fromIntegral nsmps0 / fromIntegral nstms :: Double)
        k = natVal prxk
        expmnt = Experiment prjnm expnm

    rlkl <- realize $ rlklr lgmu lgvr lpmu lpvr

    mcnvlkl <- realize $ mcnvlklr lgmu lgvr pmu

    let nrlkl :: Mean #> Natural # Neurons k <* VonMises
        nrlkl = normalizeLikelihood rlkl

    let cnvlkln = cnvlkl gmu pmu

    let dsts = ["convolutional","modulated-convolutional","random","random-normalized"]

    goalWriteDatasetsCSV expmnt dsts

    let tcgpi = "population-parameters/tuning-curves.gpi"
        ppgpi = "population-parameters/population-parameter-histogram.gpi"

    sequence_ $ do

        (lkl,dst) <- zip [cnvlkln,mcnvlkl,rlkl,nrlkl] dsts

        return $ do

            putStrLn $ concat ["\nDataset: ", dst, "\n"]

            infs <- realize (estimateInformations 1000 1000 1000 Nothing lkl)
            print infs

            (zss :: [[Response k]]) <- realize (mapM (sample nsmps) $ lkl >$>* stms)

            let zxs :: [([Int], Double)]
                zxs = combineStimuli zss

            goalWriteDataset expmnt dst $ show (k,zxs)

            let msbexpt = Just $ SubExperiment "true-tuning-curves" dst

            goalExport True expmnt msbexpt $ analyzeTuningCurves xsmps lkl

            runGnuplot expmnt msbexpt defaultGnuplotOptions tcgpi

            let msbexph = Just $ SubExperiment "true-histograms" dst

            let (hstcsv:hstcsvs,ppds) = unzip . fst $ populationParameters 20 lkl

            goalExportNamed True expmnt msbexph hstcsv
            mapM_ (goalExportNamed False expmnt msbexph) hstcsvs
            mapM_ (goalExportNamed False expmnt msbexph) ppds

            runGnuplot expmnt msbexph defaultGnuplotOptions ppgpi


--- CLI ---


data SyntheticOpts = SyntheticOpts String NatNumber NatNumber Double Double Double Double

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> strArgument
        ( metavar "PROJ"
        <> help "Name of this synthetic experiment."
        <> showDefault
        <> value "synthetic" )
    <*> option auto
        ( short 'k'
        <> long "k-neurons"
        <> help "Number of neurons in the model population."
        <> showDefault
        <> value 20 )
    <*> option auto
        ( short 'n'
        <> long "n-samples"
        <> help "Number of samples to generate from the model."
        <> showDefault
        <> value 400 )
    <*> option auto
        ( short 'g'
        <> long "log-mu-gain"
        <> help "The mu parameter of the gain log-normal."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'G'
        <> long "log-sd-gain"
        <> help "The sd parameter the gain log-normal."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'p'
        <> long "log-mu-precision"
        <> help "The mu parameter of the precision log-normal."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'P'
        <> long "log-sd-precision"
        <> help "The sd parameter of the precision log-normal."
        <> showDefault
        <> value 1 )

runOpts :: SyntheticOpts -> IO ()
runOpts (SyntheticOpts expnm k nsmps gmu lgsd pmu lpsd) =
    case someNatVal k of
      SomeNat prxk -> synthesizeData expnm prxk nsmps gmu lgsd pmu lpsd


--- Main ---


main :: IO ()
main = do

    let hdrstr = "Generate synthetic datasets from populations of Poisson neurons."
        prgstr =
            "Generate synthetic data from a model population of \
            \stimulus-dependent, Poisson neurons. Model generation parameters \
            \can be specified with command line arguments."
        opts = info (syntheticOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)

    runOpts =<< execParser opts

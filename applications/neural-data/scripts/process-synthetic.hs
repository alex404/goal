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


cmus :: KnownNat k => S.Vector k Double
cmus = S.init $ S.range mnx mxx

cnps :: KnownNat k => Double -> S.Vector k (Natural # VonMises)
cnps pmu = S.map (toNatural . Point @ Source . flip S.doubleton pmu) cmus

clkl :: KnownNat k => Double -> Double -> Mean #> Natural # Neurons k <* VonMises
clkl gmu pmu = vonMisesPopulationEncoder True (Left (log gmu)) (cnps pmu)

-- Random --

rmus :: KnownNat k => Random r (S.Vector k Double)
rmus = S.replicateM $ uniformR (mnx,mxx)

rkps :: KnownNat k => Double -> Double -> Random r (S.Vector k Double)
rkps lpmu lpvr = S.replicateM $ exp <$> samplePoint generator
    where generator :: Source # Normal
          generator = Point $ S.fromTuple (lpmu,lpvr)

rnps :: KnownNat k => Double -> Double -> Random r (S.Vector k (Natural # VonMises))
rnps lpmu lpvr = do
    mus <- rmus
    kps <- rkps lpmu lpvr
    let mukps = S.zipWith S.doubleton mus kps
    return $ S.map (toNatural . Point @ Source) mukps

rgns :: KnownNat k => Double -> Double -> Random r (Natural # Neurons k)
rgns lgmu lgvr = initialize generator
    where generator :: Source # Normal
          generator = Point $ S.fromTuple (lgmu,lgvr)

rlklr
    :: KnownNat k
    => Double
    -> Double
    -> Double
    -> Double
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
rlklr lgmu lgvr lpmu lpvr = do
    gns <- rgns lgmu lgvr
    nps <- rnps lpmu lpvr
    return $ vonMisesPopulationEncoder True (Right gns) nps

mclklr
    :: KnownNat k
    => Double
    -> Double
    -> Double
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
mclklr lgmu lgvr pmu = do
    gns <- rgns lgmu lgvr
    return . vonMisesPopulationEncoder True (Right gns) $ cnps pmu

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
synthesizeData expnm prxk nsmps0 gmu lgsd pmu lpsd = do

    let lgvr = square lgsd
        lgmu = log gmu - lgvr/2
        lpvr = square lpsd
        lpmu = log pmu - lgvr/2

    let nsmps = round $ (fromIntegral nsmps0 / fromIntegral nstms :: Double)
        k = natVal prxk
        expmnt = Experiment prjnm expnm

    rlkl <- realize $ rlklr lgmu lgvr lpmu lpvr

    mclkl <- realize $ mclklr lgmu lgvr pmu

    let nrlkl :: Mean #> Natural # Neurons k <* VonMises
        nrlkl = normalizeLikelihood rlkl

    let clkln = clkl gmu pmu

    let dsts = ["convolutional","modulated-convolutional","random","random-normalized"]

    goalWriteDatasetsCSV expmnt dsts

    let tcgpi = "population-parameters/tuning-curves.gpi"
        ppgpi = "population-parameters/population-parameter-histogram.gpi"

    sequence_ $ do

        (lkl,dst) <- zip [clkln,mclkl,rlkl,nrlkl] dsts

        return $ do

            (zss :: [[Response k]]) <- realize (mapM (sample nsmps) $ lkl >$>* stms)

            let zxs :: [([Int], Double)]
                zxs = combineStimuli zss

            goalWriteDataset expmnt dst $ show (k,zxs)

            let msbexpt = Just $ SubExperiment "true-tuning-curves" dst

            goalWriteAnalysis True expmnt msbexpt $ analyzeTuningCurves xsmps lkl

            runGnuplot expmnt msbexpt defaultGnuplotOptions tcgpi

            let msbexph = Just $ SubExperiment "true-histograms" dst

            let (prmcsv:prmcsvs) = populationParameterHistogram nbns lkl

            goalWriteNamedAnalysis True expmnt msbexph prmcsv
            mapM_ (goalWriteNamedAnalysis False expmnt msbexph) prmcsvs

            runGnuplot expmnt msbexph defaultGnuplotOptions ppgpi

--- CLI ---


data SyntheticOpts = SyntheticOpts String NatNumber NatNumber Double Double Double Double

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> option auto
        ( short 'e'
        <> long "experiment-name"
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
        <> long "mean-gain"
        <> help "The average gain of the neurons."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 'G'
        <> long "sd-log-gain"
        <> help "The standard deviation of the log-gains."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'p'
        <> long "mean-precision"
        <> help "The average precision of the tuning curves."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'P'
        <> long "sd-log-precision"
        <> help "The standard deviation of the log tuning-precision."
        <> showDefault
        <> value 0.5 )

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

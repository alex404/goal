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


-- Population Generation --


-- Convolutional --

convolutionalMus :: KnownNat k => S.Vector k Double
convolutionalMus = S.init $ S.range mnx mxx

convolutionalTuningCurves
    :: KnownNat k
    => Double -- ^ Tuning Curve Precision
    -> S.Vector k (Natural # VonMises)
convolutionalTuningCurves prcs = S.map (toNatural . Point @ Source . flip S.doubleton prcs) convolutionalMus

convolutionalLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural # Dirichlet (m+1) -- ^ Mixture parameters
    -> S.Vector (m+1) Double -- ^ Global gains
    -> Double -- ^ Global precision
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) Int m VonMises)
convolutionalLikelihood drch gns prcs = do
    cts <- samplePoint drch
    let nctgl = toNatural $ Point @ Source cts
        tcs = convolutionalTuningCurves prcs
        lgns = S.map (Point . S.replicate . log) gns
    return $ vonMisesMixturePopulationEncoder True nctgl lgns tcs

-- Modulated Convolutional --

modulatedConvolutionalLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural # Categorical Int m
    -> S.Vector (m+1) (Natural # Normal) -- ^ Log-Gain Distributions
    -> Double -- ^ Global precision
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) Int m VonMises)
modulatedConvolutionalLikelihood nctgl rgns prcs = do
    lgns <- S.mapM (fmap Point . S.replicateM . samplePoint) rgns
    return . vonMisesMixturePopulationEncoder True nctgl lgns $ convolutionalTuningCurves prcs

-- Random --

randomLikelihood
    :: KnownNat k
    => Double
    -> Double
    -> Double
    -> Double
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
randomLikelihood lgmu lgvr lpmu lpvr = do
    let rgns = Point $ S.doubleton lgmu lgvr
        rprcs = Point $ S.doubleton lpmu lpvr
    randomIPLikelihood rgns zero rprcs

-- Normalized --

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


--- Functions ---


-- Utility --

combineStimuli :: [[Response k]] -> [([Int],Double)]
combineStimuli zss =
    concat $ zipWith (\zs x -> zip (toList <$> zs) $ repeat x) zss stms

-- IO --

synthesizeData :: forall k m . (KnownNat k , KnownNat m)
               => Proxy m -- ^ Number of mixers
               -> Proxy k -- ^ Population size
               -> [Double] -- ^ Gains
               -> [Double] -- ^ Concetrantions
               -> String -- ^ Experiment name
               -> NatNumber -- ^ Number of population samples
               -> Double -- ^ log-gain sd
               -> Double -- ^ log-precision mean
               -> Double -- ^ log-precision sd
               -> IO ()
synthesizeData prxm prxk gs alphs expnm nsmps0 lgsd lpmu lpsd = do

    let lgvr = square lgsd
        gmus = [ exp $ lgmu + lgvr/2 | lgmu <- gs ]
        lpvr = square lpsd
        pmu = exp $ lpmu + lgvr/2

    putStrLn "\nMean Gains:"
    print gmus
    putStrLn "\nMean Precision:"
    print pmu

    let nsmps = round (fromIntegral nsmps0 / fromIntegral nstms :: Double)
        k = natVal prxk
        expmnt = Experiment prjnm expnm

    rlkl <- realize $ randomLikelihood lgmu lgvr lpmu lpvr

    mcnvlkl <- realize $ modulatedConvolutionalLikelihood lgmu lgvr pmu

    let nrlkl :: Mean #> Natural # Neurons k <* VonMises
        nrlkl = normalizeLikelihood rlkl

    let cnvlkln = convolutionalLikelihood gmu pmu

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


data SyntheticOpts = SyntheticOpts [String] String NatNumber NatNumber Double Double Double

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> many (strArgument
        ( metavar "GAIN,CONCENTRATION"
        <> showDefault
        <> value "10,1" ) )
    <*> option auto
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
        ( short 'G'
        <> long "gain-sd"
        <> help "The sd of the log-gains."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'p'
        <> long "precision-mu"
        <> help "The mean precision."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 'P'
        <> long "log-precision-sd"
        <> help "The standard deviation of the log-precisions."
        <> showDefault
        <> value 1 )

runOpts :: SyntheticOpts -> IO ()
runOpts (SyntheticOpts galphstrs expnm k nsmps lgsd pmu lpsd) =
    let (gs,alphs) = unzip [ read $ '(' ++ galphstr ++ ')' | galphstr <- galphstrs ]
        m = length gs
     in case someNatVal m of
          SomeNat prxm ->
              case someNatVal k of
                  SomeNat prxk -> do
                      print gmus
                      synthesizeData prxm prxk gs alphs expnm nsmps gmu lgsd pmu lpsd


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

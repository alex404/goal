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
import NeuralData.Mixture

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L


--- Globals ---


nbns :: Int
nbns = 10

xsmps :: [Double]
xsmps = init $ range mnx mxx 101

-- Population Generation --

randomCategorical :: KnownNat m => Natural # Dirichlet (m+1) -> Random r (Natural # Categorical Int m)
randomCategorical drch = do
    cts <- samplePoint drch
    return $ toNatural . Point @ Source $ S.tail cts

-- Convolutional --

convolutionalMus :: KnownNat k => S.Vector k Double
convolutionalMus = S.init $ S.range mnx mxx

convolutionalTuningCurves
    :: KnownNat k
    => Double -- ^ Tuning Curve Precision
    -> S.Vector k (Natural # VonMises)
convolutionalTuningCurves prcs =
    S.map (toNatural . Point @ Source . flip S.doubleton prcs) convolutionalMus

convolutionalLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural # Dirichlet (m+1) -- ^ Mixture parameters
    -> S.Vector (m+1) Double -- ^ Global gains
    -> Double -- ^ Global precision
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
convolutionalLikelihood drch gns prcs = do
    nctgl <- randomCategorical drch
    let tcs = convolutionalTuningCurves prcs
        lgns = S.map (Point . S.replicate . log) gns
    return $ joinVonMisesMixturePopulationEncoder nctgl lgns tcs

-- Modulated Convolutional --

modulatedConvolutionalLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural # Dirichlet (m+1) -- ^ Mixture parameters
    -> S.Vector (m+1) (Natural # Normal) -- ^ Log-Gain Distributions
    -> Double -- ^ Global precision
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
modulatedConvolutionalLikelihood drch rngnss prcs = do
    nctgl <- randomCategorical drch
    ngnss <- S.mapM initialize rngnss
    return . joinVonMisesMixturePopulationEncoder nctgl ngnss $ convolutionalTuningCurves prcs

-- Random --

somewhatRandomTuningCurves
    :: KnownNat k
    => Natural # LogNormal -- ^ Tuning Curve Precision Distribution
    -> Random r (S.Vector k (Natural # VonMises))
somewhatRandomTuningCurves rprc = do
    prcs <- S.replicateM $ samplePoint rprc
    return $ S.zipWith tcFun convolutionalMus prcs
        where tcFun mu prc = toNatural . Point @ Source $ S.doubleton mu prc

randomLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural # Dirichlet (m+1) -- ^ Mixture parameter distribution
    -> S.Vector (m+1) (Natural # Normal) -- ^ Log-Gain Distributions
    -> Natural # LogNormal -- ^ Precision Distribution
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) m VonMises)
randomLikelihood drch rngnss rprc = do
    nctgl <- randomCategorical drch
    ngnss <- S.mapM initialize rngnss
    ntcs <- somewhatRandomTuningCurves rprc
    return $ joinVonMisesMixturePopulationEncoder nctgl ngnss ntcs

-- Conjugated --

conjugateLikelihood
    :: (KnownNat k, KnownNat m)
    => Mean #> Natural # MixtureGLM (Neurons k) m VonMises
    -> Mean #> Natural # MixtureGLM (Neurons k) m VonMises
conjugateLikelihood lkl0 =
    let (nz,nzx) = splitBottomSubLinear lkl0
        bnd = 0.001
        eps = -0.05
        cauchify = last . take 10000 . cauchySequence euclideanDistance bnd
        rho0 = average $ potential <$> lkl0 >$>* xsmps
        diff = conditionalHarmoniumConjugationDifferential rho0 zero xsmps nzx
        nz' = cauchify $ vanillaGradientSequence diff eps defaultAdamPursuit nz
     in joinBottomSubLinear nz' nzx



--- Functions ---


-- Utility --


-- IO --

synthesizeData
    :: forall k m . (KnownNat k , KnownNat m)
    => String -- ^ Experiment name
    -> Proxy k -- ^ Population size
    -> Int -- ^ Number of unique stimuli
    -> S.Vector (m+1) Double -- ^ Gains
    -> S.Vector (m+1) Double -- ^ Concentrations
    -> Double -- ^ log-gain sd
    -> Double -- ^ precision mean
    -> Double -- ^ log-precision sd
    -> NatNumber -- ^ Number of population samples
    -> IO ()
synthesizeData expnm prxk nstms gnmus alphs lgnsd prcmu lprcsd nsmps0 = do

    let stms = init $ range mnx mxx (nstms + 1)

        combineStimuli :: [[Response k]] -> [([Int],Double)]
        combineStimuli zss =
            concat $ zipWith (\zs x -> zip (toList <$> zs) $ repeat x) zss stms

        combineStimuli' :: [[Response k]] -> [(Response k,Double)]
        combineStimuli' zss =
            concat $ zipWith (\zs x -> zip zs $ repeat x) zss stms


    let lgnvr = square lgnsd
        lprcvr = square lprcsd

    let lgnmus = S.map (\gnmu -> log gnmu - lgnvr/2) gnmus
        lprcmu = log prcmu - lprcvr/2

    let drch = Point @ Natural alphs
        lgnnrms :: S.Vector (m+1) (Natural # Normal)
        lgnnrms = S.map (toNatural . Point @ Source . flip S.doubleton lgnvr) lgnmus
        prclnrm :: Natural # LogNormal
        prclnrm = toNatural . Point @ Source $ S.doubleton lprcmu lprcvr

    let nsmps = round (fromIntegral nsmps0 / fromIntegral nstms :: Double)
        k = natVal prxk
        m = natVal (Proxy @ m)
        expmnt = Experiment prjnm expnm

    rlkl <- realize $ randomLikelihood drch lgnnrms prclnrm

    mcnvlkl <- realize $ modulatedConvolutionalLikelihood drch lgnnrms prcmu

    let nrlkl = conjugateLikelihood rlkl

    cnvlkln <- realize $ convolutionalLikelihood drch gnmus prcmu

    let dsts = ["random","convolutional","modulated-convolutional","conjugated-random"]

    goalWriteDatasetsCSV expmnt dsts

    let mgndstss = [Just $ S.toList lgnnrms,Nothing,Just $ S.toList lgnnrms,Nothing]
        mprcsdsts = [Just prclnrm,Nothing,Nothing,Just prclnrm]

    sequence_ $ do

        (lkl,dst,mgndsts,mprcsdst) <- L.zip4 [rlkl,cnvlkln,mcnvlkl,nrlkl] dsts mgndstss mprcsdsts

        return $ do

            putStrLn $ concat ["\nDataset: ", dst, "\n"]

            --infs <- realize (estimateInformations 1000 1000 1000 Nothing lkl)
            --print infs

            (zss :: [[Response k]]) <- realize $ do
                smps <- mapM (sample nsmps) $ lkl >$>* stms
                return $ map hHead <$> smps

            let zxs = combineStimuli zss

            let zxs' = combineStimuli' zss

            goalWriteDataset expmnt dst $ show (k,zxs)
            goalWriteDataset expmnt (dst ++ "-parameters") $ show (k,m,listCoordinates lkl)

            let mgndsts' = map breakPoint <$> mgndsts

            runPopulationParameterAnalyses expmnt dst xsmps nbns "population-parameters/" "true" mprcsdst mgndsts' zxs' lkl



--- CLI ---


data SyntheticOpts = SyntheticOpts [Double] String Double NatNumber NatNumber Int Double Double Double

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> many (argument auto
        ( metavar "GAINS..."))
    <*> option str
        ( short 'e'
        <> long "experiment-name"
        <> help "Name of this synthetic experiment."
        <> showDefault
        <> value "synthetic" )
    <*> option auto
        ( short 'c'
        <> long "concentration"
        <> help "Gain profile probability concentration (higher means less variable)."
        <> showDefault
        <> value 5 )
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
        ( short 's'
        <> long "n-stimuli"
        <> help "Number of unique stimuli to generate synthetic data from."
        <> showDefault
        <> value 8 )
    <*> option auto
        ( short 'G'
        <> long "gain-sd"
        <> help "The sd of the log-gains."
        <> showDefault
        <> value 0.5 )
    <*> option auto
        ( short 'p'
        <> long "precision-mu"
        <> help "The mean precision."
        <> showDefault
        <> value 0.5 )
    <*> option auto
        ( short 'P'
        <> long "log-precision-sd"
        <> help "The standard deviation of the log-precisions."
        <> showDefault
        <> value 1 )

runOpts :: SyntheticOpts -> IO ()
runOpts (SyntheticOpts gmus0 expnm cnc k nsmps nstms lgsd pmu lpsd) = do
    let m1 = fromIntegral $ length gmus0 - 1
    case someNatVal m1 of
      SomeNat (_ :: Proxy m) ->
          let gmus :: S.Vector (m+1) Double
              gmus = fromJust $ S.fromList gmus0
              alphs :: S.Vector (m+1) Double
              alphs = S.replicate cnc
           in case someNatVal k of
                SomeNat prxk ->
                    synthesizeData expnm prxk nstms gmus alphs lgsd pmu lpsd nsmps


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

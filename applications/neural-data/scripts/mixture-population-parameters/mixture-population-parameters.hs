{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData
import NeuralData.VonMises

import Paths_neural_data

import Goal.Core
import Goal.Probability


--- Globals ---


ananm :: String
ananm = "tuning-curves"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000


--- CLI ---


data AnalysisOpts = AnalysisOpts String String

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument ( help "Which data collection to analyze" )
    <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
--    <*> option auto
--        (long "tck" <> help "subsampled population size" <> short 'k' <> value 0)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg) = do

    dsts <- if null dstarg
               then fromJust <$> goalReadDatasetsCSV (Experiment prjnm expnm)
               else return [dstarg]

    tcgpi <- getDataFileName "tuning-curves.gpi"
    ppgpi <- getDataFileName "population-parameter-histogram.gpi"

    let expmnt = Experiment prjnm expnm


    forM_ dsts $ \dst -> do

        let msbexpt = Just $ SubExperiment "tuning-curves" dst
            msbexph = Just $ SubExperiment "histograms" dst

        (k,zxs0 :: [([Int], Double)]) <- getNeuralData expnm dst

        (tcss,hstcsv:hstcsvs) <- realize $ case someNatVal k of
            SomeNat (Proxy :: Proxy k) -> do
                let zxs :: [(Response k, Double)]
                    zxs = strengthenNeuralData zxs0
                lkl <- fitIPLikelihood zxs
                return (analyzeTuningCurves xsmps lkl,populationParameterHistogram 10 lkl)

        goalWriteAnalysis True expmnt msbexpt tcss

        runGnuplot expmnt msbexpt defaultGnuplotOptions tcgpi

        goalWriteNamedAnalysis True expmnt msbexph hstcsv
        mapM_ (goalWriteNamedAnalysis False expmnt msbexph) hstcsvs

        runGnuplot expmnt msbexph defaultGnuplotOptions ppgpi



--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts


--- Graveyard ---


{-

analyzeMixtureTuningCurves
    :: forall k . KnownNat k
    => [([Int],Double)]
    -> Proxy k
    -> [[Double]]
analyzeMixtureTuningCurves zxs0 _ =
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
        ppc = fitIPLikelihood zxs
        (rho0,rprms) = regressRectificationParameters ppc xsmps
        rcrv = rectificationCurve rho0 rprms xsmps
        tcs = listCoordinates . dualTransition <$> ppc >$>* xsmps
        mxtcs = maximum <$> tcs
--        (mupr,vrpr) = estimateMeanVariance . map (head . tail . listCoordinates . toSource)
--            . S.toList . toRows . snd $ splitAffine ppc
        stdoustr = concat ["mupr: ", show mupr, "; sdpr: ", show $ sqrt vrpr]
    in trace stdoustr $ zipWith (++) (L.transpose [xsmps,potential <$> ppc >$>* xsmps,rcrv,mxtcs]) tcs

analyzeMixtureTuningCurves0
    :: (KnownNat k, KnownNat n)
    => [Natural # Harmonium Tensor (Neurons k) (Categorical Int n)]
    -> Int
    -> ([Double],[[Double]])
analyzeMixtureTuningCurves0 nzks i =
    let nzs = (>.>* i) . fst . splitBottomHarmonium <$> nzks
        tcs = listCoordinates . dualTransition <$> nzs
    in (potential <$> nzs, zipWith (:) xsmps tcs)

analyzeMixtureTuningCurves1
    :: forall n k . (KnownNat k, KnownNat n)
    => Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
    -> Proxy k
    -> Proxy n
    -> ([[Double]],[[Double]],[[[Double]]])
analyzeMixtureTuningCurves1 mppc _ prxn =
    let nzks = mppc >$>* xsmps
        idxs = [0.. natValInt prxn - 1]
        (stcss,tcss) = unzip $ analyzeMixtureTuningCurves0 nzks <$> idxs
        (rho0,rprms) = regressRectificationParameters mppc xsmps
        rcrv = rectificationCurve rho0 rprms xsmps
        stcs' = potential <$> nzks
        wghts = [ density cat <$> idxs | cat <- snd . splitMixtureModel <$> nzks ]
     in (zipWith (:) xsmps wghts,L.transpose $ xsmps:rcrv:stcs':stcss, tcss)

analyzeMixtureTuningCurves
    :: forall n k . (KnownNat k, KnownNat n)
    => [Double]
    -> Proxy k
    -> Proxy n
    -> ([[Double]],[[Double]],[[[Double]]])
analyzeMixtureTuningCurves cs prxk prxn =
    let mppc = strengthenMixtureLikelihood cs
     in analyzeMixtureTuningCurves1 mppc prxk prxn

fitAnalyzeMixtureTuningCurves
    :: forall n k r . (KnownNat k, KnownNat n)
    => [([Int], Double)]
    -> Proxy k
    -> Proxy n
    -> Random r ([[Double]],[[Double]],[[[Double]]],[Double])
fitAnalyzeMixtureTuningCurves zxs0 prxk prxn = do
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
    mppc <- fitMixtureLikelihood zxs
    let cvrs = estimateCorrelations $ fst <$> zxs
    let (wghts,stcs,tcss) = analyzeMixtureTuningCurves1 mppc prxk prxn
    return (wghts,stcs,tcss,cvrs)
-}

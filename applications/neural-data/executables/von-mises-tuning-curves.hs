{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Data.Semigroup ((<>))
import qualified Data.List as L


--- Globals ---


ananm :: String
ananm = "tuning-curves"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000


--- Functions ---


analyzeTuningCurves0
    :: forall k . KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Proxy k
    -> [[Double]]
analyzeTuningCurves0 lkl _ =
    let nzs = lkl >$>* xsmps
        tcs = listCoordinates . dualTransition <$> nzs
        stcs = potential <$> nzs
        (rho0,rprms) = regressRectificationParameters lkl xsmps
        rcrv = rectificationCurve rho0 rprms xsmps
     in L.transpose $ xsmps:rcrv:stcs:tcs

analyzeTuningCurves
    :: forall k . KnownNat k
    => [Double]
    -> Proxy k
    -> [[Double]]
analyzeTuningCurves cs prxk =
    let lkl = strengthenIPLikelihood cs
     in analyzeTuningCurves0 lkl prxk

fitAnalyzeTuningCurves
    :: forall k r . KnownNat k
    => [([Int], Double)]
    -> Proxy k
    -> Random r [[Double]]
fitAnalyzeTuningCurves zxs0 prxk = do
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
    lkl <- fitIPLikelihood zxs
    return $ analyzeTuningCurves0 lkl prxk

analyzeMixtureTuningCurves0
    :: (KnownNat k, KnownNat n, 1 <= n)
    => [Natural # Harmonium Tensor (Neurons k) (Categorical Int n)]
    -> Int
    -> ([Double],[[Double]])
analyzeMixtureTuningCurves0 nzks i =
    let nzs = (>.>* i) . fst . splitBottomHarmonium <$> nzks
        tcs = listCoordinates . dualTransition <$> nzs
    in (potential <$> nzs, zipWith (:) xsmps tcs)

analyzeMixtureTuningCurves1
    :: forall n k . (KnownNat k, KnownNat n, 1 <= n)
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
    :: forall n k . (KnownNat k, KnownNat n, 1 <= n)
    => [Double]
    -> Proxy k
    -> Proxy n
    -> ([[Double]],[[Double]],[[[Double]]])
analyzeMixtureTuningCurves cs prxk prxn =
    let mppc = strengthenMixtureLikelihood cs
     in analyzeMixtureTuningCurves1 mppc prxk prxn

fitAnalyzeMixtureTuningCurves
    :: forall n k r . (KnownNat k, KnownNat n, 1 <= n)
    => [([Int], Double)]
    -> Proxy k
    -> Proxy n
    -> Random r ([[Double]],[[Double]],[[[Double]]])
fitAnalyzeMixtureTuningCurves zxs0 prxk prxn = do
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
    mppc <- fitMixtureLikelihood zxs
    return $ analyzeMixtureTuningCurves1 mppc prxk prxn


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument ( help "Which data collection to plot" )
    <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nmixers" <> short 'm' <> help "Number of mixers to use. Only valid for response data, not true tuning curves" <> value 4)
--    <*> option auto
--        (long "tck" <> help "subsampled population size" <> short 'k' <> value 0)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg m) = do

    dsts <- if dstarg == ""
               then fromJust <$> goalReadDatasetsCSV prjnm expnm
               else return [Dataset dstarg]

    if take 4 expnm == "true"

        then if last expnm == 'n'

           then forM_ dsts $ \dst -> do

                    (k,n,cs) <- getFittedMixtureLikelihood expnm dst

                    let wghts :: [[Double]]
                        stcs :: [[Double]]
                        tcss :: [[[Double]]]
                        (wghts,stcs,tcss) = withNat1 n (withNat k (analyzeMixtureTuningCurves cs))

                    goalWriteAnalysis prjnm expnm ananm (Just dst) wghts
                    goalAppendAnalysis prjnm expnm ananm (Just dst) stcs
                    mapM_ (goalAppendAnalysis prjnm expnm ananm (Just dst)) tcss

           else forM_ dsts $ \dst -> do

                    (k,cs) <- getFittedIPLikelihood expnm dst

                    let tcss = withNat k (analyzeTuningCurves cs)

                    goalWriteAnalysis prjnm expnm ananm (Just dst) tcss

        else if last expnm == 'n'

           then forM_ dsts $ \dst -> do

                    (k,(zxs :: [([Int], Double)])) <- getNeuralData expnm dst

                    stctcss <- realize $ withNat1 m (withNat k (fitAnalyzeMixtureTuningCurves zxs))

                    let wghts :: [[Double]]
                        stcs :: [[Double]]
                        tcss :: [[[Double]]]
                        (wghts,stcs,tcss) = stctcss

                    goalWriteAnalysis prjnm expnm ananm (Just dst) wghts
                    goalAppendAnalysis prjnm expnm ananm (Just dst) stcs
                    mapM_ (goalAppendAnalysis prjnm expnm ananm (Just dst)) tcss

           else forM_ dsts $ \dst -> do

                    (k,(zxs :: [([Int], Double)])) <- getNeuralData expnm dst

                    tcss <- realize $ withNat k (fitAnalyzeTuningCurves zxs)

                    goalWriteAnalysis prjnm expnm ananm (Just dst) tcss



--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts


--- Graveyard ---


--analyzeMixtureTuningCurves
--    :: forall k . KnownNat k
--    => [([Int],Double)]
--    -> Proxy k
--    -> [[Double]]
--analyzeMixtureTuningCurves zxs0 _ =
--    let zxs :: [(Response k, Double)]
--        zxs = strengthenNeuralData zxs0
--        ppc = fitIPLikelihood zxs
--        (rho0,rprms) = regressRectificationParameters ppc xsmps
--        rcrv = rectificationCurve rho0 rprms xsmps
--        tcs = listCoordinates . dualTransition <$> ppc >$>* xsmps
--        mxtcs = maximum <$> tcs
----        (mupr,vrpr) = estimateMeanVariance . map (head . tail . listCoordinates . toSource)
----            . S.toList . toRows . snd $ splitAffine ppc
--        stdoustr = concat ["mupr: ", show mupr, "; sdpr: ", show $ sqrt vrpr]
--    in trace stdoustr $ zipWith (++) (L.transpose [xsmps,potential <$> ppc >$>* xsmps,rcrv,mxtcs]) tcs

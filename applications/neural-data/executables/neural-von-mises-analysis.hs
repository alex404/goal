{-# LANGUAGE TupleSections,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Storable as S

import qualified Data.Map as M


--- Globals ---


eps :: Double
eps = -0.2

nepchs :: Int
nepchs = 100

pltsmps :: [Double]
pltsmps = range 0 (2*pi) 500

--- Functions ---


-- Von Mises Statistics --

fitPPC
    :: forall k . KnownNat k
    => M.Map Double [Response k]
    -> ([Double],[Response k],[Mean #> Natural # Neurons k <* VonMises])
fitPPC xzmp =
    let sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        ppc0 = vonMisesPopulationEncoder True (Left 1) sps
        (xs,zs) = unzip $ concat [ (x,) <$> zs0 | (x,zs0) <- M.toList xzmp ]
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
     in (xs,zs,take nepchs $ vanillaGradientSequence backprop eps defaultAdamPursuit ppc0)

ppcStimulusDerivatives
    :: KnownNat k => Mean #> Natural # Neurons k <* VonMises
    -> SamplePoint VonMises
    -> S.Vector k Double
ppcStimulusDerivatives ppc x =
    let fxs = coordinates . dualTransition $ ppc >.> mx
        tcs = toRows . snd $ splitAffine ppc
     in S.zipWith zipper fxs tcs
    where mx = sufficientStatistic x
          (cx,sx) = S.toPair $ coordinates mx
          zipper fx (Point cs) =
              let (tht1,tht2) = S.toPair cs
               in fx*(cx * tht2 - sx * tht1)

fisherInformation
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> SamplePoint VonMises
    -> Double
fisherInformation ppc x =
    let fxs2' = S.map square $ ppcStimulusDerivatives ppc x
        fxs = coordinates . dualTransition $ ppc >.>* x
     in S.sum $ S.zipWith (/) fxs2' fxs

-- Under the assumption of a flat prior
numericalVonMisesPPCPosterior0
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
numericalVonMisesPPCPosterior0 lkl z =
    let tcs = snd $ splitAffine lkl
        nxs = z *<.< tcs
        uldns x = nxs <.> sufficientStatistic x - S.sum (coordinates . dualTransition $ lkl >.>* x)
        avg = integrate 1e-500 uldns 0 (2*pi) / (2*pi)
        udns x = exp $ uldns x - avg
        nrm = trace (show $ udns <$> range 0 (2*pi) 8) $ integrate 1e-5000 udns 0 (2*pi)
     in (/nrm) . udns

-- Under the assumption of a flat prior
approximateVonMisesPPCPosterior0
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
approximateVonMisesPPCPosterior0 lkl z =
    let tcs = snd $ splitAffine lkl
     in density (z *<.< tcs)


--- Plots ---


tuningCurveLayout
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Layout Double Double
tuningCurveLayout lkl = execEC $ do

    goalLayout
    radiansAbscissa

    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Rate"
    --layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tuningCurves pltsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ zip pltsmps $ sumOfTuningCurves lkl pltsmps ]


decodingLayout
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Layout Double Double
decodingLayout lkl z = execEC $ do

    goalLayout
    radiansAbscissa

    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Posterior Density"

    let pstnm = numericalVonMisesPPCPosterior0 lkl z
        pstht = approximateVonMisesPPCPosterior0 lkl z

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip pltsmps $ pstnm <$> pltsmps]
        plot_lines_title .= "Numerical Posterior"

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [zip pltsmps $ pstht <$> pltsmps]
        plot_lines_title .= "Approximate Posterior"

generateDecodingLayouts
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> [Double]
    -> Random r [Layout Double Double]
generateDecodingLayouts lkl xs = do
    zs <- mapM samplePoint $ lkl >$>* xs
    return $ decodingLayout lkl <$> zs

--responseStatistics
--    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
--    => M.Map s [Response k]
--    -> Int
--    -> Proxy k'
--    -> Random r ((Double,Double,Double,Double),(Double,Double,Double))
--responseStatistics zxmp n _ = do
--    let mzxmp = empiricalTuningCurves zxmp
--    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
--    let subrs = subSampleResponses zxmp <$> idxss
--        subtcss = subSampleTuningCurves mzxmp . G.convert <$> idxss
--        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
--        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcss
--        tccvs = zip subtcss scvs
--        (mntcs,mncv) = minimumBy (comparing snd) tccvs
--        (mxtcs,mxcv) = maximumBy (comparing snd) tccvs
--    mis <- mapM (estimateEmpiricalPPCMutualInformation0 1) subtcss
--    mnmi <- estimateEmpiricalPPCMutualInformation0 nsmps mntcs
--    mxmi <- estimateEmpiricalPPCMutualInformation0 nsmps mxtcs
--    return ((average rcvs, mncv, average scvs, mxcv),(mnmi, average mis, mxmi))
--

--- Main ---

vonMisesAnalysis
    :: forall k
    . (KnownNat k, 1 <= k)
    => NeuralData k Double
    -> IO ()
vonMisesAnalysis nd = do

    let sgdpth = neuralDataProject nd ++ "/sgd-analysis"
        k = natValInt (Proxy :: Proxy k)
        vmpth = neuralDataProject nd ++ "/vm-analysis" ++ "/" ++ neuralDataOutputFile nd
            ++ "/" ++ show k ++ "-neurons"

    zxs <- getNeuralData nd

    let zxmp = head $ stimulusResponseMap <$> zxs

    idxs <- realize $ generateIndices (Proxy :: Proxy k)
    let zxmp' = subSampleResponses zxmp idxs

    let ppcs :: [Mean #> Natural # Neurons k <* VonMises]
        (xs,zs,ppcs) = fitPPC zxmp'

    let nlllyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "Epochs"
            layout_y_axis . laxis_title .= "Negative Log-Likelihood"

            plot . liftEC $ do

                plot_lines_title .= "SGD"
                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int)..]
                    $ stochasticConditionalCrossEntropy xs zs <$> ppcs ]

    let sgdflnm = "cross-entropy-descent-" ++ neuralDataOutputFile nd ++ "-" ++ show k ++ "-neurons"
    goalRenderableToPNG sgdpth sgdflnm 600 300 $ toRenderable nlllyt

    let xs0 = M.keys zxmp
        ppc1 = last ppcs

    dcdlyts <- realize $ generateDecodingLayouts ppc1 xs0

    goalRenderableToPNG vmpth "tuning-curves" 600 300 . toRenderable $ tuningCurveLayout ppc1

    sequence_ $ do
        (x,dcdlyt) <- zip xs0 dcdlyts
        let dcdflnm = "decoding-stimulus-" ++ show (roundSD 2 x)
        return . goalRenderableToPNG vmpth dcdflnm 600 300 $ toRenderable dcdlyt


main :: IO ()
main = do

    vonMisesAnalysis pattersonSmallPooled
    vonMisesAnalysis patterson112l44
    vonMisesAnalysis patterson112l45
    vonMisesAnalysis patterson112r35
    vonMisesAnalysis patterson112r36
    vonMisesAnalysis patterson105r62
    vonMisesAnalysis patterson107l114
    vonMisesAnalysis patterson112l16
    vonMisesAnalysis patterson112r32

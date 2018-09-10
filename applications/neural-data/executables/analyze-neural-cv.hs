{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import qualified Data.List as L


--- Globals ---

nsmps :: Int
nsmps = 1

--- Functions ---

-- Under the assumption of a flat prior
estimateEmpiricalPPCMutualInformation0
    :: (Ord s, KnownNat k)
    => Int
    -> M.Map s (Mean # Neurons k)
    -> Random r Double
estimateEmpiricalPPCMutualInformation0 nsmps' xzmp = do
    zss <- mapM (sample nsmps') xzmp
    let scl = fromIntegral . length $ M.keys xzmp
        f z = sum [ d * log (d*scl) | d <- M.elems $ empiricalPPCPosterior0 True xzmp z, d > 0 ]
    return . average $ f <$> concat zss


tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)

responseStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => M.Map s [Response k]
    -> Int
    -> Proxy k'
    -> Random r ((Double,Double,Double,Double),(Double,Double,Double))
responseStatistics zxmp n _ = do
    let mzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleTuningCurves mzxmp . G.convert <$> idxss
        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcss
        tccvs = zip subtcss scvs
        (mntcs,mncv) = minimumBy (comparing snd) tccvs
        (mxtcs,mxcv) = maximumBy (comparing snd) tccvs
    mis <- mapM (estimateEmpiricalPPCMutualInformation0 1) subtcss
    mnmi <- estimateEmpiricalPPCMutualInformation0 nsmps mntcs
    mxmi <- estimateEmpiricalPPCMutualInformation0 nsmps mxtcs
    return ((average rcvs, mncv, average scvs, mxcv),(mnmi, average mis, mxmi))


--- Plots ---


subSampleLayout :: [Int] -> [Double] -> [Double] -> [Double] -> [Double] -> Layout Int Double
subSampleLayout ks rcvs mncvs avcvs mxcvs = execEC $ do

    goalLayout

    --layout_x_axis . laxis_generate .= autoScaledLogAxis def
    --layout_y_axis . laxis_generate .= autoScaledLogAxis def
    layout_x_axis . laxis_title .= "Number of Neurons"
    layout_y_axis . laxis_title .= "Coefficient of Variation"
    --layout_y_axis . laxis_generate .= scaledAxis def (0,2)

    --let ks' = fromIntegral <$> ks
    let ks' = ks

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [zip ks' avcvs]
        plot_lines_title .= "Average CV"

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip ks' rcvs]
        plot_lines_title .= "Response CV"

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque grey)
        plot_lines_values .= [zip ks' mncvs, zip ks' mxcvs]
        plot_lines_title .= "Extremal CVs"

    plot . return $ hlinePlot "" (solidLine 3 $ opaque white) 0

mutualInformationLayout :: [Int] -> [Double] -> [Double] -> [Double] -> Layout Int Double
mutualInformationLayout ks mnmis avmis mxmis = execEC $ do

    goalLayout

    layout_x_axis . laxis_title .= "Number of Neurons"
    layout_y_axis . laxis_title .= "Mutual Information"

    --let ks' = fromIntegral <$> ks
    let ks' = ks

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [zip ks' avmis]
        plot_lines_title .= "Average MI"

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip ks' mxmis]
        plot_lines_title .= "Maximal-CV MI"

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [zip ks' mnmis]
        plot_lines_title .= "Minimal-CV MI"

    plot . return $ hlinePlot "" (solidLine 3 $ opaque white) 0


--- Main ---


analyzeCoefficientOfVariation
    :: forall k
    . (KnownNat k, 1 <= k, 5 <= k)
    => NeuralData k Double
    -> IO ()
analyzeCoefficientOfVariation nd = do

    let cvpth = neuralDataProject nd ++ "/cv-analysis"
        mipth = neuralDataProject nd ++ "/mi-analysis"
        dcdpth = neuralDataProject nd ++ "/decoder-analysis"
        --artpth = "signatures-of-bayesian-inference"

    zxs <- getNeuralData nd

    let zxmp = head $ stimulusResponseMap <$> zxs
        k = natValInt (Proxy :: Proxy k)
        ks = [1..k]

    let prxdcd = Proxy :: Proxy 5
        dcdk' = natValInt prxdcd

    allcvmis <- realize (B.generatePM' $ responseStatistics zxmp nsmps)

    let allcvs :: B.Vector k (Double,Double,Double,Double)
        allmis :: B.Vector k (Double,Double,Double)
        (allcvs,allmis) = B.unzip allcvmis
        (rcvs,mncvs,avcvs,mxcvs) = L.unzip4 $ B.toList allcvs
        (mnmis,avmis,mxmis) = L.unzip3 $ B.toList allmis

    print avmis

    --print mnmis
    --print avmis
    --print mxmis

    goalRenderableToPNG cvpth (neuralDataOutputFile nd) 600 300 . toRenderable
        $ subSampleLayout ks rcvs mncvs avcvs mxcvs

    goalRenderableToPNG mipth (neuralDataOutputFile nd) 600 300 . toRenderable
        $ mutualInformationLayout ks mnmis avmis mxmis

    --goalRenderableToPDF artpth (neuralDataOutputFile nd) 600 300 . toRenderable
    --    . (layout_legend .~ Nothing) $ subSampleLayout ks rcvs mncvs avcvs mxcvs

main :: IO ()
main = do

    --analyzeCoefficientOfVariation pattersonSmallPooled
    --analyzeCoefficientOfVariation patterson112l44
    --analyzeCoefficientOfVariation patterson112l45
    --analyzeCoefficientOfVariation patterson112r35
    --analyzeCoefficientOfVariation patterson112r36
    --analyzeCoefficientOfVariation patterson105r62
    --analyzeCoefficientOfVariation patterson107l114
    --analyzeCoefficientOfVariation patterson112l16
    --analyzeCoefficientOfVariation patterson112r32
    analyzeCoefficientOfVariation coenCagli1
    --analyzeCoefficientOfVariation coenCagli2
    --analyzeCoefficientOfVariation coenCagli3
    --analyzeCoefficientOfVariation coenCagli4
    --analyzeCoefficientOfVariation coenCagli5
    --analyzeCoefficientOfVariation coenCagli6
    --analyzeCoefficientOfVariation coenCagli7
    --analyzeCoefficientOfVariation coenCagli8
    --analyzeCoefficientOfVariation coenCagli9
    --analyzeCoefficientOfVariation coenCagli10
    --analyzeCoefficientOfVariation coenCagliPooled

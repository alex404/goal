{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Extra.Solver #-}

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

--- Functions ---

subSampleTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subSampleTuningCurves nzxmp idxs =
     Point . flip S.backpermute idxs . coordinates <$> nzxmp

subSampleResponses
    :: (Ord s, KnownNat k)
    => M.Map s [Response k]
    -> B.Vector k' Int
    -> M.Map s [Response k']
subSampleResponses zxmp idxs =
     map (`B.backpermute` idxs) <$> zxmp

tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)

--generateIndices
--    :: forall k k' d r . (KnownNat k, KnownNat k')
--    => Proxy k
--    -> Proxy k'
--    -> Random r (B.Vector (Min k k') Int)
--generateIndices _ _ = do
--    let idxs :: B.Vector k Int
--        idxs = B.generate finiteInt
--    subsampleVector idxs
--
--responseStatistics
--    :: forall s k k' d r . (Ord s, KnownNat k, KnownNat k')
--    => M.Map s [Response k]
--    -> Int
--    -> Proxy (Min k' k)
--    -> Random r (Double,Double,Double,Double)
--responseStatistics zxmp n _ = do
--    let mzxmp = empiricalTuningCurves zxmp
--    (idxss :: [B.Vector k Int]) <- replicateM n $ generateIndices (Proxy :: Proxy (Min k' k))
--    let subrs = subSampleResponses zxmp <$> idxss
--        subtcs = subSampleTuningCurves mzxmp . G.convert <$> idxss
--        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
--        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcs
--        tccvs = zip subtcs scvs
--        (_,mncv) = minimumBy (comparing snd) tccvs
--        (_,mxcv) = maximumBy (comparing snd) tccvs
--    return (average rcvs, mncv, average scvs, mxcv)
--

--- Plots ---


subSampleLayout :: [Int] -> [Double] -> [Double] -> [Double] -> [Double] -> Layout Int Double
subSampleLayout ks rcvs mncvs avcvs mxcvs = execEC $ do

    goalLayout

    --layout_x_axis . laxis_generate .= autoScaledLogAxis def
    --layout_y_axis . laxis_generate .= autoScaledLogAxis def
    layout_x_axis . laxis_title .= "Number of Neurons"
    layout_y_axis . laxis_title .= "Coefficient of Variation"
    --layout_y_axis . laxis_generate .= scaledAxis def (0,5)

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


--- Main ---


analyzeCoefficientOfVariation
    :: forall k s
    . (KnownNat k, 1 <= k, Read s, Ord s)
    => NeuralData k s
    -> IO ()
analyzeCoefficientOfVariation nd = do

    let ndpth = neuralDataProject nd ++ "/cv-analysis"
        artpth = "signatures-of-bayesian-inference"

    zxs <- getNeuralData nd

    let zxmp = head $ stimulusResponseMap <$> zxs
        k = natValInt (Proxy :: Proxy k)
        ks = [1..k]

--    let gpm :: (KnownNat k, KnownNat k', (k'+d) ~ k)
--            => Proxy k' -> Random r (Double,Double,Double,Double)
--        gpm = responseStatistics zxmp 10
--
--    (allcvs :: B.Vector k (Double,Double,Double,Double)) <- realize $ B.generatePM gpm
--
    --let (rcvs,mncvs,avcvs,mxcvs) = L.unzip4 allcvs

    --goalRenderableToPNG ndpth (neuralDataOutputFile nd) 600 300 . toRenderable
    --    $ subSampleLayout ks rcvs mncvs avcvs mxcvs

    --goalRenderableToPDF artpth (neuralDataOutputFile nd) 600 300 . toRenderable
    --    . (layout_legend .~ Nothing) $ subSampleLayout ks rcvs mncvs avcvs mxcvs
    return ()

main :: IO ()
main = undefined

    --analyzeCoefficientOfVariation pattersonSmallPooled
    --analyzeCoefficientOfVariation patterson112l44
    --analyzeCoefficientOfVariation patterson112l45
    --analyzeCoefficientOfVariation patterson112r35
    --analyzeCoefficientOfVariation patterson112r36
    --analyzeCoefficientOfVariation patterson105r62
    --analyzeCoefficientOfVariation patterson107l114
    --analyzeCoefficientOfVariation patterson112l16
    --analyzeCoefficientOfVariation patterson112r32
    --analyzeCoefficientOfVariation coenCagli1
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

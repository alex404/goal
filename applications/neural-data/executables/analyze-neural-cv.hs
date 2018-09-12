{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M


--- CSV ---


data Coefficients = Coefficients
    { responseCV :: Double
    , minimumCV :: Double
    , averageCV :: Double
    , maximumCV :: Double }
    deriving Generic

instance FromNamedRecord Coefficients
instance ToNamedRecord Coefficients
instance DefaultOrdered Coefficients
instance NFData Coefficients


--- Globals ---

nsmps :: Int
nsmps = 100

--- Functions ---


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
    -> Random r Coefficients
responseStatistics zxmp n _ = do
    let mzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleTuningCurves mzxmp . G.convert <$> idxss
        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcss
    return $ Coefficients (average rcvs) (minimum scvs) (average scvs) (maximum scvs)


--- Main ---


analyzeCoefficientOfVariation
    :: forall k
    . (KnownNat k, 1 <= k, 5 <= k)
    => NeuralData k Double
    -> IO ()
analyzeCoefficientOfVariation nd = do

    let cvpth = neuralDataProject nd ++ "/analysis/cv"
        --artpth = "signatures-of-bayesian-inference"

    zxs <- getNeuralData nd

    let zxmp = head $ stimulusResponseMap <$> zxs

    (allcvs :: B.Vector k Coefficients) <- realize (B.generatePM' $ responseStatistics zxmp nsmps)

    goalWriteCSV cvpth (neuralDataOutputFile nd) $ B.toList allcvs

main :: IO ()
main = do

    analyzeCoefficientOfVariation pattersonSmallPooled
    analyzeCoefficientOfVariation patterson112l44
    analyzeCoefficientOfVariation patterson112l45
    analyzeCoefficientOfVariation patterson112r35
    analyzeCoefficientOfVariation patterson112r36
    analyzeCoefficientOfVariation patterson105r62
    analyzeCoefficientOfVariation patterson107l114
    analyzeCoefficientOfVariation patterson112l16
    analyzeCoefficientOfVariation patterson112r32
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



--- Plots Graveyard ---


--goalRenderableToPNG cvpth (neuralDataOutputFile nd) 600 300 . toRenderable
--    $ subSampleLayout ks rcvs mncvs avcvs mxcvs
--
--goalRenderableToPNG mipth (neuralDataOutputFile nd) 600 300 . toRenderable
--    $ mutualInformationLayout ks mnmis avmis mxmis



--subSampleLayout :: [Int] -> [Double] -> [Double] -> [Double] -> [Double] -> Layout Int Double
--subSampleLayout ks rcvs mncvs avcvs mxcvs = execEC $ do
--
--    goalLayout
--
--    --layout_x_axis . laxis_generate .= autoScaledLogAxis def
--    --layout_y_axis . laxis_generate .= autoScaledLogAxis def
--    layout_x_axis . laxis_title .= "Number of Neurons"
--    layout_y_axis . laxis_title .= "Coefficient of Variation"
--    --layout_y_axis . laxis_generate .= scaledAxis def (0,2)
--
--    --let ks' = fromIntegral <$> ks
--    let ks' = ks
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque black)
--        plot_lines_values .= [zip ks' avcvs]
--        plot_lines_title .= "Average CV"
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque red)
--        plot_lines_values .= [zip ks' rcvs]
--        plot_lines_title .= "Response CV"
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque grey)
--        plot_lines_values .= [zip ks' mncvs, zip ks' mxcvs]
--        plot_lines_title .= "Extremal CVs"
--
--    plot . return $ hlinePlot "" (solidLine 3 $ opaque white) 0
--
--mutualInformationLayout :: [Int] -> [Double] -> [Double] -> [Double] -> Layout Int Double
--mutualInformationLayout ks mnmis avmis mxmis = execEC $ do
--
--    goalLayout
--
--    layout_x_axis . laxis_title .= "Number of Neurons"
--    layout_y_axis . laxis_title .= "Mutual Information"
--
--    --let ks' = fromIntegral <$> ks
--    let ks' = ks
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque black)
--        plot_lines_values .= [zip ks' avmis]
--        plot_lines_title .= "Average MI"
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque red)
--        plot_lines_values .= [zip ks' mxmis]
--        plot_lines_title .= "Maximal-CV MI"
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque blue)
--        plot_lines_values .= [zip ks' mnmis]
--        plot_lines_title .= "Minimal-CV MI"
--
--    plot . return $ hlinePlot "" (solidLine 3 $ opaque white) 0
--
--
--

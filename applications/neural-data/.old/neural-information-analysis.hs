{-# LANGUAGE TupleSections,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Storable as S

import qualified Data.Map as M
import qualified Data.List as L


--- Globals ---


--- Functions ---



--- Plots ---


tuningCurveLayout
    :: KnownNat k
    => M.Map s (Mean # Neurons k)
    -> LayoutLR Int Double Double
tuningCurveLayout lkl = execEC $ do

    goalLayoutLR

    layoutlr_x_axis . laxis_title .= "Stimulus"
    layoutlr_left_axis . laxis_title .= "Rate"
    layoutlr_right_axis . laxis_title .= "Sum"

    let tcs = L.transpose [(s,) <$> listCoordinates tcx | (s,tcx) <-  zip [0..] $ toList lkl]

    plotLeft . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= tcs

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ zip [0..] $ S.sum . coordinates <$>  toList lkl ]

    plotRight . return $ hlinePlot "" (solidLine 3 $ opaque white) 0


decodingLayout
    :: (Ord s, KnownNat k, Show s)
    => M.Map s (Mean # Neurons k)
    -> Response k
    -> Layout Int Double
decodingLayout lkl z = execEC $ do

    goalLayout

    let xs = M.keys lkl


    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Posterior Density"

    let pstnm x = empiricalPPCPosterior0 True lkl z M.! x
        pstht x = empiricalPPCPosterior0 False lkl z M.! x

    plot . liftEC $ do

        plot_points_style .= hollowCircles 8 4 (opaque red)
        plot_points_values .= (zip [0..] . traceGiven $ pstnm <$> xs)
        plot_points_title .= "Empirical Posterior"

    plot . liftEC $ do

        plot_points_style .= filledCircles 6 (opaque blue)
        plot_points_values .= (zip [0..] . traceGiven $ pstht <$> xs)
        plot_points_title .= "Approximate Posterior"

generateDecodingLayouts
    :: (KnownNat k, Ord s, Show s)
    => M.Map s (Mean # Neurons k)
    -> Random r [Layout Int Double]
generateDecodingLayouts lkl = do
    let xs = M.keys lkl
    zs <- mapM samplePoint $ (lkl M.!) <$> xs
    return $ decodingLayout lkl <$> zs

--- Main ---

informationAnalysis
    :: forall k s . (KnownNat k, 1 <= k, Ord s, Show s, Read s)
    => NeuralData k s
    -> IO ()
informationAnalysis nd = do

    let k = natValInt (Proxy :: Proxy k)
        dcdpth = neuralDataProject nd ++ "/empirical-analysis" ++ "/" ++ neuralDataOutputFile nd
            ++ "/" ++ show k ++ "-neurons"

    zxs <- getNeuralData nd

    let zxmp = head $ stimulusResponseMap <$> zxs

    idxs <- realize $ generateIndices (Proxy :: Proxy k)
    let zxmp' :: M.Map s [Response k]
        zxmp' = subSampleResponses zxmp idxs
        ppc = empiricalTuningCurves zxmp'

    let xs0 = M.keys zxmp

    dcdlyts <- realize $ generateDecodingLayouts ppc

    goalRenderableToPNG dcdpth "tuning-curves" 600 300 . toRenderable $ tuningCurveLayout ppc

    sequence_ $ do
        (x,dcdlyt) <- zip xs0 dcdlyts
        let dcdflnm = "decoding-stimulus-" ++ show x
        return . goalRenderableToPNG dcdpth dcdflnm 600 300 $ toRenderable dcdlyt


main :: IO ()
main = do

    informationAnalysis pattersonSmallPooled
    informationAnalysis patterson112l44
    informationAnalysis patterson112l45
    informationAnalysis patterson112r35
    informationAnalysis patterson112r36
    --informationAnalysis patterson105r62
    --informationAnalysis patterson107l114
    --informationAnalysis patterson112l16
    --informationAnalysis patterson112r32

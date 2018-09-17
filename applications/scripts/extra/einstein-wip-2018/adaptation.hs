{-# LANGUAGE TypeOperators,DataKinds,FlexibleContexts,TypeFamilies #-}

--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import Goal.Geometry
import Goal.Probability
import Goal.NeuralCircuits

--- Program ---


-- Globals --

-- Stimulus
x0 :: Double
x0 = pi

sx0 :: Source # VonMises
sx0 = Point $ S.doubleton pi 0

nx0 :: Natural # VonMises
nx0 = transition sx0

mix1,mix2,mix3 :: Source # VonMises
mix1 = Point $ S.doubleton (pi/2) 6
mix2 = Point $ S.doubleton pi 8
mix3 = Point $ S.doubleton (3*pi/2) 8
--
--mixt :: Natural # (Categorical Int 2 <*> VonMises)
--mixt =
--    let lb = Point $ S.replicate 0
--        ob = toNatural mix2
--        mtx = fromMatrix . S.fromRows . S.singleton . coordinates $ toNatural mix1 <-> ob
--         in joinHarmonium lb ob mtx
--
--unnormalizedMixtureDensity :: Natural # (Categorical Int 2 <*> VonMises) -> Double -> Double
--unnormalizedMixtureDensity hrm x =
--    let (lb,mx2,mtx) = splitHarmonium hrm
--        mx1 = Point (coordinates mtx) <+> mx2
--        rprms = snd $ categoricalHarmoniumRectificationParameters hrm
--        mxr = S.head . coordinates . toMean $ rprms <+> lb
--     in mxr * unnormalizedDensity mx1 x + (1 - mxr) * unnormalizedDensity mx2 x

-- Tuning Curves
type NNeurons = 10
type NNeurons' = 11

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

mus0 :: S.Vector NNeurons' Double
mus0 = S.range mnx mxx

mus :: S.Vector NNeurons Double
mus = fst $ S.splitAt mus0

--prc :: Double
--prc = 2
--
--sps :: Vector NNeurons (Source # VonMises)
--sps = Point . flip S.doubleton prc <$> mus

prcs :: S.Vector NNeurons Double
prcs = fromJust $ S.fromList [1,4,3,2,4,6,2,3,7,2]

sps :: S.Vector NNeurons (Source # VonMises)
sps = S.map Point $ S.zipWith S.doubleton mus prcs

lkl0 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl0 = vonMisesPopulationEncoder sps 1

dec :: Mean ~> Natural # VonMises <* R NNeurons Poisson
dec = joinAffine (transition sx0) . transpose . snd $ splitAffine lkl1

-- Rectification
rprms1,rprms2,rprms3 :: Natural # VonMises
rprms1 = Point $ S.doubleton 0 0
rprms2 = Point $ S.doubleton 1 0
rprms3 = Point $ S.doubleton 2 0

rho0 :: Double
rho0 = 5

rectification :: Natural # VonMises -> S.Vector NNeurons Double
rectification rprms = S.map (\x -> rprms <.> sufficientStatistic x + rho0) mus

-- Linear Least Squares

indpnds :: S.Vector NNeurons (S.Vector NNeurons Double)
indpnds = S.map (coordinates . dualTransition) $ lkl0 >$>* G.convert mus

dpnds1,dpnds2,dpnds3 :: S.Vector NNeurons Double
dpnds1 = rectification rprms1
dpnds2 = rectification rprms2
dpnds3 = rectification rprms3

gns1,gns2,gns3 :: S.Vector NNeurons Double
gns1 = linearLeastSquares indpnds dpnds1
gns2 = linearLeastSquares indpnds dpnds2
gns3 = linearLeastSquares indpnds dpnds3

lkl1,lkl2,lkl3 :: Mean ~> Natural # R NNeurons Poisson <* VonMises
lkl1 = vonMisesPopulationEncoder' sps gns1
lkl2 = vonMisesPopulationEncoder' sps gns2
lkl3 = vonMisesPopulationEncoder' sps gns3

-- Deconvolution

{-
f :: Double -> Double
f x = exp $ prc * cos x

g :: Natural # VonMises -> Double -> Double
g rprms x = sufficientStatistic x <.> rprms + rho0

gns1',gns2',gns3' :: S.Vector NNeurons Double
gns1' = deconvolve (g rprms1 <$> mus) (f <$> mus)
gns2' = deconvolve (g rprms2 <$> mus) (f <$> mus)
gns3' = deconvolve (g rprms3 <$> mus) (f <$> mus)
-}

-- Plot
xsmps :: B.Vector 500 Double
xsmps = B.range mnx mxx


-- Functions --

approximateInference :: [(Natural # VonMises, B.Vector NNeurons Int)] -> Natural # VonMises
approximateInference rprmzs =
    nx0 <+> foldl1 (<+>) [ dec >.>* z <-> rprm | (rprm,z) <- rprmzs ]

exactPosterior :: [B.Vector NNeurons Int] -> Double -> Double
exactPosterior zs y =
    let lkls x = [ density (lkl >.>* x) z  |  (lkl,z) <- zip [lkl1,lkl2,lkl3] zs]
        uposterior x = product $ unnormalizedDensity nx0 x : lkls x
        nrm = integrate 1e-5000 uposterior (-pi) pi
     in uposterior y / nrm

tclyt :: (Mean ~> Natural # R NNeurons Poisson <* VonMises) -> Natural # VonMises -> Layout Double Double
tclyt lkl rprms = execEC $ do

    goalLayout
    radiansAbscissa

    let mxplty = 10

    layout_x_axis . laxis_generate .= scaledAxis def (mnx,mxx)
    layout_y_axis . laxis_generate .= scaledAxis def (0,mxplty)
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Firing Rate"

    layout_plots
        .= [ vlinePlot "" (solidLine 3 $ opaque black) x0 ]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList . B.zip xsmps $ S.sum . coordinates . dualTransition <$> G.convert (lkl >$>* xsmps) ]

    plot . liftEC $ do
        plot_lines_style .= dashedLine 3 [8,10] (opaque blue)
        plot_lines_values .= [ toList . B.zip xsmps $ (\x -> rprms <.> sufficientStatistic x + rho0) <$> xsmps ]

prrlyt :: Layout Double Double
prrlyt = execEC $ do

    goalLayout
    radiansAbscissa

    let uprior x =
            0.3 * unnormalizedDensity (transition mix1) x
            + 0.5 * unnormalizedDensity (transition mix2) x
            + 0.2 * unnormalizedDensity (transition mix3) x
        nrm = integrate 1e-500 uprior (-pi) pi
        posterior x = uprior x / nrm

    let mx = 2

    layout_x_axis . laxis_title .= "Stimulus"

    layout_y_axis . laxis_generate .= scaledAxis def (0,mx)
    layout_y_axis . laxis_title .= "Belief"
    --layolr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [toList . B.zip xsmps $ posterior <$> xsmps]


byslyt :: [(Natural # VonMises, B.Vector NNeurons Int)] -> LayoutLR Double Int Double
byslyt rprmzs = execEC $ do

    goalLayoutLR
    radiansAbscissaLR

    let uposterior = unnormalizedDensity (approximateInference rprmzs)
        nrm = integrate 1e-2000 uposterior (-pi) pi
        posterior x = uposterior x / nrm

    let posterior0 = exactPosterior $ snd <$> rprmzs

    let mxlft = 10
        mxrght = 4

    layoutlr_x_axis . laxis_title .= "Stimulus/Preferred Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Posterior Density"
    --layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Spike Count"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    layoutlr_plots
        .= [ Left $ vlinePlot "" (solidLine 3 $ opaque black) x0 ]

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (withOpacity red 0.5)
        plot_lines_values .= [toList . B.zip xsmps $ posterior <$> xsmps]

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (withOpacity blue 0.5)
        plot_lines_values .= [toList . B.zip xsmps $ posterior0 <$> xsmps]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 4 (opaque black)
        plot_points_values .= zip (S.toList mus) (rtrns . toList . snd $ last rprmzs)

rtrns :: [a] -> [a]
rtrns rs = rs ++ [head rs]

-- Main --

main :: IO ()
main = do

    rs1 <- realize . sample $ lkl1 >.>* x0
    rs2 <- realize . sample $ lkl2 >.>* x0
    rs3 <- realize . sample $ lkl3 >.>* x0

    print gns1
    print gns2
    print gns3

    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves1" 300 200 . toRenderable $ tclyt lkl1 rprms1
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves2" 300 200 . toRenderable $ tclyt lkl2 rprms2
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves3" 300 200 . toRenderable $ tclyt lkl3 rprms3

    let rprmzs = [(rprms1,rs1),(rprms2,rs2),(rprms3,rs3)]

    void . goalRenderableToPDF "neural-circuits/adaptation" "mixture-prior" 500 200 $ toRenderable prrlyt

    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference1" 300 200 . toRenderable . byslyt $ take 1 rprmzs
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference2" 300 200 . toRenderable . byslyt $ take 2 rprmzs
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference3" 300 200 . toRenderable $ byslyt rprmzs

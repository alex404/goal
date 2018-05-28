{-# LANGUAGE TypeOperators,DataKinds,FlexibleContexts,TypeFamilies #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability

--- Program ---


-- Globals --

-- Stimulus
x0 :: Double
x0 = pi

sx0 :: Source # VonMises
sx0 = Point $ doubleton pi 0

nx0 :: Natural # VonMises
nx0 = transition sx0

-- Tuning Curves
type NNeurons = 10
type NNeurons' = 11

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

mus0 :: Vector NNeurons' Double
mus0 = rangeV mnx mxx

mus :: Vector NNeurons Double
mus = fst $ splitV mus0

--prc :: Double
--prc = 2
--
--sps :: Vector NNeurons (Source # VonMises)
--sps = Point . flip doubleton prc <$> mus

prcs :: Vector NNeurons Double
prcs = 1 & 4 & 3 & 2 & 4 & 6 & 2 & 3 & 7 & 2 & empty

sps :: Vector NNeurons (Source # VonMises)
sps = Point <$> zipWithV doubleton mus prcs

lkl0 :: Mean ~> Natural # R NNeurons Poisson <+ VonMises
lkl0 = vonMisesPopulationEncoder sps 1

dec :: Mean ~> Natural # VonMises <+ R NNeurons Poisson
dec = joinAffine (transition sx0) . transpose . snd $ splitAffine lkl1

-- Rectification
rprms1,rprms2,rprms3 :: Natural # VonMises
rprms1 = Point $ doubleton 0 0
rprms2 = Point $ doubleton 1 0
rprms3 = Point $ doubleton 2 0

rho0 :: Double
rho0 = 5

rectification :: Natural # VonMises -> Vector NNeurons Double
rectification rprms = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- Linear Least Squares

indpnds :: Vector NNeurons (Mean # R NNeurons Poisson)
indpnds = dualTransition <$> lkl0 >$>* mus

dpnds1,dpnds2,dpnds3 :: Vector NNeurons Double
dpnds1 = rectification rprms1
dpnds2 = rectification rprms2
dpnds3 = rectification rprms3

gns1,gns2,gns3 :: Vector NNeurons Double
gns1 = coordinates $ linearLeastSquares indpnds dpnds1
gns2 = coordinates $ linearLeastSquares indpnds dpnds2
gns3 = coordinates $ linearLeastSquares indpnds dpnds3

lkl1,lkl2,lkl3 :: Mean ~> Natural # R NNeurons Poisson <+ VonMises
lkl1 = vonMisesPopulationEncoder' sps gns1
lkl2 = vonMisesPopulationEncoder' sps gns2
lkl3 = vonMisesPopulationEncoder' sps gns3

-- Deconvolution

{-
f :: Double -> Double
f x = exp $ prc * cos x

g :: Natural # VonMises -> Double -> Double
g rprms x = sufficientStatistic x <.> rprms + rho0

gns1',gns2',gns3' :: Vector NNeurons Double
gns1' = deconvolve (g rprms1 <$> mus) (f <$> mus)
gns2' = deconvolve (g rprms2 <$> mus) (f <$> mus)
gns3' = deconvolve (g rprms3 <$> mus) (f <$> mus)
-}

-- Plot
xsmps :: Vector 500 Double
xsmps = rangeV mnx mxx

axs :: [Double]
axs = [0,pi/2,pi,3/2*pi,2*pi]


-- Functions --

approximateRho0 :: Natural # VonMises -> (Mean ~> Natural # R NNeurons Poisson <+ VonMises) -> Double
approximateRho0 rprms lkl =
    let lhs = sum . dualTransition <$> lkl >$>* xsmps
        rhs = (\x -> rprms <.> sufficientStatistic x) <$> xsmps
     in average $ zipWithV (-) lhs rhs

approximateInference :: [(Natural # VonMises,Vector NNeurons Int)] -> Natural # VonMises
approximateInference rprmzs =
    nx0 <+> foldl1 (<+>) [ dec >.>* z <-> rprm | (rprm,z) <- rprmzs ]

exactPosterior :: [Vector NNeurons Int] -> Double -> Double
exactPosterior zs y =
    let lkls x = [ density (lkl >.>* x) z  |  (lkl,z) <- zip [lkl1,lkl2,lkl3] zs]
        uposterior x = product $ unnormalizedDensity nx0 x : lkls x
        nrm = integrate 1e-500 uposterior (-pi) pi
     in uposterior y / nrm

tclyt :: (Mean ~> Natural # R NNeurons Poisson <+ VonMises) -> Natural # VonMises -> Layout Double Double
tclyt lkl rprms = execEC $ do

    goalLayout

    let mxplty = 10

    layout_x_axis . laxis_generate .= scaledAxis def (mnx,mxx)
    layout_y_axis . laxis_generate .= scaledAxis def (0,mxplty)
    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Firing Rate"

    layout_plots
        .= [ vlinePlot "" (solidLine 3 $ opaque black) x0 ]

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList . zipV xsmps $ sum . dualTransition <$> lkl >$>* xsmps ]

    plot . liftEC $ do
        plot_lines_style .= dashedLine 3 [8,10] (opaque red)
        plot_lines_values .= [ toList . zipV xsmps $ (\x -> rprms <.> sufficientStatistic x + rho0) <$> xsmps ]

byslyt :: [(Natural # VonMises,Vector NNeurons Int)] -> LayoutLR Double Int Double
byslyt rprmzs = execEC $ do

    goalLayoutLR

    let uposterior = unnormalizedDensity (approximateInference rprmzs)
        nrm = integrate 1e-500 uposterior (-pi) pi
        posterior x = uposterior x / nrm

    let posterior0 = exactPosterior $ snd <$> rprmzs

    let mxlft = 10
        mxrght = 4

    layoutlr_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]
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
        plot_lines_values .= [toList . zipV xsmps $ posterior <$> xsmps]

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (withOpacity blue 0.5)
        plot_lines_values .= [toList . zipV xsmps $ posterior0 <$> xsmps]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 4 (opaque black)
        plot_points_values .= zip (toList mus) (rtrns . toList . snd $ last rprmzs)

rtrns :: [a] -> [a]
rtrns rs = rs ++ [head rs]

-- Main --

main :: IO ()
main = do

    rs1 <- realize . generate $ lkl1 >.>* x0
    rs2 <- realize . generate $ lkl2 >.>* x0
    rs3 <- realize . generate $ lkl3 >.>* x0

    print gns1
    print gns2
    print gns3

    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves1" 250 200 . toRenderable $ tclyt lkl1 rprms1
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves2" 250 200 . toRenderable $ tclyt lkl2 rprms2
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-tuning-curves3" 250 200 . toRenderable $ tclyt lkl3 rprms3

    let rprmzs = [(rprms1,rs1),(rprms2,rs2),(rprms3,rs3)]

    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference1" 250 200 . toRenderable . byslyt $ take 1 rprmzs
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference2" 250 200 . toRenderable . byslyt $ take 2 rprmzs
    void . goalRenderableToPDF "neural-circuits/adaptation" "adaptive-inference3" 250 200 . toRenderable $ byslyt rprmzs

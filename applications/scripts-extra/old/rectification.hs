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

prc :: Double
prc = 2

prcs :: Vector NNeurons Double
prcs = (0.5*) <$> (1 & 4 & 3 & 2 & 4 & 6 & 2 & 3 & 7 & 2 & empty)

--prcs :: Vector NNeurons Double
--prcs = replicateV 1

sps :: Vector NNeurons (Source # VonMises)
sps = Point <$> zipWithV doubleton mus prcs

gns :: Vector NNeurons Double
gns = (0.1*) <$> (5 & 2 & 1 & 1 & 3 & 1 & 5 & 4 & 3 & 2 & empty)

--gns :: Vector NNeurons Double
--gns = replicateV 1

lkl :: Mean ~> Natural # R NNeurons Poisson <+ VonMises
lkl = vonMisesPopulationEncoder' sps gns

dec :: Mean ~> Natural # VonMises <+ R NNeurons Poisson
dec = joinAffine (transition sx0) . transpose . snd $ splitAffine lkl

-- Linear Least Squares

indpnds :: Vector NNeurons (Cartesian # Euclidean 3)
indpnds = (\mu -> Point $ 1 & cos mu & sin mu & empty) <$> mus

dpnds :: Vector NNeurons Double
dpnds = sum . dualTransition <$> lkl >$>* mus

rprm0s :: (Double, Natural # VonMises)
rprm0s =
    let (rho0',rhos') = headTail . coordinates $ linearLeastSquares indpnds dpnds
     in (rho0', Point rhos')

rho0 :: Double
rho0 = fst rprm0s

rprms :: Natural # VonMises
rprms = snd rprm0s

-- Plot
xsmps :: Vector 500 Double
xsmps = rangeV mnx mxx

axs :: [Double]
axs = [0,pi/2,pi,3/2*pi,2*pi]


-- Functions --

approximateInference :: [Vector NNeurons Int] -> Natural # VonMises
approximateInference zs =
    nx0 <+> foldl1 (<+>) [ dec >.>* z <-> rprms | z <- zs ]

exactPosterior :: [Vector NNeurons Int] -> Double -> Double
exactPosterior zs y =
    let lkls x = [ density (lkl >.>* x) z  | z <- zs]
        uposterior x = product $ unnormalizedDensity nx0 x : lkls x
        nrm = integrate 1e-5000 uposterior (-pi) pi
     in uposterior y / nrm

tclyt :: Layout Double Double
tclyt = execEC $ do

    goalLayout

    --let mxplty = 20

    layout_x_axis . laxis_generate .= scaledAxis def (mnx,mxx)
    --layout_y_axis . laxis_generate .= scaledAxis def (0,mxplty)
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

byslyt :: [Vector NNeurons Int] -> LayoutLR Double Int Double
byslyt zs = execEC $ do

    goalLayoutLR

    let uposterior = unnormalizedDensity (approximateInference zs)
        nrm = integrate 1e-500 uposterior (-pi) pi
        posterior x = uposterior x / nrm

    let posterior0 = exactPosterior zs

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
        plot_points_values .= zip (toList mus) (rtrns . toList $ last zs)

rtrns :: [a] -> [a]
rtrns rs = rs ++ [head rs]

-- Main --

main :: IO ()
main = do

    rs1 <- realize . generate $ lkl >.>* x0
    rs2 <- realize . generate $ lkl >.>* x0
    rs3 <- realize . generate $ lkl >.>* x0

    void . goalRenderableToPDF "neural-circuits/rectification" "rectified-tuning-curves" 250 200 . toRenderable $ tclyt

    let zs = [rs1,rs2,rs3]

    void . goalRenderableToPDF "neural-circuits/rectification" "rectified-inference1" 250 200 . toRenderable . byslyt $ take 1 zs
    void . goalRenderableToPDF "neural-circuits/rectification" "rectified-inference2" 250 200 . toRenderable . byslyt $ take 2 zs
    void . goalRenderableToPDF "neural-circuits/rectification" "rectified-inference3" 250 200 . toRenderable $ byslyt zs

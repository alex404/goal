{-# LANGUAGE TypeOperators,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability


--- Program ---

axs :: [Double]
axs = [pi,pi/2,0,-pi/2,-pi]

-- Globals --

mnx,mxx :: Double
mnx = -pi
mxx = pi

xsmps :: Vector 200 Double
xsmps = rangeV mnx mxx

mus1 :: Vector 5 Double
mus1 = rangeV mnx mxx

mus2 :: Vector 9 Double
mus2 = rangeV mnx mxx

kp :: Double
kp = 1.5

sps1 :: Vector 4 (Source # VonMises)
sps1 = tailV $ Point . flip doubleton kp <$> mus1

sps2 :: Vector 8 (Source # VonMises)
sps2 = tailV $ Point . flip doubleton kp <$> mus2

gn :: Double
gn = 1

lkl1 :: Mean ~> Natural # R 4 Poisson <+ VonMises
lkl1 = vonMisesPopulationEncoder sps1 gn

lkl2 :: Mean ~> Natural # R 8 Poisson <+ VonMises
lkl2 = vonMisesPopulationEncoder sps2 gn

-- Functions --

tcrnbl :: KnownNat k => Mean ~> Natural # R k Poisson <+ VonMises -> Layout Double Double
tcrnbl lkl = execEC $ do

    goalLayout

    let mxlft = 15

    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(pi,"π"),(pi/2,"π/2"),(0,"0"),(-pi/2,"-π/2"),(-pi,"-π")]
    layout_y_axis . laxis_generate .= scaledAxis def (0,mxlft)
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Rate"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList . zipV xsmps $ sum <$> coordinates . dualTransition <$> lkl >$>* xsmps ]

-- Main --

main :: IO ()
main = do

    void . goalRenderableToSVG "neural-circuits/tuning-curves-sum" "tuning-curves-sum1" 250 200 . toRenderable $ tcrnbl lkl1
    void . goalRenderableToSVG "neural-circuits/tuning-curves-sum" "tuning-curves-sum2" 250 200 . toRenderable $ tcrnbl lkl2

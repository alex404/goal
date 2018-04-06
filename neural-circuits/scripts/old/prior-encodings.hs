{-# LANGUAGE TypeOperators,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability


--- Program ---

axs :: [Double]
axs = [0,pi/2,pi,3*pi/2,2*pi]

-- Globals --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

xsmps :: Vector 200 Double
xsmps = rangeV mnx mxx

mu :: Double
mu = 3/4*pi

mus0 :: Vector 11 Double
mus0 = rangeV mnx mxx

mus :: Vector 10 Double
mus = fst $ splitV mus0

kp :: Double
kp = 2.5

sp :: Source # VonMises
sp = Point $ doubleton mu 6

sps :: Vector 10 (Source # VonMises)
sps = Point . flip doubleton kp <$> mus

gns1 :: Vector 10 Double
gns1 = 1.2 & 1.4 & 1.2 & 1 & 1 & 1 & 1.2 & 1.4 & 1.2 & singleton 1

lkl0 :: Mean ~> Natural # R 10 Poisson <+ VonMises
lkl0 = vonMisesPopulationEncoder sps 1

rprms :: Natural # VonMises
rprms = Point $ doubleton 4 0

rho0 :: Double
rho0 = 10

-- Linear System of Equations
tcmtx :: Natural ~> Cartesian # Euclidean 10 * R 10 Poisson
tcmtx = fromMatrix . fromRows $ coordinates . dualTransition <$> lkl0 >$>* mus

rectification :: Natural # VonMises -> Cartesian # Euclidean 10
rectification rprms = Point $ (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

gns2 :: Vector 10 Double
gns2 = coordinates $ fromJust (inverse tcmtx) >.> rectification rprms

lkl :: Mean ~> Natural # R 10 Poisson <+ VonMises
lkl = vonMisesPopulationEncoder' sps gns2

-- Functions --

tclyt :: Layout Double Double
tclyt = execEC $ do

    goalLayout

    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Rate"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= toList (toList <$> tuningCurves xsmps lkl)

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList . zipV xsmps $ sum <$> coordinates . dualTransition <$> lkl >$>* xsmps ]

vmlyt :: Layout Double Double
vmlyt = execEC $ do

    goalLayout

    let vm0 = unnormalizedDensity $ transition sp
        nrm = integrate 1e-5000 vm0 (-pi) pi
        vmdensity x = vm0 x / nrm

    layout_y_axis . laxis_generate .= scaledAxis def (0,1)
    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Density"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [ toList . zipV xsmps $ vmdensity <$> xsmps ]

prlyt :: Layout Double Double
prlyt = execEC $ do

    goalLayout

    let pr0 x = unnormalizedDensity (transition sp) x * exp (sum . coordinates . dualTransition $ lkl >.>* x)
        nrm = integrate 1e-500 pr0 (-pi) pi
        prdensity x = pr0 x / nrm

    layout_y_axis . laxis_generate .= scaledAxis def (0,1)
    layout_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_x_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]
    layout_x_axis . laxis_title .= "Stimulus"
    layout_y_axis . laxis_title .= "Density"

    plot . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [ toList . zipV xsmps $ prdensity <$> xsmps ]


-- Main --

main :: IO ()
main = do

    void . goalRenderableToPDF "neural-circuits/prior-encodings" "tuning-curves2" 200 200 $ toRenderable tclyt
    void . goalRenderableToPDF "neural-circuits/prior-encodings" "von-mises2" 200 200 $ toRenderable vmlyt
    void . goalRenderableToPDF "neural-circuits/prior-encodings" "prior2" 200 200 $ toRenderable prlyt

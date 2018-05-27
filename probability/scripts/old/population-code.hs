{-# LANGUAGE DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry
import Goal.Probability

--- Program ---


-- Globals --

x0 :: Double
x0 = 0

sp0 :: Source # VonMises
sp0 = Point $ doubleton 0 0.001

--- Program ---

-- Globals --

mnx,mxx :: Double
mnx = -pi
mxx = pi

xsmp :: Int
xsmp = 200

xsmps :: [Double]
xsmps = range mnx mxx xsmp

mus :: Vector 7 Double
mus = rangeV mnx mxx

kp :: Double
kp = 1.5

sps :: Vector 6 (Source # VonMises)
sps = tailV $ Point . flip doubleton kp <$> mus

gn1,gn2 :: Double
gn1 = 0.5
gn2 = 1
lkl1,lkl2 :: Mean ~> Natural # R 6 Poisson <* VonMises
lkl1 = vonMisesPopulationEncoder sps gn1
lkl2 = vonMisesPopulationEncoder sps gn2

dec :: Mean ~> Natural # VonMises <* R 6 Poisson
dec = joinAffine (transition sp0) . transpose . snd $ splitAffine lkl1

-- Functions --

rtrns :: [a] -> [a]
rtrns rs = last rs : rs

byslyt :: Vector 6 Int -> LayoutLR Double Int Double
byslyt rs = execEC $ do

    goalLayoutLR

    let uposterior = unnormalizedDensity (dec >.>* rs)
        nrm = integrate 1e-6 uposterior (-pi) pi
        posterior x = uposterior x / nrm


    let axs = [-pi,-pi/2,0,pi/2,pi]
        mxlft = 10
        mxrght = 2

    layoutlr_x_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(-pi,"-π"),(-pi/2,"-π/2"),(0,"0"),(pi/2,"π/2"),(pi,"π")]
    layoutlr_x_axis . laxis_title .= "Stimulus"

    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,mxrght)
    layoutlr_right_axis . laxis_title .= "Posterior Density"
    --layoutlr_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]

    layoutlr_left_axis . laxis_title .= "Response"
    layoutlr_left_axis . laxis_generate .= scaledIntAxis defaultIntAxis (0,mxlft)

    layoutlr_plots
        .= [ Left $ vlinePlot "" (solidLine 3 $ opaque black) x0 ]

    plotLeft . liftEC $ do
        plot_points_style .= filledCircles 3 (opaque black)
        plot_points_values .= zip (toList mus) (rtrns $ toList rs)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip xsmps $ posterior <$> xsmps]



-- Main --

main :: IO ()
main = do

    let rs0 = replicateV 0
    rs1 <- realize . generate $ lkl1 >.>* x0
    rs2 <- realize . generate $ lkl2 >.>* x0
    let rs3 = zipWithV (+) rs1 rs2

    let [rnbl0,rnbl1,rnbl2,rnbl3] = toRenderable . byslyt <$> [rs0,rs1,rs2,rs3]

    void $ goalRenderableToSVG "probability/population-code" "population-code-inference0" 250 150 rnbl0
    void $ goalRenderableToSVG "probability/population-code" "population-code-inference1" 250 150 rnbl1
    void $ goalRenderableToSVG "probability/population-code" "population-code-inference2" 250 150 rnbl2
    void $ goalRenderableToSVG "probability/population-code" "population-code-inference3" 250 150 rnbl3

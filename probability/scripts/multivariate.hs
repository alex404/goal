--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Data.Vector.Storable as C


--- Globals ---


nsmps = 10
tru = Source # fromList (MultivariateNormal 2) [0,0.5,1,0.5,0,1]

rng = (-4,4,200)
niso = 10

axprms = LinearAxisParams (show . round <$>) 5 5

vectorToPair xs = (xs C.! 0, xs C.! 1)
pairToVector (x,y) = C.fromList [x,y]

--- Main ---


main :: IO ()
main = do

    smps <- runWithSystemRandom . replicateM nsmps $ generate tru

    let mlenrm = Source # mle (MultivariateNormal 2) smps
        --efnrm = chart Natural $ transition mlenrm
        efnrm = transitionTo Natural $ Source # mle (MultivariateNormal 2) smps

        truf x y = density tru $ pairToVector (x,y)
        mlef x y = density mlenrm $ pairToVector (x,y)
        eff x y = density efnrm $ pairToVector (x,y)

        trucntrs = contours rng rng niso truf
        mlecntrs = contours rng rng niso mlef
        efcntrs = contours rng rng niso eff

        truclrs = rgbaGradient (1,0,0,0.5) (1,0,0,1) niso
        mleclrs = rgbaGradient  (0,0,1,0.5) (0,0,1,1) niso
        efclrs = rgbaGradient  (0,1,0,0.5) (0,1,0,1) niso

        bls = True : repeat False

        rnbl = toRenderable . execEC $ do

            layout_title .= ("Contours of True and Approximate Multivariate Normals" ++ "; Divergence from True: " ++ showFFloat (Just 3) (klDivergence mlenrm tru) "")

            layout_x_axis . laxis_generate .= scaledAxis axprms (-4,4)
            layout_x_axis . laxis_override .= axisGridHide
            layout_x_axis . laxis_title .= "x"
            layout_y_axis . laxis_generate .= scaledAxis axprms (-4,4)
            layout_y_axis . laxis_override .= axisGridHide
            layout_y_axis . laxis_title .= "y"

            sequence_ $ do

                ((_,cntr),clr,bl) <- zip3 trucntrs truclrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "True"
                    plot_lines_style .= solidLine 3 clr
                    plot_lines_values .= cntr

            sequence_ $ do

                ((_,cntr),clr,bl) <- zip3 mlecntrs mleclrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "MLE"
                    plot_lines_style .= dashedLine 3 [4,2,4] clr
                    plot_lines_values .= cntr

            sequence_ $ do

                ((_,cntr),clr,bl) <- zip3 efcntrs efclrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "EF-MLE"
                    plot_lines_style .= dashedLine 3 [2,4,2] clr
                    plot_lines_values .= cntr

            plot . liftEC $ do
                plot_points_title .= "Samples"
                plot_points_values .= map vectorToPair smps
                plot_points_style .= filledCircles 5 (opaque black)

    void $ goalRenderableToSVG "probability" "multivariate" 900 900 rnbl

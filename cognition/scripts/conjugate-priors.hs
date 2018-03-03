{-# LANGUAGE TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Cognition


--- Globals ---

zstp = 1
var = 1
x0 = 0.5
xs = range (-4) 4 1000

p = Standard # fromList Normal [x0,var]

--- Main ---


main :: IO ()
main = do

    zs <- runWithSystemRandom (replicateM 3 $ generate p)
    print zs
    let q0 = Standard # fromList Normal [0,1000]
        qs = scanl (flip $ kalmanInference1D zstp var) q0 zs
        clrs = rgbaGradient (0,0,1,1) (1,0,0,1) 4

    let lyt = execEC $ do

            goalLayout
            layout_x_axis . laxis_title .= "x"
            layout_y_axis . laxis_title .= "Probability Density"

            plot . liftCState . return $ vlinePlot "" (solidLine 3 (opaque black)) x0

            sequence_ $ do

                (q,clr) <- zip qs clrs

                return $ plot . liftEC $ do

                    plot_lines_style .= solidLine 3 clr
                    plot_lines_values .= [zip xs (density q <$> xs)]

    goalRenderableToPDF "cognition" "conjugate-priors" 250 200 $ toRenderable lyt

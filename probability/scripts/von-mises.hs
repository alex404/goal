{-# LANGUAGE DataKinds,ScopedTypeVariables,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

--- Program ---


-- Globals --

type SampleSize = 100000

mu,kap :: Double
mu = -2
kap = 2

tru :: Source # VonMises
tru = Point $ S.doubleton mu kap

-- Plot

mn,mx :: Double
(mn,mx) = (0,2*pi)

xs :: [Double]
xs = range mn mx 200

nb :: Int
nb = 25


-- Main --

main :: IO ()
main = do

    (smps :: Sample SampleSize VonMises) <- realize $ sample tru
    let cosht = average $ cos <$> smps
        sinht = average $ sin <$> smps

    let (cosht',sinht') = S.toPair . coordinates . dualTransition $ toNatural tru

    putStrLn "Expected Value of Cos (Samples):"
    print cosht
    putStrLn "Expected Value of Cos (Bessel Approx.):"
    print cosht'

    putStrLn "Expected Value of Sin (Samples):"
    print sinht
    putStrLn "Expected Value of Sin (Bessel Approx.):"
    print sinht'

    let lyt = execEC $ do

            histogramLayoutLR nb mn mx

            layoutlr_title .= "Sampling Algorithm vs . Unnormalized von Mises Density"
            layoutlr_left_axis . laxis_title .= "Sample Count"
            layoutlr_left_axis . laxis_override .= axisGridHide

            layoutlr_right_axis . laxis_title .= "Unnormalized Probability Mass"
            layoutlr_right_axis . laxis_override .= axisGridHide

            layoutlr_x_axis . laxis_title .= "Value"
            layoutlr_x_axis . laxis_override .= axisGridHide

            plotLeft . fmap plotBars . liftEC $ do
                void $ histogramPlot nb mn mx [B.toList smps]
                plot_bars_titles .= ["Samples"]
                plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

            plotRight . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_title .= "Density"
                plot_lines_values .= [ [(x,unnormalizedDensity (transition tru) x)  | x <- xs] ]

    goalRenderableToSVG "probability" "von-mises" 800 400 $ toRenderable lyt


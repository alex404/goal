{-# LANGUAGE DataKinds,ScopedTypeVariables,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

--- Program ---


-- Globals --

nsmps :: Int
nsmps = 100000

mu,kap :: Double
mu = 2
kap = 2

tru :: Source # VonMises
tru = Point $ S.doubleton mu kap

-- Plot

mn,mx :: Double
(mn,mx) = (0,2*pi)

xs :: [Double]
xs = range mn mx 200

nb :: Int
nb = 50

normalize :: [(Double, [Int])] -> [(Double,[Double])]
normalize xys = do
    (x,[y]) <- xys
    return (x,[fromIntegral y / nrm])
        where nrm = (2*pi / fromIntegral nb) * fromIntegral nsmps





-- Main --

main :: IO ()
main = do

    smps <- realize $ sample nsmps tru
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

    let (xys,_,_) = histogram nb mn mx [smps]
        xfs = normalize xys
        l2 = sqrt $ sum [square (f - density tru x) | (x,[f]) <- xfs]

    putStrLn "L2:"
    print l2

    let lyt = execEC $ do

            goalLayout
            radiansAbscissa

            layout_title .=  "Sampling Algorithm vs von Mises Density"

            plot . fmap plotBars . liftEC $ do
                plot_bars_values .= xfs
                plot_bars_titles .= ["Samples"]
                plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_title .= "Density"
                plot_lines_values .= [ [(x,density tru x)  | x <- xs] ]

    goalRenderableToSVG "probability" "von-mises" 800 400 $ toRenderable lyt


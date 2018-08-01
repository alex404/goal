{-# LANGUAGE DataKinds,TypeOperators #-}
--- Imports ---


import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

--- Globals ---


-- Functions --

f :: RealFrac x => B.Vector 2 x -> x
f xs =
    let (x,y) = G.toPair xs
        two :: Int
        two = 2
     in x^two + y^two + (x-y)^two

-- Plot --

niso :: Int
niso = 10

cntrf :: Double -> Double -> Double
cntrf x y = f $ G.doubleton x y

rng :: (Double,Double,Int)
rng = (-4,4,400)

clrs :: [AlphaColour Double]
clrs = rgbaGradient (0.9,0,0,1) (0,0,0,1) niso

-- Gradient Descent --

p0 :: Cartesian # Euclidean 2
p0 = Point $ G.doubleton (-4) 2

bnd,eps :: Double
bnd = 0.0001
eps = -0.05

cauchify :: [Cartesian # Euclidean 2] -> [Cartesian # Euclidean 2]
cauchify = cauchySequence euclideanDistance bnd

grds,mtms,adms :: [Cartesian # Euclidean 2]
grds = cauchify $ gradientSequence (differential' f) eps Classic p0
--nwts = cauchify $ newtonSequence (differential' f) p0
mtms = cauchify $ gradientSequence (differential' f) eps (defaultMomentumPursuit 0.9) p0
adms = cauchify $ gradientSequence (differential' f) eps defaultAdamPursuit p0


--- Main ---


main :: IO ()
main = do

    -- Contour plots
    let rnbl = toRenderable . execEC $ do

            layout_title .= "Gradient Descent on a Convex Function"
            layout_x_axis . laxis_title .= "x"
            layout_y_axis . laxis_title .= "y"


            let cntrs = contours rng rng niso cntrf

            sequence_ $ do

                ((_,cntr),clr) <- zip cntrs clrs

                return . plot . liftEC $ do

                    plot_lines_style .= solidLine 3 clr
                    plot_lines_values .= cntr

            plot . liftEC $ do
                plot_points_style .= filledCircles 5 (opaque red)
                plot_points_values .= [(0,0)]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_title .= "Gradient Sequence"
                plot_lines_values .= [G.toPair . coordinates <$> grds]

--            plot . liftEC $ do
--                plot_lines_style .= solidLine 3 (opaque blue)
--                plot_lines_title .= "Newton Method"
--                plot_lines_values .= [G.toPair . coordinates <$> nwts]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque green)
                plot_lines_title .= "Momentum"
                plot_lines_values .= [G.toPair . coordinates <$> mtms]

            plot . liftEC $ do
                plot_lines_style .= solidLine 3 (opaque purple)
                plot_lines_title .= "Adam"
                plot_lines_values .= [G.toPair . coordinates <$> adms]


    putStrLn "Gradient Descent Steps:"
    print $ length grds - 1

    putStrLn "Momentum Steps:"
    print $ length mtms - 1

    putStrLn "Adam Steps:"
    print $ length adms - 1

--    putStrLn "Newton Method Steps:"
--    print $ length nwts - 1

    void $ goalRenderableToSVG "geometry" "gradient-descent" 600 600 rnbl

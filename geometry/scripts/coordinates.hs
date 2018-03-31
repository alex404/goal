{-# LANGUAGE TypeOperators,DataKinds #-}
--- Imports ---


import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

--- Program ---


-- Globals --

mxx,mnx,mxy,mny :: Double
mxx = pi
mnx = -pi
mxy = pi
mny = -pi

nstps,npnts :: Int
nstps = 11
npnts = 50

hlns,vlns,lns0 :: [[S.Vector 2 Double]]
hlns = [ [ S.doubleton x y | x <- range mnx mxx npnts ] | y <- range mnx mxx nstps ]
vlns = [ [ S.doubleton x y | y <- range mny mxy npnts ] | x <- range mny mxy nstps ]
lns0 = hlns ++ vlns

eclds :: [[Cartesian # Euclidean 2]]
eclds = map Point <$> lns0

plrs :: [[Polar # Euclidean 2]]
plrs = map Point <$> lns0

layoutMaker :: [[Cartesian # Euclidean 2]] -> Layout Double Double
layoutMaker lns = execEC $ do

    goalLayout
    layout_x_axis . laxis_title .= "x"
    layout_y_axis . laxis_title .= "y"

    plot . liftEC $ do

        plot_lines_values .= (map (G.toPair . coordinates) <$> lns)
        plot_lines_style .= solidLine 3 (opaque black)


-- Main --

main :: IO ()
main = do
    let rnbl1 = toRenderable $ layoutMaker eclds
        rnbl2 = toRenderable . layoutMaker $ map transition <$> plrs
    goalRenderableToSVG "geometry" "cartesian-coordinates" 250 250 rnbl1
    goalRenderableToSVG "geometry" "polar-coordinates" 250 250 rnbl2

{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

nsmps :: Int
nsmps = 10

mux,muy,vrx,cvr,vry :: Double
mux = 1
muy = -1
vrx = 1
cvr = 0.7
vry = 2

tru :: Source # MultivariateNormal 2
tru = Point $ S.fromTuple (mux,muy,vrx,cvr,vry)

rng :: (Double,Double,Int)
rng = (-8,8,1000)
niso :: Int
niso = 10


--- Main ---


main :: IO ()
main = do

    smps <- realize $ sample nsmps tru

    let mnrm :: Mean # MultivariateNormal 2
        mnrm = sufficientStatisticT smps
        snrm = toSource mnrm
        nnrm = toNatural snrm

        truf x y = density tru $ S.fromTuple (x,y)
        sf x y = density snrm $ S.fromTuple (x,y)
        nf x y = density nnrm $ S.fromTuple (x,y)

        trucntrs = contours rng rng niso truf
        scntrs = contours rng rng niso sf
        ncntrs = contours rng rng niso nf

        bls = True : repeat False

        rnbl = toRenderable . execEC $ do

            goalLayout

            sequence_ $ do

                ((_,cntr),bl) <- zip trucntrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "True"
                    plot_lines_style .= dashedLine 3 [20,20] (opaque green)
                    plot_lines_values .= cntr

            sequence_ $ do

                ((_,cntr),bl) <- zip scntrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "Source"
                    plot_lines_style .= solidLine 6 (opaque blue)
                    plot_lines_values .= cntr

            sequence_ $ do

                ((_,cntr),bl) <- zip ncntrs bls

                return . plot . liftEC $ do

                    when bl $ plot_lines_title .= "Natural"
                    plot_lines_style .= solidLine 2 (opaque red)
                    plot_lines_values .= cntr

            plot . liftEC $ do
                plot_points_title .= "Samples"
                plot_points_values .= map S.toPair smps
                plot_points_style .= filledCircles 8 (opaque black)

    void $ goalRenderableToSVG "probability" "multivariate" 600 600 rnbl

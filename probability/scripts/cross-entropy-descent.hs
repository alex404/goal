{-# LANGUAGE TypeOperators, TypeFamilies, FlexibleContexts #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

--- Globals ---

-- True Normal --

sp1 :: Source # Normal
sp1 = Point $ S.doubleton 2 3

-- Gradient Descent --

meps,neps,geps :: Double
meps = -0.1
neps = -0.01
geps = -0.1

mbnd,nbnd,gbnd :: Double
mbnd = 1e-10
nbnd = 1e-10
gbnd = 1e-10

sp0 :: Source # Normal
sp0 = Point $ S.doubleton 0.5 1.5

np0 :: Natural # Normal
np0 = transition sp0

mp0 :: Mean # Normal
mp0 = transition sp0

-- Plot --

mnmu,mxmu,mnvr,mxvr :: Double
mnmu = 0
mxmu = 4
mnvr = 1
mxvr = 5

murng,vrrng :: (Double,Double,Int)
murng = (mnmu,mxmu,1000)
vrrng = (mnvr,mxvr,1000)

niso :: Int
niso = 20

clrs :: [AlphaColour Double]
clrs = rgbaGradient (0,0,0,0.8) (0,0,0,0.1) niso

naturalDifferentials
     :: Point Natural Normal
     -> CotangentPair Natural Normal
naturalDifferentials q = joinTangentPair q $ crossEntropyDifferential (transition sp1) q

-- Functions --

axprms :: LinearAxisParams Double
axprms = LinearAxisParams (show . round <$>) 4 4

-- Layout --

main :: IO ()
main = do

    --let mps = cauchySequence relativeEntropy mbnd $ vanillaGradientSequence meps meanDifferentials mp0
        --nmps = take stps $ gradientSequence meps mixtureDifferentials mp0

    let nps = cauchySequence (transition2 relativeEntropy) nbnd
            $ vanillaGradientSequence naturalDifferentials neps Classic np0
        gps = cauchySequence (transition2 relativeEntropy) gbnd
            $ gradientSequence naturalDifferentials geps Classic np0

    --putStrLn "Mean Coordinate Descent Steps:"
    --print $ length mps

    putStrLn "Natural Coordinate Descent Steps:"
    print $ length nps

    putStrLn "Geometric Gradient Descent Steps:"
    print $ length gps

    let rnbl = toRenderable . execEC $ do

            goalLayout

            let f mu vr =
                    let p :: Source # Normal
                        p = Point $ S.doubleton mu vr
                     in transition2 relativeEntropy sp1 p
                cntrs = contours murng vrrng niso f

            layout_x_axis . laxis_generate .= scaledAxis axprms (mnmu,mxmu)
            layout_x_axis . laxis_title .= "μ"
            --layout_x_axis . laxis_title_style .= (font_size .~ 14 $ def)
            layout_y_axis . laxis_generate .= scaledAxis axprms (mnvr,mxvr)
            layout_y_axis . laxis_title .= "σ²"
            --layout_y_axis . laxis_title_style .= (font_size .~ 14 $ def)

            sequence_ $ do

                ((_,cntr),clr) <- zip cntrs clrs

                return . plot . liftEC $ do

                    plot_lines_style .= solidLine 2 clr
                    plot_lines_values .= cntr

            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque red)
                plot_lines_values .= [S.toPair . coordinates . toSource <$> nps]

--            plot . liftEC $ do
--                plot_lines_style .= solidLine 2 (opaque blue)
--                plot_lines_values .= [S.toPair . coordinates . toSource <$> mps]

{-
            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque purple)
                plot_lines_values .= [B.toPair <$> nmps]
                -}

            plot . liftEC $ do
                plot_lines_style .= solidLine 2 (opaque purple)
                plot_lines_values .= [S.toPair . coordinates . toSource <$> gps]

            plot . liftEC $ do
                plot_points_style .= filledCircles 4 (opaque black)
                plot_points_values .= [S.toPair $ coordinates sp1]

    void $ goalRenderableToSVG "probability" "cross-entropy-descent" 500 300 rnbl

{-# LANGUAGE TypeOperators, DataKinds #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S


--- Globals ---

prjdr :: String
prjdr = "probability/expectation-maximization"

-- Manifolds --

type Latent = Categorical Int 3
type Observable = (Normal, Normal)
type Harmonium' = Harmonium Tensor Observable Latent

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 2
muy1 = 0.5
mux2 = -2
muy2 = 1
mux3 = 0.5
muy3 = -2

vrx1,vrx2,vrx3,vry1,vry2,vry3 :: Double
vrx1 = 1.5
vry1 = 0.5
vrx2 = 0.5
vry2 = 1.5
vrx3 = 1
vry3 = 0.5

nrm1,nrm2,nrm3 :: Source # Observable
nrm1 = Point $ S.fromTuple (mux1, vrx1, muy1, vry1)
nrm2 = Point $ S.fromTuple (mux2, vrx2, muy2, vry2)
nrm3 = Point $ S.fromTuple (mux3, vrx3, muy3, vry3)

nrms :: S.Vector 3 (Natural # Observable)
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Harmonium Tensor Observable Latent
truhrm = buildCategoricalHarmonium zero nrms $ toNatural wghts

-- Mixture Distributions --

nrm1',nrm2',nrm3' :: Source # Observable
nrm1' = Point $ S.fromTuple (1, 2, 0, 2)
nrm2' = Point $ S.fromTuple (0, 2, 0, 2)
nrm3' = Point $ S.fromTuple (-1, 2, 0, 2)

nrms' :: S.Vector 3 (Natural # Observable)
nrms' = S.fromTuple (toNatural nrm1',toNatural nrm2',toNatural nrm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Latent
wghts' = Point $ S.doubleton mix1' mix2'

hrm0 :: Natural # Harmonium Tensor Observable Latent
hrm0 = buildCategoricalHarmonium zero nrms' $ toNatural wghts'

-- Training --

w0 :: Source # Normal
w0 = Point $ S.doubleton 0 0.001

type SampleSize = 100

-- Functions --

filterCat :: [SamplePoint Harmonium'] -> Int -> [SamplePoint Observable]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

elipse :: Source # Observable -> [[(Double,Double)]]
elipse nrm =
    let [mux,vrx,muy,vry] = listCoordinates nrm
        f t = (mux + sqrt vrx * cos t, muy + sqrt vry * sin t)
     in [f <$> range 0 (2*pi) 100]

plotCluster
    :: AlphaColour Double
    -> Sample SampleSize Harmonium'
    -> Int
    -> Source # Observable
    -> EC (Layout Double Double) ()
plotCluster clr cxys cat nrm = do

    let [mux,_,muy,_] = listCoordinates nrm

    plot . liftEC $ do

        plot_points_values .= [(mux,muy)]
        plot_points_style .= filledCircles 7 clr

    plot . liftEC $ do

        plot_points_values .= filterCat (toList cxys) cat
        plot_points_style .= filledCircles 3 clr

    plot . liftEC $ do

        plot_lines_values .= elipse nrm
        plot_lines_style .= solidLine 3 clr


clusterLayout :: Natural # Harmonium' -> Sample SampleSize Harmonium' -> Layout Double Double
clusterLayout hrm1 cxys = execEC $ do

    let lkl = fst $ splitBottomHarmonium hrm1

        [mux1',_,muy1',_] = listCoordinates . toSource $ lkl >.>* 0
        [mux2',_,muy2',_] = listCoordinates . toSource $ lkl >.>* 1
        [mux3',_,muy3',_] = listCoordinates . toSource $ lkl >.>* 2

    goalLayout

    layout_x_axis . laxis_title .= "x"
    layout_y_axis . laxis_title .= "y"

    let def' = def {_la_nLabels = 3}
    layout_x_axis . laxis_generate .= scaledAxis def' (-5,5)
    layout_y_axis . laxis_generate .= scaledAxis def' (-5,5)

    plotCluster (opaque red) cxys 0 nrm1
    plotCluster (opaque green) cxys 1 nrm2
    plotCluster (opaque blue) cxys 2 nrm3

    plot . liftEC $ do

        plot_points_values .= [(mux1',muy1'),(mux2',muy2'),(mux3',muy3')]
        plot_points_style .= filledCircles 6 (opaque black)


--- Main ---


main :: IO ()
main = do

    vcxys <- realize $ sample truhrm
    tcxys <- realize $ sample truhrm

    let vxys,txys :: Sample SampleSize Observable
        vxys = hHead <$> vcxys
        txys = hHead <$> tcxys

    let hrms = take 20 $ iterate (categoricalHarmoniumExpectationMaximization txys) hrm0

    let anlls = [ average $ categoricalHarmoniumNegativeLogLikelihood hrm <$> vxys | hrm <- hrms ]

    let anllrnbl = toRenderable . execEC $ do

            goalLayout

            layout_x_axis . laxis_title .= "Epoch"
            layout_y_axis . laxis_title .= "-Log-Likelihood"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [zip [(0 :: Int)..] anlls]

    sequence_ $ do

        (hrm,k) <- zip hrms [(0 :: Int)..]
        return . goalRenderableToSVG prjdr ("em-step-" ++ show k) 800 800
            . toRenderable $ clusterLayout hrm vcxys

    goalRenderableToSVG prjdr "cross-entropy-descent" 1200 600 anllrnbl

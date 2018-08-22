{-# LANGUAGE GADTs,FlexibleContexts,TypeOperators,DataKinds #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

prjdr :: String
prjdr = "probability/gaussian-mixture"

-- Manifolds --

type Latent = Categorical Int 3
type Observable = MultivariateNormal 2
type Harmonium' = Harmonium Tensor Observable Latent

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 2
muy1 = 0.5
mux2 = -2
muy2 = 1
mux3 = 0.5
muy3 = -2

vrx1,vrx2,vrx3,cvr1,cvr2,cvr3,vry1,vry2,vry3 :: Double
vrx1 = 1.5
cvr1 = 0.5
vry1 = 0.5
vrx2 = 0.5
cvr2 = 0.2
vry2 = 1.5
vrx3 = 1
cvr3 = -0.2
vry3 = 0.5

nrm1,nrm2,nrm3 :: Source # Observable
nrm1 = Point $ S.fromTuple (mux1, muy1, vrx1, cvr1, vry1)
nrm2 = Point $ S.fromTuple (mux2, muy2, vrx2, cvr2, vry2)
nrm3 = Point $ S.fromTuple (mux3, muy3, vrx3, cvr3, vry3)

nrms :: S.Vector 3 (Natural # Observable)
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Harmonium Tensor Observable Latent
truhrm = buildMixtureModel nrms $ toNatural wghts

-- Mixture Distributions --

nrm1',nrm2',nrm3' :: Source # Observable
nrm1' = Point $ S.fromTuple (1, 0, 2, 0, 2)
nrm2' = Point $ S.fromTuple (0, 0, 2, 0, 2)
nrm3' = Point $ S.fromTuple (-1, 0, 2, 0, 2)

nrms' :: S.Vector 3 (Natural # Observable)
nrms' = S.fromTuple (toNatural nrm1',toNatural nrm2',toNatural nrm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Latent
wghts' = Point $ S.doubleton mix1' mix2'

hrm0 :: Natural # Harmonium Tensor Observable Latent
hrm0 = buildMixtureModel nrms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 100

-- EM
emepchs :: Int
emepchs = 50

-- Functions --

filterCat :: [SamplePoint Harmonium'] -> Int -> [SamplePoint Observable]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

elipse :: Source # MultivariateNormal 2 -> [[(Double,Double)]]
elipse nrm =
    let (mu,sgma) = splitMultivariateNormal nrm
        rtsgma = S.matrixRoot sgma
        f t = S.toPair . S.add mu . S.matrixVectorMultiply rtsgma $ S.fromTuple (cos t, sin t)
     in [f <$> range 0 (2*pi) 100]

plotNormal
    :: AlphaColour Double
    -> Source # Observable
    -> EC (Layout Double Double) ()
plotNormal clr nrm = do

    let mux:muy:_ = listCoordinates nrm

    plot . liftEC $ do

        plot_points_values .= [(mux,muy)]
        plot_points_style .= filledCircles 7 clr

    plot . liftEC $ do

        plot_lines_values .= elipse nrm
        plot_lines_style .= solidLine 3 clr

plotSample
    :: AlphaColour Double
    -> Sample Harmonium'
    -> Int
    -> EC (Layout Double Double) ()
plotSample clr cxys cat =

    plot . liftEC $ do

        plot_points_values .= (S.toPair <$> filterCat cxys cat)
        plot_points_style .= filledCircles 3 clr


clusterLayout :: Natural # Harmonium' -> Sample Harmonium' -> Layout Double Double
clusterLayout hrm1 cxys = execEC $ do

    let lkl = fst $ splitBottomHarmonium hrm1

        nrm1'' = toSource $ lkl >.>* 0
        nrm2'' = toSource $ lkl >.>* 1
        nrm3'' = toSource $ lkl >.>* 2

    goalLayout

    layout_x_axis . laxis_title .= "x"
    layout_y_axis . laxis_title .= "y"

    let def' = def {_la_nLabels = 3}
    layout_x_axis . laxis_generate .= scaledAxis def' (-5,5)
    layout_y_axis . laxis_generate .= scaledAxis def' (-5,5)

    plotNormal (opaque red) nrm1
    plotNormal (opaque green) nrm2
    plotNormal (opaque blue) nrm3
    plotSample (opaque red) cxys 0
    plotSample (opaque green) cxys 1
    plotSample (opaque blue) cxys 2

    plotNormal (opaque black) nrm1''
    plotNormal (opaque black) nrm2''
    plotNormal (opaque black) nrm3''


--- Main ---


main :: IO ()
main = do

    putStrLn "Covariance Matrices are SPD:"
    print $ S.isSemiPositiveDefinite . snd . splitMultivariateNormal <$> [nrm1,nrm2,nrm3]

    cxys <- realize $ sample nsmps truhrm

    let xys = hHead <$> cxys

    let emhrms = take emepchs $ iterate (mixtureModelExpectationMaximization xys) hrm0
        emanlls = [ average $ mixtureModelLogLikelihood hrm <$> xys | hrm <- emhrms ]


    let anllrnbl = toRenderable . execEC $ do

            goalLayout

            layout_x_axis . laxis_title .= "Epoch"
            layout_y_axis . laxis_title .= "Log-Likelihood"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [zip [(0 :: Int)..] emanlls]
                plot_lines_title .= "EM"

    sequence_ $ do

        (emhrm,k) <- zip emhrms [(0 :: Int)..]
        return . goalRenderableToSVG prjdr ("em-step-" ++ show k) 800 800
            . toRenderable $ clusterLayout emhrm cxys

    goalRenderableToSVG prjdr "cross-entropy-descent" 1200 600 anllrnbl

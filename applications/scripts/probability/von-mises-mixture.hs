{-# LANGUAGE GADTs,FlexibleContexts,TypeOperators,DataKinds #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

prjdr :: String
prjdr = "probability/von-mises-mixture"

-- Manifolds --

type Latent = Categorical Int 3
type Observable = (VonMises, VonMises)
type Harmonium' = Harmonium Tensor Observable Latent

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 1.5
muy1 = 1.5
mux2 = 3
muy2 = 3
mux3 = 4.5
muy3 = 4.5

kpx1,kpx2,kpx3,kpy1,kpy2,kpy3 :: Double
kpx1 = 2
kpy1 = 2
kpx2 = 2
kpy2 = 2
kpx3 = 2
kpy3 = 2

vm1,vm2,vm3 :: Source # Observable
vm1 = Point $ S.fromTuple (mux1, kpx1, muy1, kpy1)
vm2 = Point $ S.fromTuple (mux2, kpx2, muy2, kpy2)
vm3 = Point $ S.fromTuple (mux3, kpx3, muy3, kpy3)

vms :: S.Vector 3 (Natural # Observable)
vms = S.fromTuple (toNatural vm1,toNatural vm2,toNatural vm3)

mix1,mix2 :: Double
mix1 = 0.25
mix2 = 0.25

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Harmonium Tensor Observable Latent
truhrm = buildMixtureModel vms $ toNatural wghts

-- Mixture Distributions --

vm1',vm2',vm3' :: Source # Observable
vm1' = Point $ S.fromTuple (2, 0.5, pi, 0.5)
vm2' = Point $ S.fromTuple (3, 0.5, pi, 0.5)
vm3' = Point $ S.fromTuple (4, 0.5, pi, 0.5)

vms' :: S.Vector 3 (Natural # Observable)
vms' = S.fromTuple (toNatural vm1',toNatural vm2',toNatural vm3')

mix1',mix2' :: Double
mix1' = 0.33
mix2' = 0.33

wghts' :: Source # Latent
wghts' = Point $ S.doubleton mix1' mix2'

hrm0 :: Natural # Harmonium Tensor Observable Latent
hrm0 = buildMixtureModel vms' $ toNatural wghts'

-- Training --

nsmps :: Int
nsmps = 100

eps :: Double
eps = -0.05

bnd :: Double
bnd = 1e-5

admmlt :: Int
admmlt = 10

-- EM
nepchs :: Int
nepchs = 100

-- Functions --

vonMisesEM
    :: Sample Observable -- ^ Observations
    -> (Natural # Harmonium', S.Vector 3 (Natural # Observable))
    -> (Natural # Harmonium', S.Vector 3 (Natural # Observable))
vonMisesEM zs (hrm,nzs0) =
    let zs' = hSingleton <$> zs
        (cats,mzs) = deepMixtureModelExpectationStep zs' $ transposeHarmonium hrm
        nzs = S.zipWith cauchyFun (S.map fromOneHarmonium mzs) nzs0
     in (buildMixtureModel nzs cats, nzs)
    where diffFun mz nz = joinTangentPair nz $ crossEntropyDifferential mz nz
          cauchyFun mz nz = cauchyLimit euclideanDistance bnd
              $ vanillaGradientSequence (diffFun mz) eps defaultAdamPursuit nz

filterCat :: [SamplePoint Harmonium'] -> Int -> [SamplePoint Observable]
filterCat cxys n = hHead <$> filter ((== n) . hHead . hTail) cxys

elipse :: Source # (VonMises, VonMises) -> [[(Double,Double)]]
elipse vm =
    let [mux,kpx,muy,kpy] = listCoordinates vm
        f t = (mux + sqrt (recip kpx) * cos t, muy + sqrt (recip kpy) * sin t)
     in [f <$> range 0 (2*pi) 100]

plotVonMises
    :: AlphaColour Double
    -> Source # Observable
    -> EC (Layout Double Double) ()
plotVonMises clr vm = do

    let mux:_:muy:_ = listCoordinates vm

    plot . liftEC $ do

        plot_points_values .= [(mux,muy)]
        plot_points_style .= filledCircles 7 clr

    plot . liftEC $ do

        plot_lines_values .= elipse vm
        plot_lines_style .= solidLine 3 clr

plotSample
    :: AlphaColour Double
    -> Sample Harmonium'
    -> Int
    -> EC (Layout Double Double) ()
plotSample clr cxys cat =

    plot . liftEC $ do

        plot_points_values .= filterCat cxys cat
        plot_points_style .= filledCircles 3 clr


clusterLayout :: Natural # Harmonium' -> Sample Harmonium' -> Layout Double Double
clusterLayout hrm1 cxys = execEC $ do

    goalLayout

    radiansAbscissa

    let axs = [0,pi/2,pi,3*pi/2,2*pi]
    layout_y_axis . laxis_generate .= const (makeAxis (map show) (axs,axs,axs))
    layout_y_axis . laxis_override .=
        axisGridHide . axisLabelsOverride [(0,"0"),(pi/2,"π/2"),(pi,"π"),(3*pi/2,"3π/2"),(2*pi,"2π")]

    let lkl = fst $ splitBottomHarmonium hrm1

        vm1'' = toSource $ lkl >.>* 0
        vm2'' = toSource $ lkl >.>* 1
        vm3'' = toSource $ lkl >.>* 2

    plotVonMises (opaque red) vm1
    plotVonMises (opaque green) vm2
    plotVonMises (opaque blue) vm3
    plotSample (opaque red) cxys 0
    plotSample (opaque green) cxys 1
    plotSample (opaque blue) cxys 2

    plotVonMises (opaque black) vm1''
    plotVonMises (opaque black) vm2''
    plotVonMises (opaque black) vm3''


--- Main ---


main :: IO ()
main = do

    cxys <- realize $ sample nsmps truhrm

    let xys = hHead <$> cxys

    let emhrms = map fst . take nepchs $ iterate (vonMisesEM xys) (hrm0, vms')
        emanlls = [ average $ mixtureModelLogLikelihood hrm <$> xys | hrm <- emhrms ]

    let sgd hrm = joinTangentPair hrm $ stochasticMixtureModelDifferential xys hrm
        admhrms = takeEvery admmlt . take (admmlt*nepchs) $ vanillaGradientSequence sgd eps defaultAdamPursuit hrm0
        admanlls = [ average $ mixtureModelLogLikelihood hrm <$> xys | hrm <- admhrms ]

    let anllrnbl = toRenderable . execEC $ do

            goalLayout

            layout_x_axis . laxis_title .= "Epoch"
            layout_y_axis . laxis_title .= "Log-Likelihood"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [zip [(0 :: Int)..] emanlls]
                plot_lines_title .= "EM"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [zip [(0 :: Int)..] admanlls]
                plot_lines_title .= "Adam"

    sequence_ $ do

        (emhrm,k) <- zip emhrms [(0 :: Int)..]
        return . goalRenderableToSVG prjdr ("em-step-" ++ show k) 800 800
            . toRenderable $ clusterLayout emhrm cxys

    sequence_ $ do

        (admhrm,k) <- zip admhrms [(0 :: Int)..]
        return . goalRenderableToSVG prjdr ("adam-step-" ++ show k) 800 800
            . toRenderable $ clusterLayout admhrm cxys

    goalRenderableToSVG prjdr "cross-entropy-descent" 1200 600 anllrnbl

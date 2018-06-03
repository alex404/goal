{-# LANGUAGE ScopedTypeVariables,TypeFamilies,TypeOperators,FlexibleContexts,DataKinds,Arrows #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S

-- Qualified --

import qualified Criterion.Main as C


--- Globals ---


-- Manifolds --

type Latent = Categorical Int 3

-- Mixture Distributions --

mu1,mu2,mu3 :: Double
mu1 = -4
mu2 = 1
mu3 = 3

nrm0,nrm1,nrm2,nrm3 :: Source # Normal
nrm0 = Point $ S.doubleton 0 1000
nrm1 = Point $ S.doubleton mu1 0.5
nrm2 = Point $ S.doubleton mu2 0.5
nrm3 = Point $ S.doubleton mu3 0.5

nrms :: S.Vector 3 (Source # Normal)
nrms = S.fromTuple (nrm1,nrm2,nrm3)

mix1,mix2 :: Double
mix1 = 0.3
mix2 = 0.2

wghts :: Mean # Latent
wghts = Point $ S.doubleton mix1 mix2

hrm :: Natural # Harmonium Tensor Normal Latent
hrm = buildCategoricalHarmonium nrm0 nrms wghts

-- Training --

sx0 :: Source # Normal
sx0 = Point $ S.doubleton 0 20

nx0 :: Natural # Normal
nx0 = transition sx0

ceeps,ipeps,bt1,bt2,rg :: Double
ceeps = -0.005
ipeps = -0.01
bt1 = 0.9
bt2 = 0.999
rg = 1e-8

type NBatch = 5

nepchs :: Int
nepchs = 500

-- Plotting --

pltsmps :: B.Vector 100 Double
pltsmps = B.range (-5) 5


--- Main ---


main :: IO ()
main = do

    let rcedff nx = do
            vcxs <- sample hrm
            let vxs :: Sample NBatch Normal
                vxs = hHead <$> vcxs
            return $ stochasticCrossEntropyDifferential vxs nx

    vcxs <- realize $ sample hrm
    let vxs :: Sample NBatch Normal
        vxs = hHead <$> vcxs

    cedffcrc <- realize (accumulateRandomFunction0 rcedff)

    let ipdff nx = do
            (xs :: Sample NBatch Normal) <- sample nx
            return $ harmoniumInformationProjectionDifferential nx xs (transposeHarmonium hrm)

    ipdffcrc <- realize (accumulateRandomFunction0 ipdff)

    let trnchn :: Double
               -> Natural # Normal
               -> Circuit (Natural # Normal) (CotangentVector Natural Normal)
               -> Chain (Natural # Normal)
        trnchn eps nx0' dffcrc = accumulateCircuit0 nx0' $ proc ((),nx) -> do
            dnx <- dffcrc -< nx
            adamAscent eps bt1 bt2 rg -< sharp $ joinTangentPair nx dnx

    let cauchify :: [Natural # Normal] -> [Natural # Normal]
        cauchify = cauchySequence euclideanDistance bnd

        grds,mtms,adms :: [Cartesian # Euclidean 2]
        grds = cauchify $ gradientSequence eps (differential' f) p0
        --nwts = cauchify $ newtonSequence (differential' f) p0
        mtms = cauchify $ momentumSequence eps mu (differential' f) p0
        adms = cauchify $ adamSequence eps b1 b2 rg (differential' f) p0

    let cenx1 nx0' = cauchySequence relativeEntropy 0.01 $ streamChain (trnchn ceeps nx0' cedffcrc)
        ipnx1 nx0' = cauchySequence (flip relativeEntropy) 0.01 $ streamChain (trnchn ipeps nx0' ipdffcrc)

    putStrLn "Cross-Entropy Steps:"
    print $ length grds - 1

    putStrLn "Information Projection Steps:"
    print $ length mtms - 1

    let dnslyt = execEC $ do

            goalLayout

            layout_x_axis . laxis_title .= "Sample"
            layout_y_axis . laxis_title .= "Density"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [ B.toList . B.zip pltsmps $ mixtureDensity hrm <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= dashedLine 3 [10,5] (opaque purple)
                plot_lines_values .= [ B.toList . B.zip pltsmps $ density nx0 <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ B.toList . B.zip pltsmps $ density (cenx1 nx0) <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ B.toList . B.zip pltsmps $ density (ipnx1 nx0) <$> pltsmps ]

            plot . liftEC $ do

                plot_points_style .= filledCircles 5 (opaque black)
                plot_points_values .= zip (B.toList vxs) (repeat 0.2)

    goalRenderableToSVG "simulation" "projection" 1200 800 $ toRenderable dnslyt

    C.defaultMain
       [ C.bench "cross-entropy-minimization" $ C.nf cenx1 nx0
       , C.bench "information-projection" $ C.nf ipnx1 nx0 ]

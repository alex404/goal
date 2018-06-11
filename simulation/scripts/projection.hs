{-# LANGUAGE ScopedTypeVariables,TypeFamilies,TypeOperators,FlexibleContexts,DataKinds,Arrows #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

-- Qualified --

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S


--- Globals ---


-- Manifolds --

type Latent = Categorical Int 3

-- Mixture Distributions --

mu1,mu2,mu3 :: Double
mu1 = -4
mu2 = 0.5
mu3 = 3

nrm0,nrm1,nrm2,nrm3 :: Source # Normal
nrm0 = Point $ S.doubleton 0 1000
nrm1 = Point $ S.doubleton mu1 0.6
nrm2 = Point $ S.doubleton mu2 0.4
nrm3 = Point $ S.doubleton mu3 0.5

nrms :: S.Vector 3 (Source # Normal)
nrms = S.fromTuple (nrm1,nrm2,nrm3)

mix1,mix2 :: Double
mix1 = 0.3
mix2 = 0.3

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

type NBatch = 50

-- Plotting --

pltsmps :: B.Vector 100 Double
pltsmps = B.range (-7) 7


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

    let cenxs = cauchySequence relativeEntropy 1e-6 $ streamChain (trnchn ceeps nx0 cedffcrc)
        ipnxs = cauchySequence (flip relativeEntropy) 1e-6 $ streamChain (trnchn ipeps nx0 ipdffcrc)

    putStrLn "Cross-Entropy Steps:"
    print $ length cenxs - 1

    putStrLn "Information Projection Steps:"
    print $ length ipnxs - 1

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
                plot_lines_values .= [ B.toList . B.zip pltsmps $ density (last cenxs) <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ B.toList . B.zip pltsmps $ density (last ipnxs) <$> pltsmps ]

            plot . liftEC $ do

                plot_points_style .= filledCircles 5 (opaque black)
                plot_points_values .= zip (B.toList vxs) (repeat 0.2)

    goalRenderableToSVG "simulation" "projection" 1200 800 $ toRenderable dnslyt

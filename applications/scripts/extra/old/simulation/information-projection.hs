{-# LANGUAGE ScopedTypeVariables,TypeFamilies,TypeOperators,FlexibleContexts,DataKinds,Arrows #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

-- Qualified --

import qualified Goal.Core.Vector.Storable as S


--- Globals ---


-- Manifolds --

type Latent = Categorical Int 3

-- Mixture Distributions --

mu1,mu2,mu3 :: Double
mu1 = -4
mu2 = 0.5
mu3 = 3

nrm1,nrm2,nrm3 :: Source # Normal
nrm1 = Point $ S.doubleton mu1 0.6
nrm2 = Point $ S.doubleton mu2 0.4
nrm3 = Point $ S.doubleton mu3 0.5

nrms :: S.Vector 3 (Natural # Normal)
nrms = S.fromTuple (toNatural nrm1,toNatural nrm2,toNatural nrm3)

mix1,mix2 :: Double
mix1 = 0.3
mix2 = 0.3

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

hrm :: Natural # Harmonium Tensor Normal Latent
hrm = buildMixtureModel nrms $ toNatural wghts

-- Training --

sx0 :: Source # Normal
sx0 = Point $ S.doubleton 0 20

nx0 :: Natural # Normal
nx0 = transition sx0

ceeps,ipeps :: Double
ceeps = -0.005
ipeps = -0.01

nbtch :: Int
nbtch = 50

-- Plotting --

pltsmps :: [Double]
pltsmps = range (-7) 7 100


--- Main ---


main :: IO ()
main = do

    let rcedff nx = do
            vcxs <- sample nbtch hrm
            let vxs = hHead <$> vcxs
            return $ stochasticCrossEntropyDifferential vxs nx

    vcxs <- realize $ sample nbtch hrm
    let vxs = hHead <$> vcxs

    cedffcrc <- realize (accumulateRandomFunction0 rcedff)

    let ipdff nx = do
            xs <- sample nbtch nx
            return $ harmoniumInformationProjectionDifferential nx xs (transposeHarmonium hrm)

    ipdffcrc <- realize (accumulateRandomFunction0 ipdff)

    let trnchn :: Double
               -> Natural # Normal
               -> Circuit (Natural # Normal) (CotangentVector Natural Normal)
               -> Chain (Natural # Normal)
        trnchn eps nx0' dffcrc = accumulateCircuit0 nx0' $ proc ((),nx) -> do
            dnx <- dffcrc -< nx
            gradientCircuit eps defaultAdamPursuit -< sharp $ joinTangentPair nx dnx

    let relativeEntropy' p = relativeEntropy (dualTransition p)
        relativeEntropy'' q p = relativeEntropy (dualTransition p) q

    let cenxs = cauchySequence relativeEntropy' 1e-6 $ streamChain (trnchn ceeps nx0 cedffcrc)
        ipnxs = cauchySequence relativeEntropy'' 1e-6 $ streamChain (trnchn ipeps nx0 ipdffcrc)

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
                plot_lines_values .= [ zip pltsmps $ mixtureDensity hrm <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= dashedLine 3 [10,5] (opaque purple)
                plot_lines_values .= [ zip pltsmps $ density nx0 <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip pltsmps $ density (last cenxs) <$> pltsmps ]

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip pltsmps $ density (last ipnxs) <$> pltsmps ]

            plot . liftEC $ do

                plot_points_style .= filledCircles 5 (opaque black)
                plot_points_values .= zip vxs (repeat 0.2)

    goalRenderableToSVG "simulation" "projection" 1200 800 $ toRenderable dnslyt

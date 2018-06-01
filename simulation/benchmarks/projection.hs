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
type Observable = MeanNormal (1/2)

-- Mixture Distributions --

mu1,mu2,mu3 :: Double
mu1 = -4
mu2 = 0
mu3 = 3

nrm0,nrm1,nrm2,nrm3 :: Source # Observable
nrm0 = zero
nrm1 = Point $ S.singleton mu1
nrm2 = Point $ S.singleton mu2
nrm3 = Point $ S.singleton mu3

nrms :: S.Vector 3 (Source # Observable)
nrms = S.fromTuple (nrm1,nrm2,nrm3)

mix1,mix2 :: Double
mix1 = 0.2
mix2 = 0.5

wghts :: Mean # Latent
wghts = Point $ S.doubleton mix1 mix2

hrm :: Natural # Harmonium Tensor Observable Latent
hrm = buildCategoricalHarmonium nrm0 nrms wghts

-- Training --

sx0 :: Source # Normal
sx0 = Point $ S.doubleton (-3) 8

nx0 :: Natural # Normal
nx0 = transition sx0

eps,bt1,bt2,rg :: Double
eps = -0.01
bt1 = 0.9
bt2 = 0.999
rg = 1e-8

type NBatch = 10

nepchs :: Int
nepchs = 100

-- Plotting --

pltsmps :: B.Vector 100 Double
pltsmps = B.range (-5) 5


--- Main ---


main :: IO ()
main = do

    let rscdff nx = do
            vcxys <- sample hrm
            let vxys :: Sample NBatch Observable
                vxys = hHead <$> vcxys
            return $ stochasticCrossEntropyDifferential vxys nx

    scdffcrc <- realize (accumulateRandomFunction0 rscdff)

    let ipdff nx = do
            (xs :: Sample NBatch Normal) <- sample nx
            return $ harmoniumInformationProjectionDifferential nx xs hrm

    ipdffcrc <- realize (accumulateRandomFunction0 ipdff)

    let trnchn :: Natural # Normal
               -> Circuit (Natural # Normal) (CotangentVector Natural Normal)
               -> Chain (Natural # Normal)
        trnchn nx0' dffcrc = accumulateCircuit0 nx0' $ proc ((),nx) -> do
            dnx <- dffcrc -< nx
            adamAscent eps bt1 bt2 rg -< sharp $ joinTangentPair nx dnx

    let cenx1 nx0' = streamChain (trnchn nx0' scdffcrc) !! nepchs
        ipnx1 nx0' = streamChain (trnchn nx0' ipdffcrc) !! nepchs

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

    goalRenderableToSVG "simulation" "projection" 1200 800 $ toRenderable dnslyt

    C.defaultMain
       [ C.bench "cross-entropy-minimization" $ C.nf cenx1 nx0
       , C.bench "information-projection" $ C.nf ipnx1 nx0 ]

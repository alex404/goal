{-# LANGUAGE FlexibleContexts,Arrows #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Program ---


-- Globals --

-- True Population Code
truvr0 = 0.5
trumu0 = 1
trusp0 = Standard # fromList VonMises [trumu0,recip truvr0]

truvr = 2
trumn = -pi
trumx = pi
trunkrns = 10
trumus = range trumn trumx trunkrns
trusps = tail [ fromList VonMises [mu,truvr] | mu <- trumus]
trugn = 2
truppc = vonMisesPopulationEncoder trusps trugn

truhrm = affineToHarmonium (transition trusp0) truppc

{-
-- Initial Population Code
trnvr0 = 10
trnmu0 = 0
trnsp0 = Standard # fromList VonMises [trnmu0,recip trnvr0]

trnvr = 10
trnmn = -1
trnmx = 1
trnnkrns = 10
trnmus = range trnmn trnmx trnnkrns
trnsps = tail [ fromList VonMises [mu,trnvr] | mu <- trnmus]
trngn = 1
trnppc0 = vonMisesPopulationEncoder trnsps trngn

--trnhrm0 = affineToHarmonium (transition trnsp0) trnppc0
-}

-- Training
--nsmps = 10^7
blkn = 20
nstpstrn = 2*10^4
cdn = 1
eps = 0.00005
bt1 = 0.9
bt2 = 0.999
rg = 1e-8

-- Plot
(pltmn,pltmx) = (trumn, trumx)
pltxs = range pltmn pltmx 200
avg = 100

-- Main --

main = do

        {-
    x0s <- runWithSystemRandom . replicateM nsmps $ generate trusp0
    n0s <- runWithSystemRandom . mapM generate $ truppc >$>* x0s
    let responseGenerator = generator $ breakEvery blkn n0s
        -}
    responseGenerator <- runWithSystemRandom (accumulateRandomFunction0
        (\() -> replicateM blkn (generate trusp0) >>= (\xs -> mapM standardGenerate (truppc >$>* xs))))
    trnhrm0 <- runWithSystemRandom . initialize (Standard # fromList Normal [0,1]) $ manifold truhrm
    let ls0 = zero VonMises
        --trnhrm0' = let (_,os,tns) = splitHarmonium trnhrm0 in joinHarmonium ls0 os tns
        trnhrm0' = trnhrm0
    contrastor <- runWithSystemRandom $ accumulateRandomFunction0 (uncurry (bulkContrastiveDivergence cdn))
    --rectifier <- runWithSystemRandom $ accumulateRandomFunction0 (uncurry rectificationDifferentials)

    let trnchn = accumulateMealy (trnhrm0',0) $ proc ((),(hrm,k)) -> do
            ns <- responseGenerator -< ()
            dhrm <- contrastor -< (ns,hrm)
            --dhrm' <- rectifier -< (ns,hrm)
            hrm' <- adamDescent eps bt1 bt2 rg -< dhrm --averagePoint [dhrm,dhrm']
            --hrm' <- vanillaGradientDescent eps -< averagePoint [dhrm,dhrm']
            {-
            let (lbs,_,_) = splitHarmonium hrm
                (_,obs',tns') = splitHarmonium hrm'
                hrm'' = joinHarmonium lbs obs' tns'
                -}
            returnA -< (hrm',(hrm',k+1))

        trnhrms = streamChain trnchn
        trnhrm1 = trnhrms !! nstpstrn
        (trno1,trnl1,trnaff1) = splitHarmonium $ harmoniumTranspose trnhrm1
        trnppc1 = joinAffine trno1 trnaff1

    let tcrnbl = toRenderable . execEC $ do

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= (map (zip pltxs) . transpose $ listCoordinates . transitionTo Standard <$> truppc >$>* pltxs)
                plot_lines_title .= "True Tuning Curves"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= (map (zip pltxs) . transpose $ listCoordinates . transitionTo Standard <$> trnppc1 >$>* pltxs)
                plot_lines_title .= "Learned Tuning Curves"

    let mnrnbl = toRenderable . execEC $ do

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque red)
                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average $ breakEvery avg [ let (trnl,_,_) = splitHarmonium trnhrm in coordinate 0 $ transitionTo Standard trnl | trnhrm <- take nstpstrn trnhrms ]  ]
                plot_lines_title .= "Latent Mean"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average $ breakEvery avg [ let (trnl,_,_) = splitHarmonium trnhrm in coordinate 1 $ transitionTo Standard trnl | trnhrm <- take nstpstrn trnhrms ]  ]
                plot_lines_title .= "Latent Variance"

    let msernbl = toRenderable . execEC $

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average
                    $ breakEvery avg
                    [ let (_,obs,iprms) = splitHarmonium trnhrm
                          ppc' = joinAffine obs (matrixTranspose iprms)
                       in average [ divergence (ppc' >.>* x) (transition $ truppc >.>* x) | x <- range (-pi) pi 100 ]
                    | trnhrm <- take nstpstrn trnhrms ] ]
                plot_lines_title .= "Latent Mean"


    goalRenderableToSVG "simulation" "contrastive-divergence" 800 800 . gridToRenderable . weights (1,1) $ tval mnrnbl ./. tval tcrnbl ./. tval msernbl

    print trusp0
    print . transitionTo Standard $ trnl1
    goalWriteFile "simulation" "contrastive-divergence" . show $ listCoordinates trnhrm1

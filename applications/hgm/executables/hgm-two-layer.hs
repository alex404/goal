{-# LANGUAGE DataKinds,FlexibleContexts,Arrows,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

--- Program ---


-- Globals --

type NNeurons = 10

mn,mx :: Double
mn = 0
mx = 2*pi

-- True PPC Harmonium
truvm0 :: Source # VonMises
truvm0 = Point $ S.doubleton 1 2

rndvm0 :: Random s (Source # VonMises)
rndvm0 = uniformInitialize $ B.doubleton (0,2*pi) (0,10)

truprs :: Double
truprs = 2

trutcs :: S.Vector NNeurons (Source # VonMises)
trutcs = S.init . S.map (\mu -> Point $ S.doubleton mu truprs) $ S.range mn mx

rndtcs :: Random s (S.Vector NNeurons (Source # VonMises))
rndtcs = S.replicateM (uniformInitialize $ B.doubleton (0,2*pi) (0,5))

trugn :: Double
trugn = 2

truppc :: Mean ~> Natural # Affine Tensor (Replicated NNeurons Poisson) VonMises
truppc = vonMisesPopulationEncoder trutcs trugn

truhrm :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
truhrm = let (nrnb,mtx) = splitAffine truppc
          in joinBottomHarmonium nrnb (transpose mtx) . toOneHarmonium $ transition truvm0

rndhrm :: Random s (Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises)
rndhrm = do
    vm0 <- rndvm0
    tcs <- rndtcs
    let (nrnb,mtx) = splitAffine $ vonMisesPopulationEncoder tcs 1
    return $ joinBottomHarmonium nrnb (transpose mtx) . toOneHarmonium $ transition vm0

-- Training --

cdn :: Int
cdn = 1

tprxy :: Proxy 100
tprxy = Proxy

vprxy :: Proxy 100
vprxy = Proxy

eps,bt1,bt2,rg :: Double
eps = -0.01
bt1 = 0.9
bt2 = 0.999
rg = 1e-8

type NBatch = 10

nepchs,trnepchn :: Int
nepchs = 100
trnepchn = 1000


--- Plot ---

type NPlot = 100

pltsmps :: Sample NPlot VonMises
pltsmps = B.range mn mx

harmoniumTuningCurves
    :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
    -> Layout Double Double
harmoniumTuningCurves hrm =
    let (nb,mtx,_) = splitBottomHarmonium hrm
        tcs = tuningCurves pltsmps $ joinAffine nb (transpose mtx)
     in execEC $ do

            goalLayout
            radiansAbscissa

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= B.toList (B.toList <$> tcs)



--- Main ---


main :: IO ()
main = do

    hrm0 <- realize rndhrm

    dffcrc <- realize (accumulateRandomFunction0 $ uncurry (contrastiveDivergence cdn))

    smpchn <- realize $ accumulateRandomFunction0 (const $ fmap hHead <$> sampleRectified zero truhrm)

    let trncrc
            :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
            -> Circuit (Sample NBatch (Replicated NNeurons Poisson))
                   (Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises)
        trncrc hrm0' = accumulateCircuit0 hrm0' $ proc (xs,hrm) -> do
            dhrm <- dffcrc -< (xs,hrm)
            let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
            adamAscent eps bt1 bt2 rg -< dhrmpr

    let hrm1 = last . take nepchs . takeEvery trnepchn . streamChain $ trncrc hrm0 <<< smpchn

        trurnbl = toRenderable $ harmoniumTuningCurves truhrm
        rnbl0 = toRenderable $ harmoniumTuningCurves hrm0
        rnbl1 = toRenderable $ harmoniumTuningCurves hrm1

    goalRenderableToSVG "hgm/two-layer" "true-tuning-curves" 1200 800 trurnbl
    goalRenderableToSVG "hgm/two-layer" "initial-tuning-curves" 1200 800 rnbl0
    goalRenderableToSVG "hgm/two-layer" "final-tuning-curves" 1200 800 rnbl1


---- Plot
--(pltmn,pltmx) = (trumn, trumx)
--pltxs = range pltmn pltmx 200
--avg = 100
--
---- Main --
--
--main = do
--
--        {-
--    x0s <- runWithSystemRandom . replicateM nsmps $ generate trusp0
--    n0s <- runWithSystemRandom . mapM generate $ truppc >$>* x0s
--    let responseGenerator = generator $ breakEvery blkn n0s
--        -}
--    responseGenerator <- runWithSystemRandom (accumulateRandomFunction0
--        (\() -> replicateM blkn (generate trusp0) >>= (\xs -> mapM standardGenerate (truppc >$>* xs))))
--    trnhrm0 <- runWithSystemRandom . initialize (Standard # fromList Normal [0,1]) $ manifold truhrm
--    let ls0 = zero VonMises
--        --trnhrm0' = let (_,os,tns) = splitHarmonium trnhrm0 in joinHarmonium ls0 os tns
--        trnhrm0' = trnhrm0
--    contrastor <- runWithSystemRandom $ accumulateRandomFunction0 (uncurry (bulkContrastiveDivergence cdn))
--    --rectifier <- runWithSystemRandom $ accumulateRandomFunction0 (uncurry rectificationDifferentials)
--
--    let trnchn = accumulateMealy (trnhrm0',0) $ proc ((),(hrm,k)) -> do
--            ns <- responseGenerator -< ()
--            dhrm <- contrastor -< (ns,hrm)
--            --dhrm' <- rectifier -< (ns,hrm)
--            hrm' <- adamDescent eps bt1 bt2 rg -< dhrm --averagePoint [dhrm,dhrm']
--            --hrm' <- vanillaGradientDescent eps -< averagePoint [dhrm,dhrm']
--            {-
--            let (lbs,_,_) = splitHarmonium hrm
--                (_,obs',tns') = splitHarmonium hrm'
--                hrm'' = joinHarmonium lbs obs' tns'
--                -}
--            returnA -< (hrm',(hrm',k+1))
--
--        trnhrms = streamChain trnchn
--        trnhrm1 = trnhrms !! nstpstrn
--        (trno1,trnl1,trnaff1) = splitHarmonium $ harmoniumTranspose trnhrm1
--        trnppc1 = joinAffine trno1 trnaff1
--
--    let tcrnbl = toRenderable . execEC $ do
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque blue)
--                plot_lines_values .= (map (zip pltxs) . transpose $ listCoordinates . transitionTo Standard <$> truppc >$>* pltxs)
--                plot_lines_title .= "True Tuning Curves"
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque red)
--                plot_lines_values .= (map (zip pltxs) . transpose $ listCoordinates . transitionTo Standard <$> trnppc1 >$>* pltxs)
--                plot_lines_title .= "Learned Tuning Curves"
--
--    let mnrnbl = toRenderable . execEC $ do
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque red)
--                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average $ breakEvery avg [ let (trnl,_,_) = splitHarmonium trnhrm in coordinate 0 $ transitionTo Standard trnl | trnhrm <- take nstpstrn trnhrms ]  ]
--                plot_lines_title .= "Latent Mean"
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque blue)
--                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average $ breakEvery avg [ let (trnl,_,_) = splitHarmonium trnhrm in coordinate 1 $ transitionTo Standard trnl | trnhrm <- take nstpstrn trnhrms ]  ]
--                plot_lines_title .= "Latent Variance"
--
--    let msernbl = toRenderable . execEC $
--
--            plot . liftEC $ do
--
--                plot_lines_style .= solidLine 3 (opaque black)
--                plot_lines_values .= [ zip [(0 :: Int),avg..] . map average
--                    $ breakEvery avg
--                    [ let (_,obs,iprms) = splitHarmonium trnhrm
--                          ppc' = joinAffine obs (matrixTranspose iprms)
--                       in average [ divergence (ppc' >.>* x) (transition $ truppc >.>* x) | x <- range (-pi) pi 100 ]
--                    | trnhrm <- take nstpstrn trnhrms ] ]
--                plot_lines_title .= "Latent Mean"
--
--
--    goalRenderableToSVG "simulation" "contrastive-divergence" 800 800 . gridToRenderable . weights (1,1) $ tval mnrnbl ./. tval tcrnbl ./. tval msernbl
--
--    print trusp0
--    print . transitionTo Standard $ trnl1
--    goalWriteFile "simulation" "contrastive-divergence" . show $ listCoordinates trnhrm1

{-# LANGUAGE TypeFamilies,ScopedTypeVariables,DataKinds,FlexibleContexts,Arrows,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B


--- Program ---


-- Variables --

-- General
type NNeurons = 10

mn,mx :: Double
mn = 0
mx = 2*pi

prc :: Double
prc = 2

tcs :: S.Vector NNeurons (Source # VonMises)
tcs = S.init . S.map (\mu -> Point $ S.doubleton mu prc) $ S.range mn mx

type NSamples = 100

xsmp :: Sample NSamples VonMises
xsmp = B.init $ B.range mn mx

-- Initial
vm0 :: Natural # VonMises
vm0 = zero

gn0 :: Double
gn0 = 1

ppc0 :: Mean ~> Natural # Affine Tensor (Replicated NNeurons Poisson) VonMises
ppc0 = vonMisesPopulationEncoder tcs gn0

hrm0 :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
hrm0 = joinBottomHarmonium ppc0 $ toOneHarmonium vm0

rho00 :: Double
rprms0 :: Natural # VonMises
(rho00,rprms0) = populationCodeRectificationParameters ppc0 xsmp

-- True
truvm0 :: Source # VonMises
truvm0 = Point $ S.doubleton 1 2

trurho0 :: Double
trurho0 = 2 * realToFrac (natVal (Proxy :: Proxy NNeurons))

trurprms :: Natural # VonMises
trurprms = Point $ S.doubleton 2 4

truppc :: Mean ~> Natural # R NNeurons Poisson <* VonMises
truppc = rectifyPopulationCode trurho0 trurprms xsmp $ vonMisesPopulationEncoder tcs 1

truhrm :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
truhrm = joinBottomHarmonium truppc . toOneHarmonium $ transition truvm0

-- Training --

cdn :: Int
cdn = 1

eps,bt1,bt2,rg :: Double
eps = -0.01
bt1 = 0.9
bt2 = 0.999
rg = 1e-8

type TBatch = 10
type IPBatch = 10
type RBatch = 10

nepchs,trnepchn :: Int
nepchs = 10
trnepchn = 10


--- Plot ---

type NPlot = 100

pltsmps :: Sample NPlot VonMises
pltsmps = B.range mn mx

harmoniumTuningCurves
    :: Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
    -> Layout Double Double
harmoniumTuningCurves hrm =
    let tcs' = tuningCurves pltsmps . fst $ splitBottomHarmonium hrm
     in execEC $ do

            goalLayout
            radiansAbscissa

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque black)
                plot_lines_values .= B.toList (B.toList <$> tcs')


--- Main ---


main :: IO ()
main = do

    -- Generic
    smpchn <- realize $ accumulateRandomFunction0 (const $ fmap hHead <$> sampleRectified zero truhrm)

    -- Contrastive Divergence
    cdcrc <- realize (accumulateRandomFunction0 $ uncurry (contrastiveDivergence cdn))

    let cdtrncrc
            :: Circuit (Sample TBatch (Replicated NNeurons Poisson))
                   (Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises)
        cdtrncrc = accumulateCircuit0 hrm0 $ proc (xs,hrm) -> do
            dhrm <- cdcrc -< (xs,hrm)
            let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
            adamAscent eps bt1 bt2 rg -< dhrmpr

    -- Rectified

    rccrc <- realize (accumulateRandomFunction0 (uncurry (uncurry estimateRectifiedHarmoniumDifferentials)))
    (smpcrc :: Circuit (Natural # VonMises) (Sample RBatch VonMises)) <- realize (accumulateRandomFunction0 sample)
--    (ipcrc :: Circuit (Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises, Natural # VonMises)
--        (Natural # VonMises)) <- realize (accumulateRandomFunction0 (uncurry $ harmoniumInformationProjection (Proxy :: Proxy 10) eps bt1 bt2 rg 10))

    let rctrncrc
            :: Circuit (Sample TBatch (Replicated NNeurons Poisson))
                   (Natural # VonMises, Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises,Int)
        rctrncrc = accumulateCircuit0 (rprms0,hrm0,0) $ proc (xs,(rprms,hrm,k)) -> do
            dhrm <- rccrc -< ((xs,rprms),hrm)
            let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
            hrm' <- adamAscent eps bt1 bt2 rg -< dhrmpr
            let (lkl',xp') = splitBottomHarmonium hrm'
            --rsmp <- smpcrc -< rprms
            let (rho0,rprms') = populationCodeRectificationParameters lkl' xsmp
                lkl'' = rectifyPopulationCode rho0 rprms' xsmp lkl'
            returnA -< (traceGiven rprms',joinBottomHarmonium lkl'' xp',k+1)

    -- Simulation

    let cdhrm1 = last . take nepchs . takeEvery trnepchn . streamChain $ cdtrncrc <<< smpchn
        (rprms1,rchrm1,n) = last . take nepchs . takeEvery trnepchn . streamChain $ rctrncrc <<< smpchn

        trurnbl = toRenderable $ harmoniumTuningCurves truhrm
        rnbl0 = toRenderable $ harmoniumTuningCurves hrm0
        rnbl1 = toRenderable $ harmoniumTuningCurves cdhrm1
        rnbl2 = toRenderable $ harmoniumTuningCurves rchrm1

    print n

    goalRenderableToSVG "hgm/two-layer" "true-tuning-curves" 1200 800 trurnbl
    goalRenderableToSVG "hgm/two-layer" "initial-tuning-curves" 1200 800 rnbl0
    goalRenderableToSVG "hgm/two-layer" "cd-tuning-curves" 1200 800 rnbl1
    goalRenderableToSVG "hgm/two-layer" "rectification-tuning-curves" 1200 800 rnbl2

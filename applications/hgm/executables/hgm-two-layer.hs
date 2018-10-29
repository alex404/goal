{-# LANGUAGE TypeFamilies,ScopedTypeVariables,DataKinds,FlexibleContexts,Arrows,TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Goal.Core.Vector.Storable as S


--- Program ---


-- Variables --

-- General
type NNeurons = 10

mn,mx :: Double
mn = 0
mx = 2*pi

prc :: Double
prc = 2

mus :: S.Vector NNeurons Double
mus = S.init $ S.range mn mx

tcs :: S.Vector NNeurons (Source # VonMises)
tcs = S.map (\mu -> Point $ S.doubleton mu prc) mus

nsmps :: Int
nsmps = 100

type Harmonium' = Harmonium Tensor (Replicated NNeurons Poisson) VonMises

xsmp :: Sample VonMises
xsmp = init $ range mn mx nsmps

-- Initial
vm0 :: Natural # VonMises
vm0 = zero

gn0 :: Double
gn0 = 1

ppc0 :: Mean #> Natural # Affine Tensor (Replicated NNeurons Poisson) VonMises
ppc0 = vonMisesPopulationEncoder False (Left gn0) tcs

hrm0 :: Natural # Harmonium'
hrm0 = joinBottomHarmonium ppc0 $ toOneHarmonium vm0

rho00 :: Double
rprms0 :: Natural # VonMises
(rho00,rprms0) = populationCodeRectificationParameters ppc0 xsmp

-- True
truvm0 :: Source # VonMises
truvm0 = Point $ S.doubleton 2 1

trurho0 :: Double
trurho0 = 2 * realToFrac (natVal (Proxy :: Proxy NNeurons))

trurprms1 :: Natural # VonMises
trurprms1 = Point $ S.doubleton 2 4

truppc1 :: Mean #> Natural # R NNeurons Poisson <* VonMises
truppc1 = rectifyPopulationCode trurho0 trurprms1 xsmp $ vonMisesPopulationEncoder False (Left 1) tcs

truhrm1 :: Natural # Harmonium'
truhrm1 = joinBottomHarmonium truppc1 . toOneHarmonium $ transition truvm0

truppc2 :: Mean #> Natural # Affine Tensor (Replicated NNeurons Poisson) VonMises
truppc2 = vonMisesPopulationEncoder False (Right $ S.fromTuple (1,2,3,2,1,2,3,2,1,2)) tcs

truhrm2 :: Natural # Harmonium'
truhrm2 = joinBottomHarmonium truppc2 . toOneHarmonium $ transition truvm0

trurprms2 :: Natural # VonMises
trurprms2 = snd $ populationCodeRectificationParameters truppc2 xsmp

truhrm :: Natural # Harmonium'
truhrm = truhrm2

trurprms :: Natural # VonMises
trurprms = trurprms2

-- Training --

cdn :: Int
cdn = 10

eps,ipeps :: Double
eps = -0.005
ipeps = -0.05

tbtch,ipbtch,rbtch :: Int
tbtch = 10
ipbtch = 10
rbtch = 10

nepchs,trnepchn :: Int
nepchs = 100
trnepchn = 100

ngstps :: Int
ngstps = 10

-- Circuits --

rectifiedSampleChain :: s ~> Chain (Sample (Replicated NNeurons Poisson))
rectifiedSampleChain =
    accumulateRandomChain (fmap hHead <$> sampleRectifiedHarmonium tbtch (toSingletonSum trurprms) truhrm)

gibbsSample
    :: Natural # VonMises
    -> Natural # Harmonium'
    -> Int
    -> s ~> Sample (Replicated NNeurons Poisson)
gibbsSample rprms hrm n = do
    zxs0 <- sampleRectifiedHarmonium tbtch (toSingletonSum rprms) hrm
    gchn <- bulkGibbsChain hrm zxs0
    return . fmap hHead $ streamChain gchn !! n


cdTrainingCircuit
    :: s ~> Sample (Replicated NNeurons Poisson)
        >>> Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
cdTrainingCircuit = do
    cdcrc <- accumulateRandomFunction0 $ uncurry (contrastiveDivergence cdn)
    return . accumulateCircuit0 hrm0 $ proc (xs,hrm) -> do
        dhrm <- cdcrc -< (xs,hrm)
        let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
        gradientCircuit eps defaultAdamPursuit -< dhrmpr

rcTrainingCircuit
    :: s ~> Sample (Replicated NNeurons Poisson)
        >>> (Natural # VonMises, Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises)
rcTrainingCircuit = do
    rccrc <- accumulateRandomFunction0 (\(x,y,z) -> stochasticRectifiedHarmoniumDifferential x y z)
    return . accumulateCircuit0 (rprms0,hrm0) $ proc (xs,(rprms,hrm)) -> do
        dhrm <- rccrc -< (xs,rprms,hrm)
        let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
        hrm' <- gradientCircuit eps defaultAdamPursuit -< dhrmpr
        let (lkl',xp') = splitBottomHarmonium hrm'
            (rho0,rprms') = populationCodeRectificationParameters lkl' xsmp
            lkl'' = rectifyPopulationCode rho0 rprms' xsmp lkl'
        returnA -< (rprms',joinBottomHarmonium lkl'' xp')

ipTrainingCircuit
    :: s ~> Sample (Replicated NNeurons Poisson)
           >>> (Natural # VonMises, Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises)
ipTrainingCircuit = do
    rccrc <- accumulateRandomFunction0 (\(zs,rprms,hrm) -> stochasticRectifiedHarmoniumDifferential zs rprms hrm)
    ipcrc <- accumulateRandomFunction0 (\(rprms,hrm) ->
        harmoniumInformationProjection ipbtch ipeps defaultAdamPursuit 20 hrm rprms)
    --dcdcrc <- accumulateRandomFunction0 (uncurry (dualContrastiveDivergence cdn (Proxy :: Proxy IPBatch)))
    return . accumulateCircuit0 (rprms0,hrm0) $ proc (zs,(rprms,hrm)) -> do
        dhrm <- rccrc -< (zs,rprms,hrm)
        let dhrmpr = joinTangentPair hrm (breakPoint dhrm)
        hrm' <- gradientCircuit eps defaultAdamPursuit -< dhrmpr
        rprms' <- ipcrc -< (rprms,hrm')
        --dhrm' <- dcdcrc -< (rprms',hrm')
        --let dhrmpr' = joinTangentPair hrm' (breakPoint dhrm')
        --hrm'' <- gradientCircuit eps defaultAdamPursuit -< dhrmpr'
        returnA -< (rprms',hrm')

--- Plot ---

nplt :: Int
nplt = 100

pltsmp :: Sample VonMises
pltsmp = range mn mx nplt

numericalPrior :: Natural # Harmonium' -> SamplePoint VonMises -> Double
numericalPrior hrm x =
    let uprior = unnormalizedHarmoniumObservableDensity (transposeHarmonium hrm)
        nrm = integrate 1e-500 uprior (-pi) pi
     in uprior x / nrm

populationSampleLayout
    :: [String]
    -> [AlphaColour Double]
    -> [Sample (R NNeurons Poisson)]
    -> Layout Double Double
populationSampleLayout ttls clrs zsmps = execEC $ do

    goalLayout

    sequence_ $ do

        (ttl,clr,zsmp) <- zip3 ttls clrs zsmps

        let mz :: Mean # R NNeurons Poisson
            mz = sufficientStatisticT zsmp

        return $ do

            plot . liftEC $ do

                plot_points_style .= filledCircles 5 clr
                plot_points_values .= zip (S.toList mus) (listCoordinates mz)

            plot . liftEC $ do

                plot_lines_title .= ttl
                plot_lines_style .= solidLine 3 clr
                plot_lines_values .= [zip (S.toList mus) (listCoordinates mz)]

harmoniumTuningCurves
    :: Maybe (Natural # VonMises)
    -> Natural # Harmonium Tensor (Replicated NNeurons Poisson) VonMises
    -> LayoutLR Double Double Double
harmoniumTuningCurves mrprms hrm = execEC $ do

    let tcs' = tuningCurves pltsmp . fst $ splitBottomHarmonium hrm
        prr0 = numericalPrior hrm

    goalLayoutLR
    radiansAbscissaLR
    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,1.5)

    plotLeft . liftEC $ do


        plot_lines_title .= "Tuning Curves"
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= tcs'

    plotRight . liftEC $ do

        plot_lines_title .= "Numerical Prior"
        plot_lines_style .= solidLine 6 (opaque blue)
        plot_lines_values .= [toList . zip pltsmp $ prr0 <$> pltsmp]

    case mrprms of

      Just rprms -> plotRight . liftEC $ do

          let prr1 = density (fromOneHarmonium . snd $ marginalizeRectifiedHarmonium (toSingletonSum rprms) hrm)
          plot_lines_title .= "Rectification Marginal"
          plot_lines_style .= solidLine 3 (opaque skyblue)
          plot_lines_values .= [zip pltsmp $ prr1 <$> pltsmp]

      Nothing -> return ()


--- Main ---


main :: IO ()
main = do


    -- Generic
    smpchn <- realize $ accumulateRandomChain (gibbsSample trurprms truhrm ngstps)

    cdtrncrc <- realize cdTrainingCircuit
    --rctrncrc <- realize rcTrainingCircuit
    iptrncrc <- realize ipTrainingCircuit

    -- Simulation

    let cdhrm1 = last . take nepchs . takeEvery trnepchn . streamChain $ cdtrncrc <<< smpchn
        (iprprms1,iphrm1) = last . take nepchs . takeEvery trnepchn . streamChain $ iptrncrc <<< smpchn

        trurnbl = toRenderable $ harmoniumTuningCurves (Just trurprms) truhrm
        rnbl0 = toRenderable $ harmoniumTuningCurves (Just zero) hrm0
        cdrnbl = toRenderable $ harmoniumTuningCurves Nothing cdhrm1
        iprnbl = toRenderable $ harmoniumTuningCurves (Just iprprms1) iphrm1

    truzsmp <- realize $ gibbsSample trurprms truhrm ngstps
    let cdrprms1 = snd $ populationCodeRectificationParameters (fst $ splitBottomHarmonium cdhrm1) xsmp
    cdzsmp <- realize $ gibbsSample cdrprms1 cdhrm1 ngstps
    ipzsmp <- realize $ gibbsSample iprprms1 iphrm1 ngstps

    let ttls = ["True Gibbs", "Contrastive Divergence", "Information Projection"]
        clrs = [opaque red, opaque green, opaque blue,opaque purple]
        smprnbl = toRenderable $ populationSampleLayout ttls clrs [truzsmp,cdzsmp,ipzsmp]

    zxs0 <- realize $ sampleRectifiedHarmonium nplt (toSingletonSum trurprms) truhrm
    gchn <- realize $ bulkGibbsChain truhrm zxs0

    let stps = [0,10,100,1000]
        prgsmps = ((fmap hHead <$> streamChain gchn) !!) <$> stps
        prgttls = (++ " Steps") . show <$> stps
        prgrnbl = toRenderable $ populationSampleLayout prgttls clrs prgsmps


    goalRenderableToSVG "hgm/two-layer" "true-tuning-curves" 800 600 trurnbl
    goalRenderableToSVG "hgm/two-layer" "initial-tuning-curves" 800 600 rnbl0
    goalRenderableToSVG "hgm/two-layer" "cd-tuning-curves" 800 600 cdrnbl
    goalRenderableToSVG "hgm/two-layer" "ip-tuning-curves" 800 600 iprnbl
    goalRenderableToSVG "hgm/two-layer" "sample-comparison" 800 600 smprnbl
    goalRenderableToSVG "hgm/two-layer" "progressive-sampling" 800 600 prgrnbl


{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import KohnLab

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import Data.List
import qualified Data.Map as M

--- Globals ---


-- Training --

nepchs :: Int
nepchs = 10000

eps :: Double
eps = -0.005

-- Adam
b1,b2,rg :: Double
b1 = 0.9
b2 = 0.999
rg = 1e-8


--- Functions ---


stimulusSpikeRate
    :: Maybe NeuronID
    -> M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])
    -> [(Stimulus,Double)]
stimulusSpikeRate (Just nrn) stmttls = do
    (stm,(k,nmp)) <- M.toList stmttls
    return (stm, (/ fromIntegral k) . genericLength $ nmp M.! nrn)
stimulusSpikeRate Nothing stmttls = do
    (stm,(k,nmp)) <- M.toList stmttls
    return (stm, (/ fromIntegral k) . sum . M.elems $ genericLength <$> nmp)


--- Plots ---


-- Raw Data --

stimulusTuningCurves
    :: Double
    -> M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])
    -> M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])
    -> Maybe NeuronID
    -> LayoutLR Double Double Double
stimulusTuningCurves adpt stmttls0 stmttls1 mnrn = execEC $ do

    goalLayoutLR
    radiansAbscissaLR

    let preln = stimulusSpikeRate mnrn stmttls0
        pstln = stimulusSpikeRate mnrn stmttls1

    let adptpi0 = 2*pi*adpt/360
        adptpi = 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

    let dffln = zipWith (\(x1,y1) (_,y2) -> (x1,y1 - y2)) preln pstln

    plotRight . return $ vlinePlot "adaptor" (solidLine 4 $ opaque black) adptpi

    plotRight . liftEC $ do
        plot_lines_values .= [[(0,0),(1,0)]]
        plot_lines_style .= solidLine 4 (opaque white)

    plotLeft . liftEC $ do
        plot_lines_values .= [[(0,0),(1,0)]]
        plot_lines_style .= solidLine 4 (opaque white)

    plotRight . liftEC $ do
        plot_lines_title .= "diff"
        plot_lines_values .= [loopRadiansPlotData dffln]
        plot_lines_style .= dashedLine 4 [4,4] (opaque blue)

    plotLeft . liftEC $ do
        plot_lines_title .= "pre"
        plot_lines_values .= [loopRadiansPlotData preln]
        plot_lines_style .= solidLine 4 (opaque black)

    plotLeft . liftEC $ do
        plot_lines_title .= "post"
        plot_lines_values .= [loopRadiansPlotData pstln]
        plot_lines_style .= solidLine 4 (opaque red)

stimulusFanoFactor
    :: Double
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> Maybe NeuronID
    -> Layout Double Double
stimulusFanoFactor adpt stmstrm0 stmstrm1 mnrn = execEC $ do

    let grps0 = streamToSpikeGroups mnrn stmstrm0
        grps1 = streamToSpikeGroups mnrn stmstrm1

    let ffs0 = estimateFanoFactor . fmap snd <$> grps0
        ffs1 = estimateFanoFactor . fmap snd <$> grps1

    let adptpi = adaptorToRads adpt

    goalLayout
    radiansAbscissa

    layout_title .= "Fano Factors"

    plot . liftEC $ do

        plot_lines_title .= "pre"
        plot_lines_style .= solidLine 3 (opaque blue)
        plot_lines_values .= [ loopRadiansPlotData $ zip (fst . head <$> grps0) ffs0 ]

    plot . liftEC $ do

        plot_lines_title .= "post"
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [ loopRadiansPlotData $ zip (fst . head <$> grps1) ffs1 ]

    plot . return $ vlinePlot "" (solidLine 3 $ opaque black) adptpi

rawStimulusSpikes
    :: forall t
    . (KnownNat t, 1 <= t)
    => Proxy t
    -> Double
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> (Layout Double Double, Layout Int Double)
rawStimulusSpikes _ adpt stmstrm =

    let adptpi = adaptorToRads adpt
        spks = streamToSpikeCounts Nothing stmstrm
        (xs,ys) = B.unzip . fromJust $ (B.fromList spks :: Maybe (B.Vector t (Stimulus, Int)))
        --(lm0,aic0,bic0) = fourierFit sinusoid0 xs ys
        ys' = realToFrac <$> ys
        nrm0 = mle ys' :: Source # Normal
        aic0 = akaikesInformationCriterion nrm0 ys'
        bic0 = bayesianInformationCriterion nrm0 ys'
        rng = range 0 (2*pi) 100
        mu = fst . S.toPair $ coordinates nrm0
        ln0 = zip rng $ repeat mu
        (lm1,aic1,bic1) = fourierFit sinusoid1 xs ys
        (lm2,aic2,bic2) = fourierFit sinusoid2 xs ys
        (_,aic3,bic3) = fourierFit sinusoid3 xs ys
        (_,aic4,bic4) = fourierFit sinusoid4 xs ys
        --(_,aic5,bic5) = fourierFit sinusoid5 xs ys
        (_,ln1,_) = fourierFitToLines sinusoid1 lm1
        (_,ln2,_) = fourierFitToLines sinusoid2 lm2

        grps = streamToSpikeGroups Nothing stmstrm
        mus = average . map (realToFrac . snd) <$> grps

        spklyt = execEC $ do

            goalLayout
            radiansAbscissa
            layout_x_axis . laxis_title .= "Stimulus"
            layout_y_axis . laxis_title .= "Spike Count"

            plot . return $ vlinePlot "adaptor" (solidLine 2 $ opaque black) adptpi

            plot . liftEC $ do

                plot_points_style .= filledCircles 4 (opaque black)
                plot_points_values .= zip (fst . head <$> grps) mus

            plot . liftEC $ do

                plot_points_title .= "spike totals"
                plot_points_style .= filledCircles 1 (black `withOpacity` 0.5)
                plot_points_values .= [ (x,realToFrac k) | (x,k) <- spks ]

            plot . liftEC $ do

                plot_lines_title .= "alt. fits"
                plot_lines_style .= solidLine 2 (green `withOpacity` 0.6)
                plot_lines_values .= [ln0]

            plot . liftEC $ do

                plot_lines_title .= "alt. fits"
                plot_lines_style .= solidLine 2 (red `withOpacity` 0.5)
                plot_lines_values .= [ln2]

            plot . liftEC $ do

                plot_lines_title .= "theory fit"
                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [ln1]

        iclyt = execEC $ do

            goalLayout

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque green)
                plot_lines_values .= [zip [0,2..] [aic0,aic1,aic2,aic3,aic4]]
                plot_lines_title .= "aic"

            plot . liftEC $ do

                plot_lines_style .= solidLine 3 (opaque blue)
                plot_lines_values .= [zip [0,2..] [bic0,bic1,bic2,bic3,bic4]]
                plot_lines_title .= "bic"

        in (spklyt,iclyt)


fanoFactorScatter
    :: Double
    -> M.Map Stimulus (Int,M.Map NeuronID [SpikeTime])
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> [(Stimulus,M.Map NeuronID [SpikeTime])]
    -> Layout Double Double
fanoFactorScatter adpt stmttls0 stmstrm0 stmstrm1 = execEC $ do

    let adptpi00 = 2*pi*(adpt+90)/360
        adptpi0 = toPi adptpi00
        adptpi = 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

        spks0 = fmap length . snd <$> filter (\(stm,_) -> roundSD 1 stm == roundSD 1 adptpi) stmstrm0
        spks1 = fmap length . snd <$> filter (\(stm,_) -> roundSD 1 stm == roundSD 1 adptpi) stmstrm1
        nrns = M.keys $ head spks0
        ffpss = [ ( (estimateFanoFactor $ (M.! nrn) <$> spks0, estimateFanoFactor $ (M.! nrn) <$> spks1)
                  , preferredStimulus nrn stmttls0 ) | nrn <- nrns ]

    goalLayout

    layout_title .= "Fano Factor Scatter"
    layout_x_axis . laxis_title .= "Pre-Adapt"
    layout_y_axis . laxis_title .= "Post-Adapt"

    let adaptorDistance ps =
            let d0 = (adptpi - ps) / (2*pi)
             in if d0 < 0 then d0 + 1 else d0

    sequence_ $ do

        (ff,ps) <- ffpss

        return . plot . liftEC $ do

            plot_points_style .= hollowCircles 5 4 (black `withOpacity` (1 - 0.9 * adaptorDistance ps))
            plot_points_values .= [ff]

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [[ (0,0), (10,10)]]



--- Main ---


analyzeData
    :: forall nn t1 t2
    . (KnownNat nn, KnownNat t1, KnownNat t2, 1 <= nn, 1 <= t1, 1 <= t2)
    => KohnExperiment nn t1 t2
    -> IO ()
analyzeData kxp = do

    let kpdr = kohnProjectPath kxp

    (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm0"
    (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm1"
    (stmttls0 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "stmttls0"
    (stmttls1 :: M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "stmttls1"

    let nrns = M.keys . snd . (!! 2) $ M.elems stmttls0

    renderNeuralLayouts kpdr "stimulus-tuning-curves" nrns $ stimulusTuningCurves 90 stmttls0 stmttls1

    let (prespklyt,preiclyt) = rawStimulusSpikes (Proxy :: Proxy t1) 90 stmstrm0
        (pstspklyt,psticlyt) = rawStimulusSpikes (Proxy :: Proxy t2) 90 stmstrm1
        prettl = "Pre-Adapt; Dataset: " ++ experiment kxp
        pstttl = "Post-Adapt; Dataset: " ++ experiment kxp

    goalRenderableToPDF kpdr ("stimulus-rate-pre" ++ experiment kxp) 400 200 . toRenderable
        . (layout_title .~ prettl) $ prespklyt

    goalRenderableToPDF kpdr ("stimulus-rate-post" ++ experiment kxp) 400 200 . toRenderable
        . (layout_title .~ pstttl) $ pstspklyt

    goalRenderableToPDF kpdr ("information-criteria-pre" ++ experiment kxp) 400 200 . toRenderable
        . (layout_title .~ prettl) $ preiclyt
    goalRenderableToPDF kpdr ("information-criteria-post" ++ experiment kxp) 400 200 . toRenderable
        . (layout_title .~ prettl) $ psticlyt
--    goalRenderableToPDF kpdr "fano-factor-scatter" 800 800 . toRenderable
--        $ fanoFactorScatter 90 stmttls0 stmstrm0 stmstrm1

main :: IO ()
main = do
    analyzeData experiment112l44
    analyzeData experiment112l45
    analyzeData experiment112r35
    analyzeData experiment112r36
    analyzeData experiment105r62
    analyzeData experiment107l114
    analyzeData experiment112l16
    analyzeData experiment112r32
    --analyzeData big40Pooled
    analyzeData small40Pooled

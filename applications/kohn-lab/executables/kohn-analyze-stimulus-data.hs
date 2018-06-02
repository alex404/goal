{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import KohnLab

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import Data.List

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


stimulusSpikeCount
    :: Maybe NeuronID
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> [(Stimulus,Int)]
stimulusSpikeCount (Just nrn) stmttls = [ (stm, length $ nmp M.! nrn) | (stm,nmp) <- M.toList stmttls ]
stimulusSpikeCount Nothing stmttls = [ (stm, sum . M.elems $ length <$> nmp) | (stm,nmp) <- M.toList stmttls ]


--- Plots ---


-- Raw Data --

stimulusTuningCurves
    :: Double
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
    -> Maybe NeuronID
    -> LayoutLR Double Int Int
stimulusTuningCurves adpt stmttls0 stmttls1 mnrn = execEC $ do

    goalLayoutLR
    radiansAbscissaLR

    let preln = stimulusSpikeCount mnrn stmttls0
        pstln = stimulusSpikeCount mnrn stmttls1

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
stimulusFanoFactor adpt stmttls0 stmttls1 mnrn = execEC $ do

    let stmspks0 = streamToSpikeCounts mnrn stmttls0
        stmspks1 = streamToSpikeCounts mnrn stmttls1

    let grps0 = groupBy (\(x1,_) (x2,_) -> x1 == x2) $ sortBy (comparing fst) stmspks0
        grps1 = groupBy (\(x1,_) (x2,_) -> x1 == x2) $ sortBy (comparing fst) stmspks1
        ffs0 = estimateFanoFactor . fmap snd <$> grps0
        ffs1 = estimateFanoFactor . fmap snd <$> grps1

    let adptpi0 = 2*pi*adpt/360
        adptpi = 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

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

    plot . return $ vlinePlot "" (solidLine 4 $ opaque black) adptpi

fanoFactorScatter
    :: Double
    -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
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


fitData
    :: forall nn t1 t2
    . (KnownNat nn, KnownNat t1, KnownNat t2, 1 <= nn, 1 <= t1, 1 <= t2)
    => KohnExperiment nn t1 t2
    -> IO ()
fitData kxp = do

    let kpdr = kohnProjectPath kxp

    adpt <- getAdaptor kxp

    (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm0"
    (stmstrm1 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm1"
    (stmttls0 :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "stmttls0"
    (stmttls1 :: M.Map Stimulus (M.Map NeuronID [SpikeTime])) <- read <$> goalReadFile kpdr "stmttls1"

    let nrns = M.keys . (!! 2) $ M.elems stmttls0

    --renderNeuralLayouts kpdr "stimulus-total-spikes" nrns $ stimulusTuningCurves adpt stmttls0 stmttls1

    --renderNeuralLayouts kpdr "raw-fano-factors" nrns $ stimulusFanoFactor adpt stmstrm0 stmstrm1

    goalRenderableToSVG kpdr "fano-factor-scatter" 800 800 . toRenderable
        $ fanoFactorScatter adpt stmttls0 stmstrm0 stmstrm1

main :: IO ()
main = do
    fitData experiment112l44
    fitData experiment112l45
    fitData experiment112r35
    fitData experiment112r36
    fitData experiment105r62
    fitData experiment107l114
    fitData experiment112l16
    fitData experiment112r32


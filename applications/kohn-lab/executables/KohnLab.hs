{-# LANGUAGE KindSignatures,DataKinds,TypeOperators,TypeFamilies #-}

module KohnLab
    (
    -- * Experiments
      KohnExperiment (KohnExperiment,protocol,experiment)
    , kohnProjectPath
    , experiment112l44
    , experiment112l45
    , experiment112r35
    , experiment112r36
    , experiment105r62
    , experiment107l114
    , experiment112l16
    , experiment112r29
    , experiment112r32
    -- * Functions
    , blockStream
    , blockToStimulusStream
    , blockIDTotals
    , blockIDToStimulusTotals
    -- ** Type Synonyms
    , BlockID
    , BlockEvent
    , NeuronID
    , SpikeTime
    , SpikeInterval
    , Stimulus
    , NStimuli
    -- ** Von Mises
    , vonMisesFits
    -- * IO
    , getBIDs
    , getSpikes
    , getChannels
    , getAdaptor
    -- * Miscellaneous
    , sinusoid
    , pltsmps
    , renderNeuralLayouts
    , streamToSpikeCounts
    , preferredStimulus
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

-- Other --

import qualified Data.Map as M
import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V
import Data.List


--- Experiments ---


-- Type Synonyms --

type BlockID = Int
type BlockEvent = (BlockID,SpikeTime)
type NeuronID = (Int,Int)
type Stimulus = Double
type SpikeTime = Double
type SpikeInterval = Double
type NStimuli = 8


-- Kohn Lab Record --

data KohnExperiment (nn :: Nat) (t0 :: Nat) (t1 :: Nat) = KohnExperiment
    { protocol :: String
    , experiment :: String }

-- Experiments --

kohnProjectPath :: KohnExperiment nn t0 t1 -> FilePath
kohnProjectPath kd = "kohn-data/" ++ protocol kd ++ "/" ++ experiment kd

experiment112l44 :: KohnExperiment 55 400 320
experiment112l44 = KohnExperiment "small40" "112l44"

experiment112l45 :: KohnExperiment 42 400 320
experiment112l45 = KohnExperiment "small40" "112l45"

experiment112r35 :: KohnExperiment 11 400 320
experiment112r35 = KohnExperiment "small40" "112r35"

experiment112r36 :: KohnExperiment 13 400 320
experiment112r36 = KohnExperiment "small40" "112r36"

experiment105r62 :: KohnExperiment 81 211 240
experiment105r62 = KohnExperiment "big40" "105r62"

experiment107l114 :: KohnExperiment 126 211 240
experiment107l114 = KohnExperiment "big40" "107l114"

experiment112l16 :: KohnExperiment 118 400 320
experiment112l16 = KohnExperiment "big40" "112l16"

experiment112r29 :: KohnExperiment 121 400 320
experiment112r29 = KohnExperiment "big40" "112r29"

experiment112r32 :: KohnExperiment 126 400 320
experiment112r32 = KohnExperiment "big40" "112r32"


--- Functions ---


-- Sample Stream Builder --

nullNeuronMap :: Monoid a => [(Int,Int,Double)] -> M.Map NeuronID a
nullNeuronMap ecss =
    let ecs = [ (e,c) | (e,c,_) <- ecss, e /= 2000 ]
     in M.fromList . zip (nub ecs) $ repeat mempty

nullBlockIDMap :: Monoid a => [BlockID] -> M.Map BlockID a
nullBlockIDMap bids =
     M.fromList . zip bids $ repeat mempty

-- | Filter out weird channels, and cut off the initial (meaningless?) part of the data stream.
preFilterSpikes :: Maybe [Int] -> [(Int,Int,Double)] -> [(Int,Int,Double)]
preFilterSpikes mchns ecss =
    let ecss' = dropWhile (\(e,_,_) -> e /= 2000) ecss
     in filter (\(e,c,_) -> e == 2000 || (maybe True (elem e) mchns && 0 < c && c < 10)) ecss'

-- | Breaks the data stream into a stream of pairs of blockIDs + times, and the stream of spikes
-- that happened in that block.
breakSpikeBlocks :: [BlockID] -> [(Int,Int,Double)] -> [(BlockEvent,[(NeuronID,[SpikeTime])])]
breakSpikeBlocks [] _ = []
breakSpikeBlocks (bid:bids) ((_,_,s0):ecss) =
    let (ecss1,ecss') = span (\(e,_,_) -> e /= 2000) ecss
     in ( (bid,s0), [ ((e,c),[s]) | (e,c,s) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: Maybe [Int] -> [Int] -> [(Int,Int,Double)] -> [(BlockEvent,M.Map NeuronID [SpikeTime])]
blockStream mchns bids ecss =
    let ecss' = preFilterSpikes mchns ecss
        mp0 = nullNeuronMap ecss'
        bnspks = breakSpikeBlocks bids ecss'
     in [ (bspk, flip M.union mp0 $ M.fromListWith (++) nspks) | (bspk,nspks) <- bnspks ]

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockToStimulusStream :: [(BlockEvent,M.Map NeuronID [SpikeTime])] -> [(Stimulus,M.Map NeuronID [SpikeTime])]
blockToStimulusStream = mapMaybe patternMatch
    where patternMatch ((k,tm),mp)
            | k == 0 || k == 1 = Nothing
            | k <= 9 = Just ((fromIntegral k-2) * 2*pi/8,mp)
            | k <= 17 = patternMatch ((k-8,tm),mp)
            | otherwise = Nothing

blockIDTotals :: [BlockID] -> [(BlockEvent,M.Map NeuronID [SpikeTime])] -> M.Map BlockID (M.Map NeuronID [SpikeTime])
blockIDTotals bids bstrm =
    let bstrm' = [(bid,nmp) | ((bid,_),nmp) <- bstrm]
        --n = fromIntegral . last . map length . group . sort $ fst . fst <$> allbstrm
     in flip M.union (nullBlockIDMap bids) $ M.fromListWith (M.unionWith (++)) bstrm'

blockIDToStimulusTotals :: M.Map BlockID (M.Map NeuronID [SpikeTime]) -> M.Map Stimulus (M.Map NeuronID [SpikeTime])
blockIDToStimulusTotals bidmp =
    let bidstms = zip [2..9] $ range 0 (2*pi) 9
     in foldr foldfun mempty bidstms
    where foldfun (bid,stm) = M.insert stm (M.unionWith (++) (bidmp M.! bid) (bidmp M.! (bid + 8)))


--- Plots ---


-- Tuning Curves --

-- Fit --

pltsmps :: S.Vector 100 Double
pltsmps = S.range 0 (2*pi)

sumOfTuningCurves
    :: (KnownNat k, KnownNat j)
    => B.Vector k (B.Vector j (Double, Double))
    -> B.Vector j (Double, Double)
sumOfTuningCurves tcs =
    let tcs' = B.unzip <$> tcs
        xs' = fst . head $ B.toList tcs'
        cs' = sum $ snd <$> tcs'
     in B.zip xs' cs'

sinusoid :: Double -> S.Vector 3 Double
sinusoid x =
     fromJust $ S.fromList [cos x,sin x,1]

streamToSpikeCounts
    :: Maybe NeuronID
    -> [(x,M.Map NeuronID [SpikeTime])]
    -> [(x,Int)]
streamToSpikeCounts (Just nrn) strm = [ (x, length $ nmp M.! nrn) | (x,nmp) <- strm ]
streamToSpikeCounts Nothing strm = [ (x, sum $ length <$> M.elems nmp) | (x,nmp) <- strm ]

renderNeuralLayouts :: ToRenderable a => FilePath -> String -> [NeuronID] -> (Maybe NeuronID -> a) -> IO ()
renderNeuralLayouts pdr flnm nrns rf = do
    goalRenderableToSVG pdr flnm 1200 800 . toRenderable $ rf Nothing
    sequence_ $ do
        nrn <- nrns
        return . goalRenderableToSVG
            (pdr ++ "/neuron-" ++ show nrn) flnm 1200 800 . toRenderable $ rf (Just nrn)

preferredStimulus :: NeuronID -> M.Map Stimulus (M.Map NeuronID [SpikeTime]) -> Stimulus
preferredStimulus nrn stmttls =
    fst . maximumBy (comparing snd) . M.toList $ length . (M.! nrn) <$> stmttls

vonMisesFits :: KnownNat nn
      => (Mean ~> Natural # R nn Poisson <* VonMises)
      -> Double
      -> LayoutLR Double Double Double
vonMisesFits lkl adpt = execEC $ do

    let tcs = tuningCurves (G.convert pltsmps) lkl
        stcs = sumOfTuningCurves tcs
        cs0 = snd <$> stcs
        cs = G.convert cs0
        spltsmps = G.convert pltsmps
        alphs = linearLeastSquares (S.map sinusoid spltsmps) cs
        sinsmps = S.map (S.dotProduct alphs . sinusoid) spltsmps
        r2 = rSquared cs sinsmps
        amp = let x1:x2:_ = S.toList alphs
               in sqrt $ square x1 + square x2

    layoutlr_title .= ("Amplitude: " ++ take 4 (show amp) ++ ", r^2: " ++ show (roundSD 3 r2))

    goalLayoutLR
    radiansAbscissaLR

    layoutlr_x_axis . laxis_title .= "Stimulus"
    layoutlr_left_axis . laxis_title .= "Firing Rate"
    layoutlr_right_axis . laxis_title .= "Total Rate"
    layoutlr_left_axis . laxis_generate .= scaledAxis def (0,300)
    layoutlr_right_axis . laxis_generate .= scaledAxis def (0,5000)

    let adptpi0 = 2*pi*adpt/360
        adptpi =
            2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

    plotRight . return $ vlinePlot "" (solidLine 4 $ opaque black) adptpi

    plotLeft . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= toList (toList <$> tcs)

    plotRight . liftEC $ do
        plot_lines_style .= solidLine 3 (opaque black)
        plot_lines_values .= [ toList stcs ]

    plotRight . liftEC $ do
        plot_lines_style .= dashedLine 3 [8,10] (opaque blue)
        plot_lines_values .= [ zip (S.toList pltsmps) (S.toList sinsmps) ]


--- IO ---


getBIDs :: KohnExperiment nn t0 t1 -> IO [Int]
getBIDs kxp = do

    let dr = "adaptation/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    bidstr <- readFile $ csvdr ++ "/blockIDs.csv"
    return $ read <$> lines bidstr

getSpikes :: KohnExperiment nn t0 t1 -> IO [(Int,Int,Double)]
getSpikes kxp = do

    let dr = "adaptation/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    ecsstr <- BS.readFile $ csvdr ++ "/spikes.csv"
    let (Right ecssV) = C.decode C.NoHeader ecsstr
    return $ V.toList ecssV

getChannels :: KohnExperiment nn t0 t1 -> IO (Maybe [Int])
getChannels kxp = do

    let dr = "adaptation/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    bl <- doesFileExist $ csvdr ++ "/channels.csv"

    if bl
       then do
           chnstr <- readFile $ csvdr ++ "/channels.csv"
           return . Just . map read $ lines chnstr
       else return Nothing

getAdaptor :: KohnExperiment nn t0 t1 -> IO Double
getAdaptor kxp = do

    let dr = "adaptation/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    adpstr <- readFile $ csvdr ++ "/adaptor.csv"
    return . head $ read <$> lines adpstr
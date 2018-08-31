{-# LANGUAGE TupleSections,KindSignatures,DataKinds,TypeOperators,TypeFamilies #-}

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
    , big40Pooled
    , small40Pooled
    -- , experiment112r29
    , experiment112r32
    -- * Functions
    , blockStream
    , blockToStimulusStream
    , blockIDTotals
    , blockIDToStimulusTotals
    , converter
    -- ** Type Synonyms
    , BlockID
    , BlockEvent
    , NeuronID
    , SpikeTime
    , SpikeInterval
    , Stimulus
    , NStimuli
    , toPooledNeuronID
    -- * IO
    , getBIDs
    , getSpikes
    , getChannels
    , getAdaptor
    -- * Miscellaneous
    , adaptorToRads
    , bidToCentredStimulus
    , sinusoid1
    , sinusoid2
    , sinusoid3
    , sinusoid4
    , sinusoid5
    , pltsmps
    , renderNeuralLayouts
    , streamToSpikeCounts
    , streamToSpikeGroups
    , preferredStimulus
    ) where


--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

-- Other --

import qualified Data.Map as M
import qualified Data.Csv as C
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V
import qualified Data.Vector as VB
import Data.List


--- Experiments ---


-- Type Synonyms --

type BlockID = Int
type BlockEvent = (BlockID,SpikeTime)
type Stimulus = Double
type SpikeTime = Double
type SpikeInterval = Double
type NStimuli = 8

-- Types --

data NeuronID =
    NeuronID (Int,Int) | PooledNeuronID (String, Int, Int)
    deriving (Eq, Read, Show, Ord)

toPooledNeuronID
    :: String
    -> NeuronID
    -> NeuronID
toPooledNeuronID _ (PooledNeuronID _) = error "Why are you trying to pool a pooled neuron?"
toPooledNeuronID ex (NeuronID (e,c)) = PooledNeuronID (ex,e,c)

-- Kohn Lab Record --

data KohnExperiment (nn :: Nat) = KohnExperiment
    { protocol :: String
    , experiment :: String }

-- Experiments --

kohnProjectPath :: KohnExperiment nn -> FilePath
kohnProjectPath kd = "patterson-2013/" ++ protocol kd ++ "/" ++ experiment kd

experiment112l44 :: KohnExperiment 55
experiment112l44 = KohnExperiment "small40" "112l44"

experiment112l45 :: KohnExperiment 42
experiment112l45 = KohnExperiment "small40" "112l45"

experiment112r35 :: KohnExperiment 11
experiment112r35 = KohnExperiment "small40" "112r35"

experiment112r36 :: KohnExperiment 13
experiment112r36 = KohnExperiment "small40" "112r36"

experiment105r62 :: KohnExperiment 81
experiment105r62 = KohnExperiment "big40" "105r62"

experiment107l114 :: KohnExperiment 126
experiment107l114 = KohnExperiment "big40" "107l114"

experiment112l16 :: KohnExperiment 118
experiment112l16 = KohnExperiment "big40" "112l16"

--experiment112r29 :: KohnExperiment 121 400 320
--experiment112r29 = KohnExperiment "big40" "112r29"

experiment112r32 :: KohnExperiment 126
experiment112r32 = KohnExperiment "big40" "112r32"

big40Pooled :: KohnExperiment (81+126+118+126)
big40Pooled = KohnExperiment "big40" "pooled"

small40Pooled :: KohnExperiment (55+42+11+13)
small40Pooled = KohnExperiment "small40" "pooled"

--experiment112l44 :: KohnExperiment 117 400 320
--experiment112l44 = KohnExperiment "small40" "112l44"
--
--experiment112l45 :: KohnExperiment 42 400 320
--experiment112l45 = KohnExperiment "small40" "112l45"
--
--experiment112r35 :: KohnExperiment 40 400 320
--experiment112r35 = KohnExperiment "small40" "112r35"
--
--experiment112r36 :: KohnExperiment 13 400 320
--experiment112r36 = KohnExperiment "small40" "112r36"
--
--experiment105r62 :: KohnExperiment 81 211 240
--experiment105r62 = KohnExperiment "big40" "105r62"
--
--experiment107l114 :: KohnExperiment 126 211 240
--experiment107l114 = KohnExperiment "big40" "107l114"
--
--experiment112l16 :: KohnExperiment 118 400 320
--experiment112l16 = KohnExperiment "big40" "112l16"
--
----experiment112r29 :: KohnExperiment 121 400 320
----experiment112r29 = KohnExperiment "big40" "112r29"
--
--experiment112r32 :: KohnExperiment 126 400 320
--experiment112r32 = KohnExperiment "big40" "112r32"
--
--big40Pooled :: KohnExperiment (81+126+118+126) (2*400 + 2*211) (2*320 + 2*240)
--big40Pooled = KohnExperiment "big40" "pooled"
--
--small40Pooled :: KohnExperiment (117+42+40+13) 400 320
--small40Pooled = KohnExperiment "small40" "pooled"

--- Functions ---

converter :: KnownNat nn => (Stimulus,M.Map NeuronID [SpikeTime]) -> (B.Vector nn Int,Stimulus)
converter (s,mp) =
    (fromJust . B.toSized . VB.fromList $ length <$> M.elems mp,s)



-- Sample Stream Builder --

nullNeuronMap :: Monoid a => [(Int,Int,Double)] -> M.Map NeuronID a
nullNeuronMap ecss =
    let ecs = [ NeuronID (e,c) | (e,c,_) <- ecss, e /= 2000 ]
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
     in ( (bid,s0), [ (NeuronID (e,c),[s]) | (e,c,s) <- ecss1 ] ) : breakSpikeBlocks bids ecss'
breakSpikeBlocks _ _ = error "Block ID misalignment in breakSpikeBlocks"

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockStream :: Maybe [Int] -> [Int] -> [(Int,Int,Double)] -> [(BlockEvent,M.Map NeuronID [SpikeTime])]
blockStream mchns bids ecss =
    let ecss' = preFilterSpikes mchns ecss
        mp0 = nullNeuronMap ecss'
        bnspks = breakSpikeBlocks bids ecss'
     in [ (bspk, flip M.union mp0 $ M.fromListWith (++) nspks) | (bspk,nspks) <- bnspks ]

-- | Combines preFilterSpikes and breakSpikeBlocks.
blockToStimulusStream :: Double -> [(BlockEvent,M.Map NeuronID [SpikeTime])] -> [(Stimulus,M.Map NeuronID [SpikeTime])]
blockToStimulusStream adpt = mapMaybe patternMatch
    where patternMatch ((k,tm),mp)
            | k == 0 || k == 1 = Nothing
            | k <= 9 = Just (bidToCentredStimulus adpt k,mp)
            | k <= 17 = patternMatch ((k-8,tm),mp)
            | otherwise = Nothing

blockIDTotals :: [BlockID] -> [(BlockEvent,M.Map NeuronID [SpikeTime])] -> M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
blockIDTotals bids bstrm =
    let bstrm' = [(bid,(1,nmp)) | ((bid,_),nmp) <- bstrm]
        --n = fromIntegral . last . map length . group . sort $ fst . fst <$> allbstrm
     in flip M.union ((0,) <$> nullBlockIDMap bids)
            $ M.fromListWith (\(k1,nmp1) (k2,nmp2) -> (k1 + k2, M.unionWith (++) nmp1 nmp2)) bstrm'

blockIDToStimulusTotals
    :: Double
    -> M.Map BlockID (Int, M.Map NeuronID [SpikeTime])
    -> M.Map Stimulus (Int, M.Map NeuronID [SpikeTime])
blockIDToStimulusTotals adpt bidmp =
    let bidstms = zip [2..9] . traceGiven $ bidToCentredStimulus adpt <$> [2..9]
     in foldr foldfun mempty bidstms
    where foldfun (bid,stm) =
              let (k1,nmp1) = bidmp M.! bid
                  (k2,nmp2) = bidmp M.! (bid + 8)
               in  M.insert stm (k1 + k2, M.unionWith (++) nmp1 nmp2)


--- Plots ---


-- Tuning Curves --

-- Fit --

adaptorToRads :: Double -> Double
adaptorToRads adpt =
    let adptpi0 = 2*pi*adpt/360
     in 2*(if adptpi0 > pi then adptpi0 - pi else adptpi0)

bidToCentredStimulus :: Double -> BlockID -> Stimulus
bidToCentredStimulus adpt bid =
    let adpt' = circulate . round $ adpt / 22.5
        bid' = circulate $ bid - 2
     in ((pi/4)*) . fromIntegral . circulate $ bid' - adpt' + 4
    where circulate k
            | k < 0 = circulate $ k + 8
            | k > 7 = circulate $ k - 8
            | otherwise = k

pltsmps :: S.Vector 100 Double
pltsmps = S.range 0 (2*pi)

sinusoid1 :: Double -> B.Vector 2 Double
sinusoid1 x =
    fromJust $ B.fromList [cos x,sin x]

sinusoid2 :: Double -> B.Vector 4 Double
sinusoid2 x =
    fromJust $ B.fromList [cos x,sin x, cos (2*x), sin (2*x)]

sinusoid3 :: Double -> B.Vector 6 Double
sinusoid3 x =
    fromJust $ B.fromList [cos x,sin x, cos (2*x), sin (2*x), cos (3*x), sin (3*x)]

sinusoid4 :: Double -> B.Vector 8 Double
sinusoid4 x =
    fromJust $ B.fromList [cos x,sin x, cos (2*x), sin (2*x), cos (3*x), sin (3*x), cos (4*x), sin (4*x)]

sinusoid5 :: Double -> B.Vector 10 Double
sinusoid5 x =
    fromJust $ B.fromList [cos x,sin x, cos (2*x), sin (2*x), cos (3*x), sin (3*x)
      , cos (4*x), sin (4*x), cos (5*x), sin (5*x)]

streamToSpikeCounts
    :: Maybe NeuronID
    -> [(x,M.Map NeuronID [SpikeTime])]
    -> [(x,Int)]
streamToSpikeCounts (Just nrn) strm = [ (x, length $ nmp M.! nrn) | (x,nmp) <- strm ]
streamToSpikeCounts Nothing strm = [ (x, sum $ length <$> M.elems nmp) | (x,nmp) <- strm ]

streamToSpikeGroups
    :: Ord x
    => Maybe NeuronID
    -> [(x,M.Map NeuronID [SpikeTime])]
    -> [[(x,Int)]]
streamToSpikeGroups mnrn strm =
    let strm' = streamToSpikeCounts mnrn strm
     in groupBy (\(x1,_) (x2,_) -> x1 == x2) $ sortOn fst strm'

renderNeuralLayouts :: ToRenderable a => FilePath -> String -> [NeuronID] -> (Maybe NeuronID -> a) -> IO ()
renderNeuralLayouts pdr flnm nrns rf = do
    goalRenderableToSVG pdr flnm 1200 800 . toRenderable $ rf Nothing
    sequence_ $ do
        nrn0 <- nrns
        case nrn0 of
          NeuronID nrn ->
              return $ goalRenderableToSVG (pdr ++ "/neuron-" ++ show nrn) flnm 1200 800 . toRenderable $ rf (Just nrn0)
          _ -> return $ return ()

preferredStimulus :: NeuronID -> M.Map Stimulus (Int, M.Map NeuronID [SpikeTime]) -> Stimulus
preferredStimulus nrn stmttls =
    fst . maximumBy (comparing snd) . M.toList $ length . (M.! nrn) . snd <$> stmttls

--fourierFit
--    :: (1 <= t, KnownNat t, KnownNat k)
--    => (Double -> B.Vector k Double)
--    -> Sample t Normal
--    -> Sample t Poisson
--    -> (Mean #> Source # LinearModel Normal (Replicated k StandardNormal), Double, Double)
--fourierFit f xs0 ys0 =
--    let xs = f <$> xs0
--        ys = realToFrac <$> ys0
--        lm = fitLinearModel xs ys
--        aic = conditionalAkaikesInformationCriterion lm xs ys
--        bic = conditionalBayesianInformationCriterion lm xs ys
--     in (lm,roundSD 2 aic,roundSD 2 bic)
--
--fourierFitToLines
--    :: KnownNat k
--    => (Double -> B.Vector k Double)
--    -> Mean #> Source # LinearModel Normal (Replicated k StandardNormal)
--    -> ([(Double,Double)],[(Double,Double)],[(Double,Double)])
--fourierFitToLines f lm =
--    let bxs :: B.Vector 100 Double
--        bxs = B.range 0 (2*pi)
--        pys = splitReplicated $ lm >$>* (f <$> bxs)
--        xs = B.toList bxs
--        (ysmn, ys, ysmx) = unzip3 $ splitter <$> S.toList pys
--     in (zip xs ysmn, zip xs ys, zip xs ysmx)
--    where splitter py = let (mu,vr) = S.toPair $ coordinates py
--                            sd = sqrt vr
--                         in (mu - sd, mu, mu + sd)


--- IO ---


getBIDs :: KohnExperiment nn -> IO [Int]
getBIDs kxp = do

    let dr = "patterson-2013/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    bidstr <- readFile $ csvdr ++ "/blockIDs.csv"
    return $ read <$> lines bidstr

getSpikes :: KohnExperiment nn -> IO [(Int,Int,Double)]
getSpikes kxp = do

    let dr = "patterson-2013/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    ecsstr <- BS.readFile $ csvdr ++ "/spikes.csv"
    let (Right ecssV) = C.decode C.NoHeader ecsstr
    return $ V.toList ecssV

getChannels :: KohnExperiment nn -> IO (Maybe [Int])
getChannels kxp = do

    let dr = "patterson-2013/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    bl <- doesFileExist $ csvdr ++ "/channels.csv"

    if bl
       then do
           chnstr <- readFile $ csvdr ++ "/channels.csv"
           return . Just . map read $ lines chnstr
       else return Nothing

getAdaptor :: KohnExperiment nn -> IO Double
getAdaptor kxp = do

    let dr = "patterson-2013/" ++ protocol kxp
        flnm = experiment kxp

    csvdr <- goalDatasetPath dr flnm

    adpstr <- readFile $ csvdr ++ "/adaptor.csv"
    return . head $ read <$> lines adpstr

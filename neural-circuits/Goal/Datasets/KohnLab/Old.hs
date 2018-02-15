module Goal.Datasets.KohnLab.Old where

import qualified Data.Map.Strict as M


--- Globals ---


filterChannelsAndStart :: [Int] -> [(Int,Int,Double)] -> [(Int,Int,Double)]
filterChannelsAndStart chns ecss =
    let ecss' = tail $ dropWhile (\(e,_,_) -> e /= 2000) ecss
     in filter (\(e,c,_) -> e == 2000 || (elem e chns && 0 < c && c < 10)) ecss'

tuningMapper :: [((Int,Int),[(Int,[Double])])] -> M.Map (Int,Int) (M.Map Int [Double])
tuningMapper ecbss = M.fromListWith (++) <$> M.fromListWith (++) ecbss

mapListConstructor :: [Int] -> [(Int,Int,Double)] -> [((Int,Int),[(Int,[Double])])]
mapListConstructor [] _ = []
mapListConstructor (bid:bids) ecss =
    let (ecss1,ecss') = span (\(e,_,_) -> e /= 2000) ecss
     in [ ((e,c),[(bid,[s])]) | (e,c,s) <- ecss1 ] ++ mapListConstructor bids (tail ecss')

blockTimes :: [Int] -> [(Int,Int,Double)] -> [(Int,Double)]
blockTimes bids ecss = zip bids [ s | (e,_,s) <- ecss, e == 2000 ]

blockTimeTotals :: [(Int,Double)] -> M.Map Int Double
blockTimeTotals btms =
    let btms' = zipWith (\(b1,t1) (_,t2) -> (b1,t2-t1)) btms $ tail btms
     in M.fromListWith (+) btms'

tuningMap :: [Int] -> [Int] -> [(Int,Int,Double)] -> M.Map (Int,Int) (M.Map Int [Double])
tuningMap chns bids ecss =
    tuningMapper . mapListConstructor bids $ filterChannelsAndStart chns ecss

-- | First index is neuron, second index is blockID.
spikeCountHistogram :: M.Map (Int,Int) (M.Map Int [Double]) -> [[Int]]
spikeCountHistogram spkmp =
    map (length . snd) . M.toAscList . snd <$> M.toAscList spkmp

{-# LANGUAGE DataKinds,ScopedTypeVariables,FlexibleContexts #-}

import Goal.Core
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Numeric.LinearAlgebra as H
import qualified Criterion.Main as C
import qualified System.Random.MWC.Probability as P


--- Globals ---


type Rows = 1000
type Columns = 1000

n :: Int
n = 10

addMatrices :: (KnownNat n, KnownNat k) => S.Matrix n k Double -> S.Matrix n k Double -> S.Matrix n k Double
addMatrices (G.Matrix mtx1) (G.Matrix mtx2) = G.Matrix $ S.add mtx1 mtx2


--- Function ---


averageOuterProduct2 v12s =
    let ln = fromIntegral $ length v12s
     in S.withMatrix (S.scale ln) $ foldr foldfun (G.Matrix $ S.replicate 0) v12s
    where foldfun (v1,v2) mtx = addMatrices mtx $ S.outerProduct v1 v2

averageWeightedOuterProduct2 wv12s =
    foldr foldfun (G.Matrix $ S.replicate 0) wv12s
        where foldfun (w,v1,v2) mtx = addMatrices mtx . S.outerProduct v1 $ S.scale w v2

--- Main ---


main :: IO ()
main = do

    let rnd :: P.Prob IO Double
        rnd = P.uniformR (-1,1)

    v1s :: [S.Vector Rows Double]
        <- P.withSystemRandom . P.sample . replicateM n $ S.replicateM rnd
    v2s :: [S.Vector Columns Double]
        <- P.withSystemRandom . P.sample . replicateM n $ S.replicateM rnd

    let ws :: [Double]
        ws0 = [1..fromIntegral n]
        ws = (/sum ws0) <$> ws0

    let v12s = zip v1s v2s
        wv12s = zip3 ws v1s v2s
    let v11s = zip v1s v1s
        wv11s = zip3 ws v1s v1s

    C.defaultMain
       [ C.bench "average-outer-product" $ C.nf averageOuterProduct2 v12s
       , C.bench "bulk-outer-product" $ C.nf S.averageOuterProduct v12s
       , C.bench "average-weighted-outer-product" $ C.nf averageWeightedOuterProduct2 wv12s
       , C.bench "bulk-weighted-outer-product" $ C.nf S.weightedAverageOuterProduct wv12s
       , C.bench "average-outer-product0" $ C.nf averageOuterProduct2 v11s
       , C.bench "bulk-outer-product0" $ C.nf S.averageOuterProduct v11s
       , C.bench "average-weighted-outer-product0" $ C.nf averageWeightedOuterProduct2 wv11s
       , C.bench "bulk-weighted-outer-product0" $ C.nf S.weightedAverageOuterProduct wv11s ]

{-# LANGUAGE
    TypeOperators,
    DataKinds
    #-}

import Goal.Core
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Numeric.LinearAlgebra as H
import qualified Criterion.Main as C
import qualified System.Random.MWC.Probability as P



--- Globals ---


expnm :: String
expnm = "convolutions"

-- Sizes --

type KRadius = 2
type KDiameter = (2*KRadius + 1)
type KNumber = 50
type MSize = 50
type MDepth = 50

pkr :: Proxy KRadius
pkr = Proxy

pkd :: Proxy KDiameter
pkd = Proxy

pms :: Proxy MSize
pms = Proxy

kdmt,ms,krd :: Int
krd = natValInt pkr
kdmt = natValInt pkd
ms = natValInt pms

goalCorr
    :: ( S.Matrix KNumber (KDiameter*KDiameter*MDepth) Double
       , S.Matrix MDepth (MSize*MSize) Double )
    -> S.Matrix KNumber (MSize*MSize) Double
goalCorr (krns,mtx) = S.crossCorrelate2d pkr pkr pms pms krns mtx

goalConv
    :: ( S.Matrix KNumber (KDiameter*KDiameter*MDepth) Double
       , S.Matrix KNumber (MSize*MSize) Double )
    -> S.Matrix MDepth (MSize*MSize) Double
goalConv (krns,mtx) = S.convolve2d pkr pkr pms pms krns mtx

hmatrixCorr
    :: (B.Vector KNumber (B.Vector MDepth (H.Matrix Double)), B.Vector MDepth (H.Matrix Double))
    -> B.Vector KNumber (H.Matrix Double)
hmatrixCorr (krnss,mtxs) = fromJust . B.fromList
    $ [ foldr1 H.add [ H.corr2 krn mtx | (krn,mtx) <- B.toList $ B.zip krns mtxs ] | krns <- B.toList krnss ]

repadHMatrix :: H.Matrix Double -> H.Matrix Double
repadHMatrix hmtx =
    let rw' = replicate krd . H.fromList $ replicate ms 0
        hmtx' = H.fromRows . (rw' ++) . (++ rw') $ H.toRows hmtx
        cl' = replicate krd . H.fromList $ replicate (ms + 2*krd) 0
     in H.fromColumns . (cl' ++) . (++ cl') $ H.toColumns hmtx'

--depadHMatrix :: H.Matrix Double -> H.Matrix Double
--depadHMatrix hmtx =
--    let hmtx' = H.fromRows . take mn . drop kdmt $ H.toRows hmtx
--     in H.fromColumns . take mn . drop kdmt $ H.toColumns hmtx'


-- Benchmark

main :: IO ()
main = do

    let rnd :: P.Prob IO Double
        rnd = P.uniformR (-1,1)

    mtxv <- P.withSystemRandom . P.sample $ S.replicateM rnd
    krnv <- P.withSystemRandom . P.sample $ S.replicateM rnd
    mtxz <- P.withSystemRandom . P.sample $ S.replicateM rnd

    let krn = G.Matrix krnv
        mtx = G.Matrix mtxv
        crr = goalCorr (krn,mtx)

    let hmtxs = H.reshape ms . G.fromSized <$> G.convert (S.toRows mtx)
        hkrnss = fmap (H.reshape kdmt . G.fromSized . G.convert) . B.breakEvery . G.convert <$> G.convert (S.toRows krn)
        hcrr = hmatrixCorr (hkrnss,repadHMatrix <$> hmtxs)

    putStrLn "Squared Error between HMatrix and Goal solutions:"
    print . sum $
        B.zipWith (\mtx1 mtx2 -> H.sumElements . H.cmap (^(2 :: Int)) . H.add mtx1 $ H.scale (-1) mtx2)
         hcrr $ H.reshape ms . G.fromSized <$> G.convert (S.toRows crr)
    putStrLn ""

    putStrLn "Transpose Dot Product Difference:"
    let mtx' = G.Matrix mtxz
        cnv = goalConv (krn,mtx')
    print $ S.dotProduct (G.toVector crr) (G.toVector mtx')
        - S.dotProduct (G.toVector cnv) (G.toVector mtx)
    putStrLn ""

    criterionMainWithReport expnm
       [ C.bench "goal-corr" $ C.nf goalCorr (krn,mtx)
       , C.bench "goal-conv" $ C.nf goalConv (krn,mtx')
       , C.bench "hmatrix-corr" $ C.nf hmatrixCorr (hkrnss,hmtxs) ]

{-# LANGUAGE TypeOperators,DataKinds #-}

import Goal.Core
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Numeric.LinearAlgebra as H
import qualified Criterion.Main as C
import qualified System.Random.MWC.Probability as P
import qualified Data.List as L


type KRadius = 3
type KDiameter = (2*KRadius + 1)
type KNumber = 5
type MSize = 28
type MDepth = 5

pkr :: Proxy KRadius
pkr = Proxy
pkd :: Proxy KDiameter
pkd = Proxy
pkn :: Proxy KNumber
pkn = Proxy
pms :: Proxy MSize
pms = Proxy
pmd :: Proxy MDepth
pmd = Proxy

kdmt,mn,krd :: Int
krd = natValInt pkr
kdmt = natValInt pkd
mn = natValInt pms

goalCorr
    :: ( S.Vector MDepth (S.Vector KNumber (S.Matrix KDiameter KDiameter Double))
       , S.Vector MDepth (S.Matrix MSize MSize Double) )
    -> S.Vector KNumber (S.Matrix MSize MSize Double)
goalCorr (krns,mtx) = S.crossCorrelate2d krns mtx

goalConv
    :: ( S.Vector MDepth (S.Vector KNumber (S.Matrix KDiameter KDiameter Double))
       , S.Vector KNumber (S.Matrix MSize MSize Double) )
    -> S.Vector MDepth (S.Matrix MSize MSize Double)
goalConv (krns,mtx) = S.convolve2d krns mtx

goalCorr'
    :: ( S.Matrix KNumber (KDiameter*KDiameter*MDepth) Double
       , S.Vector (MSize*MSize*MDepth) Double )
    -> S.Matrix KNumber (MSize*MSize) Double
goalCorr' (krns,mtx) = S.crossCorrelate2d' pkr pkr pms pms pmd krns mtx

hmatrixCorr
    :: (B.Vector MDepth (B.Vector KNumber (H.Matrix Double)), B.Vector MDepth (H.Matrix Double))
    -> B.Vector KNumber (H.Matrix Double)
hmatrixCorr (krnss,mtxs) = fromJust . B.fromList . fmap sum
    $ L.transpose [ [ H.corr2 krn mtx | krn <- B.toList krns ] | (krns,mtx) <- B.toList $ B.zip krnss mtxs ]

repadHMatrix :: H.Matrix Double -> H.Matrix Double
repadHMatrix hmtx =
    let rw' = replicate krd . H.fromList $ replicate mn 0
        hmtx' = H.fromRows . (rw' ++) . (++ rw') $ H.toRows hmtx
        cl' = replicate krd . H.fromList $ replicate (mn + 2*krd) 0
     in H.fromColumns . (cl' ++) . (++ cl') $ H.toColumns hmtx'

depadHMatrix :: H.Matrix Double -> H.Matrix Double
depadHMatrix hmtx =
    let hmtx' = H.fromRows . take mn . drop kdmt $ H.toRows hmtx
     in H.fromColumns . take mn . drop kdmt $ H.toColumns hmtx'

kernelPrime
    :: S.Vector MDepth (S.Vector KNumber (S.Matrix KDiameter KDiameter Double))
    -> S.Matrix KNumber (KDiameter*KDiameter*MDepth) Double
kernelPrime krns =
    let krns0 = S.map (S.map S.toRows) krns
        krns' = G.toRows $ G.fromColumns krns0
        krns'' = S.map (G.toRows . G.fromColumns) krns'
        krns''' = S.map (S.map (G.toRows . G.fromColumns)) krns''
     in S.fromRows . S.map S.concat $ S.map (S.map S.concat) krns'''

matrixPrime
    :: KnownNat d
    => S.Vector d (S.Matrix MSize MSize Double)
    -> S.Vector (MSize*MSize*d) Double
matrixPrime =
    G.toVector . G.transpose . G.fromRows . S.map G.toVector

-- Benchmark

main :: IO ()
main = do

    let rnd :: P.Prob IO Double
        rnd = P.uniformR (-1,1)

    mtxv <- P.withSystemRandom . P.sample $ S.replicateM rnd
    krnv <- P.withSystemRandom . P.sample $ S.replicateM rnd
    mtxz <- P.withSystemRandom . P.sample $ S.replicateM rnd

    let mtxvs = S.breakEvery mtxv
        krnvs = S.breakEvery krnv

    let mtxs = S.map G.Matrix mtxvs
        krns = S.breakEvery $ S.map G.Matrix krnvs
        crr = goalCorr (krns,mtxs)
        krn = kernelPrime krns
        mtx = matrixPrime mtxs
        crr' = goalCorr' (krn,mtx)

    let hmtxs = S.toHMatrix <$> G.convert mtxs
        hkrns = fmap S.toHMatrix . G.convert <$> G.convert krns

    putStrLn "Squared Error between HMatrix and Goal solutions:"
    print . sum $ B.zipWith (\mtx1 mtx2 -> H.sumElements . H.cmap (^(2 :: Int)) . H.add mtx1 $ H.scale (-1) mtx2)
        (S.toHMatrix <$> G.convert crr) (hmatrixCorr (hkrns,repadHMatrix <$> hmtxs))
    putStrLn ""

    putStrLn "Squared Error between Goal and Goal' solutions:"
    print . S.sum . S.map (^2) $ G.toVector crr' - matrixPrime crr

    putStrLn "Transpose Dot Product Difference:"
    let mtxzs = S.breakEvery mtxz
        mtxs' = S.map G.Matrix mtxzs
        cnv = goalConv (krns,mtxs')
    print $ S.dotProduct (S.concatMap G.toVector crr) (S.concatMap G.toVector mtxs')
        - S.dotProduct (S.concatMap G.toVector cnv) (S.concatMap G.toVector mtxs)
    putStrLn ""

    C.defaultMain
       [ C.bench "goal-conv" $ C.nf goalCorr (krns,mtxs)
       , C.bench "goal-conv'" $ C.nf goalCorr' (krn,mtx)
       , C.bench "hmatrix-conv" $ C.nf hmatrixCorr (hkrns,hmtxs) ]

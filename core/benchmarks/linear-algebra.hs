{-# LANGUAGE DataKinds #-}

import Goal.Core
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Data.Matrix as M
import qualified Numeric.LinearAlgebra as H
import qualified Criterion.Main as C
import qualified System.Random.MWC.Probability as P


--- Globals ---


expnm :: String
expnm = "linear-algebra"

-- Sizes --

type M = 1000
type N = 10

n,m :: Int
m = 1000
n = 10

-- Matrices --

goalMatrix1 :: S.Matrix M M Double
goalMatrix1 = G.Matrix $ S.generate fromIntegral

goalMatrix2 :: S.Matrix M N Double
goalMatrix2 = G.Matrix $ S.generate fromIntegral

bGoalMatrix1 :: B.Matrix M M Double
bGoalMatrix1 = G.Matrix $ B.generate fromIntegral

bGoalMatrix2 :: B.Matrix M N Double
bGoalMatrix2 = G.Matrix $ B.generate fromIntegral

goalVal :: (S.Matrix M M Double,S.Matrix M N Double) -> Double
goalVal (m1,m2) =
    let G.Matrix v = S.matrixMatrixMultiply m1 m2
     in S.sum v

goalVal2 :: (B.Matrix M M Double, B.Matrix M N Double) -> Double
goalVal2 (m1,m2) =
    let G.Matrix v = B.matrixMatrixMultiply m1 m2
     in sum v

matrixMatrix1 :: M.Matrix Double
matrixMatrix1 = M.fromLists . take m . breakEvery m $ [0..]

matrixMatrix2 :: M.Matrix Double
matrixMatrix2 = M.fromLists . take m . breakEvery n $ [0..]

matrixVal :: (M.Matrix Double,M.Matrix Double) -> Double
matrixVal (m1,m2) =
    let m3 = M.multStd2 m1 m2
     in sum $ M.getMatrixAsVector m3

hmatrixMatrix1 :: H.Matrix Double
hmatrixMatrix1 = H.fromLists . take m . breakEvery m $ [0..]

hmatrixMatrix2 :: H.Matrix Double
hmatrixMatrix2 = H.fromLists . take m . breakEvery n $ [0..]

hmatrixVal :: (H.Matrix Double,H.Matrix Double) -> Double
hmatrixVal (m1,m2) =
    let m3 = m1 H.<> m2
     in H.sumElements m3

-- Benchmark
main :: IO ()
main = do

    let rnd :: P.Prob IO Double
        rnd = P.uniformR (-1,1)

    v1 <- P.withSystemRandom . P.sample $ S.replicateM rnd
    v2 <- P.withSystemRandom . P.sample $ S.replicateM rnd

    let m1 = G.Matrix v1
        m2 = G.Matrix v2

    let bm1 = G.Matrix $ G.convert v1
        bm2 = G.Matrix $ G.convert v2

    let m1' = M.fromLists . take m . breakEvery m $!! S.toList v1
        m2' = M.fromLists . take m . breakEvery n $!! S.toList v2

    let m1'' = H.fromLists . take m . breakEvery m $!! S.toList v1
        m2'' = H.fromLists . take m . breakEvery n $!! S.toList v2

    goalCriterionMain expnm
       [ C.bench "generative-goal" $ C.nf goalVal (goalMatrix1,goalMatrix2)
       , C.bench "generative-goal2" $ C.nf goalVal2 (bGoalMatrix1,bGoalMatrix2)
       , C.bench "generative-matrix" $ C.nf matrixVal (matrixMatrix1,matrixMatrix2)
       , C.bench "generative-hmatrix" $ C.nf hmatrixVal (hmatrixMatrix1,hmatrixMatrix2)
       , C.bench "random-goal" $ C.nf goalVal (m1,m2)
       , C.bench "random-goal2" $ C.nf goalVal2 (bm1,bm2)
       , C.bench "random-matrix" $ C.nf matrixVal (m1',m2')
       , C.bench "random-hmatrix" $ C.nf hmatrixVal (m1'',m2'') ]

-- Sanity Check
--sanityCheck :: IO ()
--sanityCheck = do
--
--    let rnd :: P.Prob IO Double
--        rnd = P.uniformR (-1,1)
--
--    putStrLn "Goal 1:"
--    print $ S.matrixMatrixMultiply goalMatrix1 goalMatrix2
--    putStrLn "Goal 2:"
--    print $ B.matrixMatrixMultiply bGoalMatrix1 bGoalMatrix2
--    putStrLn "Matrix:"
--    print $ M.multStd2 matrixMatrix1 matrixMatrix2
--    putStrLn "HMatrix:"
--    print $ hmatrixMatrix1 H.<> hmatrixMatrix2
--
--    v1 <- P.withSystemRandom . P.sample $ S.replicateM rnd
--    v2 <- P.withSystemRandom . P.sample $ S.replicateM rnd
--
--    let m1 = Matrix v1
--        m2 = Matrix v2
--
--    let bm1 = Matrix $ G.convert v1
--        bm2 = Matrix $ G.convert v2
--
--    let m1' = M.fromLists . take m . breakEvery m $!! S.toList v1
--        m2' = M.fromLists . take m . breakEvery n $!! S.toList v2
--
--    let m1'' = H.fromLists . take m . breakEvery m $!! S.toList v1
--        m2'' = H.fromLists . take m . breakEvery n $!! S.toList v2
--
--    putStrLn "Goal 1:"
--    print $ goalVal (m1,m2)
--    print $ S.matrixMatrixMultiply m1 m2
--    putStrLn "Goal 2:"
--    print $ goalVal2 (bm1,bm2)
--    print $ B.matrixMatrixMultiply bm1 bm2
--    putStrLn "Matrix:"
--    print $ matrixVal (m1',m2')
--    print $ M.multStd2 m1' m2'
--    putStrLn "HMatrix:"
--    print $ hmatrixVal (m1'',m2'')
--    print $ m1'' H.<> m2''

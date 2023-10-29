{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    TypeOperators,
    NoStarIsType,
    FlexibleContexts,
    DataKinds
    #-}

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Criterion.Main as C
import qualified System.Random.MWC as R



--- Globals ---


type M = 1000
type N = 10


-- Functions


toMatrix :: (Manifold x, Manifold y) => c # Tensor y x -> S.Matrix (Dimension y) (Dimension x) Double
toMatrix (Point xs) = G.Matrix xs


woodbury ::
    ( Cartesian # Diagonal (Euclidean M) (Euclidean M)
    , Cartesian # Tensor (Euclidean M) (Euclidean N)
    , Cartesian # Symmetric (Euclidean N) (Euclidean N) )
    -> Cartesian # Tensor (Euclidean M) (Euclidean M)
woodbury (tl,tr,br) =
    let shr = inverseSchurComplement (toTensor br) tr tl
     in woodburyMatrix tl tr shr

schur ::
    ( Cartesian # Diagonal (Euclidean M) (Euclidean M)
    , Cartesian # Tensor (Euclidean M) (Euclidean N)
    , Cartesian # Symmetric (Euclidean N) (Euclidean N) )
    -> Cartesian # Tensor (Euclidean M) (Euclidean M)
schur (tl,tr,br) = inverseSchurComplement (toTensor tl) (transpose tr) br

blockInversion ::
    ( Cartesian # Diagonal (Euclidean M) (Euclidean M)
    , Cartesian # Tensor (Euclidean M) (Euclidean N)
    , Cartesian # Symmetric (Euclidean N) (Euclidean N) )
    -> S.Matrix (N+M) (N+M) Double
blockInversion (tl,tr,br) =
    let (tl',tr',br') = blockSymmetricMatrixInversion tl tr br
        top = S.horizontalConcat (toMatrix $ toTensor tl') (toMatrix tr')
        btm = S.horizontalConcat (S.transpose $ toMatrix tr') (toMatrix $ toTensor br')
     in S.verticalConcat top btm



-- Benchmark

main :: IO ()
main = do

    g <- R.createSystemRandom

    cs1 <- S.replicateM $ R.uniformRM (-1,1) g
    cs2 <- S.replicateM $ R.uniformRM (-1,1) g
    cs3 <- S.replicateM $ R.uniformRM (-1,1) g
    cs4 <- S.replicateM $ R.uniformRM (-1,1) g

    let tl = Point cs1
        br = Point cs2
        tr = Point cs3
    let shr1 = schur (tl,tr,br)
        shr2 = woodbury (tl,tr,br)

    putStrLn "Woodbury error:"
    print $ euclideanDistance shr1 shr2

    let bigMatrix,bigMatrix2 :: S.Matrix (M+N) (M+N) Double
        bigMatrix =
            let top = S.horizontalConcat (toMatrix $ toTensor tl) (toMatrix tr)
                btm = S.horizontalConcat (S.transpose $ toMatrix tr) (toMatrix $ toTensor br)
             in S.verticalConcat top btm
        bigMatrix2 = G.Matrix cs4

    putStrLn "Inversion error:"
    let inv1 = blockInversion (tl,tr,br)
        inv2 = S.inverse bigMatrix

    print . sqrt . S.sum . S.map square $ G.toVector (inv1 - inv2)

    C.defaultMain
       [ C.bench "schur" $ C.nf schur (tl,tr,br)
       , C.bench "woodbury" $ C.nf woodbury (tl,tr,br)
       , C.bench "direct-inversion" $ C.nf S.inverse bigMatrix
       , C.bench "indirect-inversion" $ C.nf blockInversion (tl,tr,br)
       , C.bench "dense-inversion" $ C.nf S.inverse bigMatrix2 ]


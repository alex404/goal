{-# LANGUAGE
    TypeOperators,
    NoStarIsType,
    FlexibleContexts,
    DataKinds
    #-}

import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import qualified Criterion.Main as C
import qualified System.Random.MWC as R



--- Globals ---


type M = 1000
type N = 10


-- Functions


woodbury ::
    ( Cartesian # Diagonal (Euclidean M) (Euclidean M)
    , Cartesian # Tensor (Euclidean M) (Euclidean N)
    , Cartesian # Symmetric (Euclidean N) (Euclidean N) )
    -> Cartesian # Tensor (Euclidean M) (Euclidean M)
woodbury (diag,tns,trng) =
    let shr = schurComplement (toTensor trng) tns diag
     in woodburyMatrix diag tns shr

schur ::
    ( Cartesian # Diagonal (Euclidean M) (Euclidean M)
    , Cartesian # Tensor (Euclidean M) (Euclidean N)
    , Cartesian # Symmetric (Euclidean N) (Euclidean N) )
    -> Cartesian # Tensor (Euclidean M) (Euclidean M)
schur (diag,tns,trng) = schurComplement (toTensor diag) (transpose tns) trng

-- Benchmark

main :: IO ()
main = do

    g <- R.createSystemRandom

    cs1 <- S.replicateM $ R.uniformRM (-1,1) g
    cs2 <- S.replicateM $ R.uniformRM (-1,1) g
    cs3 <- S.replicateM $ R.uniformRM (-1,1) g

    let diag = Point cs1
        trng = Point cs2
        tns = Point cs3
    let shr1 = schur (diag,tns,trng)
        shr2 = woodbury (diag,tns,trng)

    putStrLn "Woodbury error:"
    print $ euclideanDistance shr1 shr2



    C.defaultMain
       [ C.bench "schur" $ C.nf schur (diag,tns,trng)
       , C.bench "woodbury" $ C.nf woodbury (diag,tns,trng) ]


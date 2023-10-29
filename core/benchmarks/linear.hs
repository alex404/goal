{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear

--- Misc

import Control.Monad (replicateM)
import Criterion.Main qualified as C
import Criterion.Types (Config (..))
import System.Random.MWC qualified as R

--- Globals ---

type N = 1000

n :: Int
n = 2000

--- Helper functions ---

conversion :: (LinearConstruct t n n) => Linear t n n -> Linear t n n
conversion = fromMatrix . toMatrix

invlndet :: (LinearConstruct t n n) => Linear t n n -> (S.Vector (LinearParameters t n n) Double, Double, Double)
invlndet f =
    let (inv, lndet, sgn) = inverseLogDeterminant f
     in (toVector inv, lndet, sgn)

vectorMultiplyTest :: (LinearConstruct t n n) => S.Vector n Double -> Linear t n n -> S.Vector n Double
vectorMultiplyTest v f = linearVectorMultiply f v

mapTest :: (LinearConstruct t n n) => [S.Vector n Double] -> Linear t n n -> [S.Vector n Double]
mapTest vs f = linearMap f vs

--- Benchmarks ---

-- Record over different test operators
data LinearTests = LinearTests
    { full :: Linear Full N N
    , symmetric :: Linear Symmetric N N
    , positiveDefinite :: Linear PositiveDefinite N N
    , diagonal :: Linear Diagonal N N
    , scale :: Linear Scale N N
    }

choleskyBenchmark :: Linear PositiveDefinite N N -> [C.Benchmark]
choleskyBenchmark f =
    [C.bench "Cholesky" $ C.nf (toVector . choleskyDecomposition) f]

conversions :: LinearTests -> [C.Benchmark]
conversions ls =
    [ C.bench "Full" $ C.nf (toVector . conversion) (full ls)
    , C.bench "Symmetric" $ C.nf (toVector . conversion) (symmetric ls)
    , C.bench "Positive Definite" $ C.nf (toVector . conversion) (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf (toVector . conversion) (diagonal ls)
    , C.bench "Scale" $ C.nf (toVector . conversion) (scale ls)
    ]

inversions :: LinearTests -> [C.Benchmark]
inversions ls =
    [ C.bench "Full" $ C.nf (toVector . inverse) (full ls)
    , C.bench "Symmetric" $ C.nf (toVector . inverse) (symmetric ls)
    , C.bench "Positive Definite" $ C.nf (toVector . inverse) (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf (toVector . inverse) (diagonal ls)
    , C.bench "Scale" $ C.nf (toVector . inverse) (scale ls)
    ]

invlndets :: LinearTests -> [C.Benchmark]
invlndets ls =
    [ C.bench "Full" $ C.nf invlndet (full ls)
    , C.bench "Symmetric" $ C.nf invlndet (symmetric ls)
    , C.bench "Positive Definite" $ C.nf invlndet (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf invlndet (diagonal ls)
    , C.bench "Scale" $ C.nf invlndet (scale ls)
    ]

determinants :: LinearTests -> [C.Benchmark]
determinants ls =
    [ C.bench "Full" $ C.nf determinant (full ls)
    , C.bench "Symmetric" $ C.nf determinant (symmetric ls)
    , C.bench "Positive Definite" $ C.nf determinant (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf determinant (diagonal ls)
    , C.bench "Scale" $ C.nf determinant (scale ls)
    ]

vectorMultiplies :: LinearTests -> S.Vector N Double -> [C.Benchmark]
vectorMultiplies ls v =
    [ C.bench "Full" $ C.nf (vectorMultiplyTest v) (full ls)
    , C.bench "Symmetric" $ C.nf (vectorMultiplyTest v) (symmetric ls)
    , C.bench "Positive Definite" $ C.nf (vectorMultiplyTest v) (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf (vectorMultiplyTest v) (diagonal ls)
    , C.bench "Scale" $ C.nf (vectorMultiplyTest v) (scale ls)
    ]

vectorMaps :: LinearTests -> [S.Vector N Double] -> [C.Benchmark]
vectorMaps ls vs =
    [ C.bench "Full" $ C.nf (mapTest vs) (full ls)
    , C.bench "Symmetric" $ C.nf (mapTest vs) (symmetric ls)
    , C.bench "Positive Definite" $ C.nf (mapTest vs) (positiveDefinite ls)
    , C.bench "Diagonal" $ C.nf (mapTest vs) (diagonal ls)
    , C.bench "Scale" $ C.nf (mapTest vs) (scale ls)
    ]

compositions :: LinearTests -> [C.Benchmark]
compositions ls =
    let fl = full ls
     in [ C.bench "Full" $ C.nf (toVector . linearCompose fl) (full ls)
        , C.bench "Symmetric" $ C.nf (toVector . linearCompose fl) (symmetric ls)
        , C.bench "PositiveDefinite" $ C.nf (toVector . linearCompose fl) (positiveDefinite ls)
        , C.bench "Diagonal" $ C.nf (toVector . linearCompose fl) (diagonal ls)
        , C.bench "Scale" $ C.nf (toVector . linearCompose fl) (scale ls)
        ]

--- Main ---

main :: IO ()
main = do

    g <- R.createSystemRandom

    vs :: [S.Vector N Double] <-
        replicateM n . S.replicateM $ R.uniformRM (0, 1) g

    let fl :: Linear Full N N
        fl = averageOuterProduct vs vs

    let sm :: Linear Symmetric N N
        sm = averageOuterProduct vs vs

    let pd :: Linear PositiveDefinite N N
        pd = averageOuterProduct vs vs

    let dg :: Linear Diagonal N N
        dg = averageOuterProduct vs vs

    let scl :: Linear Scale N N
        scl = averageOuterProduct vs vs

    let ls = LinearTests fl sm pd dg scl

    print "Positive Definite Test:"
    print . S.isSemiPositiveDefinite $ toMatrix pd

    bnchfl <- benchFilePath "linear.html"
    C.defaultMainWith
        C.defaultConfig
            { reportFile = Just bnchfl
            }
        [ C.bgroup "Cholesky" (choleskyBenchmark pd)
        , C.bgroup "Conversion" (conversions ls)
        , C.bgroup "Inversion" (inversions ls)
        , C.bgroup "Determinant" (determinants ls)
        , C.bgroup "Inverse Log-Determinant" (invlndets ls)
        , C.bgroup "Vector Multiplication" (vectorMultiplies ls $ head vs)
        , C.bgroup "Vector Maps" (vectorMaps ls vs)
        , C.bgroup "Linear Compositions" (compositions ls)
        ]

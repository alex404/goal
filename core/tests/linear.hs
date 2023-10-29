{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}

--- Imports ---

--- Goal

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Maybe (fromJust)
import Test.QuickCheck

--- Globals ---

numericalTolerance :: Double
numericalTolerance = 1e-12

--- Helpers Functions

approxEqual :: Double -> Double -> Bool
approxEqual x y = abs (x - y) < numericalTolerance

matrixApproxEqual :: S.Matrix 3 3 Double -> S.Matrix 3 3 Double -> Bool
matrixApproxEqual (G.Matrix xs) (G.Matrix ys) =
    S.and $ S.zipWith approxEqual xs ys

--- QuickCheck ---

--- Instances

instance Arbitrary (L.Linear L.PositiveDefinite 3 3) where
    arbitrary = do
        xs <- vectorOf 3 (choose (-1, 1))
        ys <- vectorOf 3 (choose (-1, 1))
        zs <- vectorOf 3 (choose (-1, 1))
        let x = fromJust $ S.fromList xs :: S.Vector 3 Double
        let y = fromJust $ S.fromList ys :: S.Vector 3 Double
        let z = fromJust $ S.fromList zs :: S.Vector 3 Double
        let pdmtx = x `S.outerProduct` x + y `S.outerProduct` y + z `S.outerProduct` z
        return . L.PositiveDefiniteLinear $ S.lowerTriangular pdmtx

--- Properties

choleskyInversion :: L.Linear L.PositiveDefinite 3 3 -> Property
choleskyInversion pdm@(L.PositiveDefiniteLinear trng) =
    let pdm' = L.SymmetricLinear trng
     in S.isPositiveDefinite (L.toMatrix pdm)
            ==> matrixApproxEqual (L.toMatrix $ L.inverse pdm) (L.toMatrix $ L.inverse pdm')

choleskyDeterminant :: L.Linear L.PositiveDefinite 3 3 -> Property
choleskyDeterminant pdm@(L.PositiveDefiniteLinear trng) =
    let pdm' = L.SymmetricLinear trng :: L.Linear L.Symmetric 3 3
     in S.isPositiveDefinite (L.toMatrix pdm)
            ==> approxEqual (L.determinant pdm) (L.determinant pdm')

--- Main ---

main :: IO ()
main = do
    -- Run tests
    putStrLn "Running tests..."
    putStrLn "\nCholesky Inversion...\n"
    result1 <- verboseCheckResult choleskyInversion
    putStrLn $ "\nCholesky Inversion Successful: " ++ show (isSuccess result1) ++ "\n\n"
    putStrLn "C\nholesky Determinant...\n"
    result2 <- verboseCheckResult choleskyDeterminant
    putStrLn $ "Cholesky Determinant Successful: " ++ show (isSuccess result2)

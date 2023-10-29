{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}

--- Imports ---

--- Goal

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Control.Monad (unless)
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
        xs <- vector 9
        let mtx0 = G.Matrix . fromJust $ S.fromList xs :: S.Matrix 3 3 Double
            mtx = mtx0 + S.matrixIdentity
            pdmtx = S.matrixMatrixMultiply mtx $ S.matrixMatrixMultiply (2 * S.matrixIdentity) (S.inverse mtx)
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
    putStrLn "Running tests..."
    result1 <- quickCheckResult choleskyInversion
    putStrLn $ "Cholesky inversion success: " ++ show (isSuccess result1)
    result2 <- quickCheckResult choleskyDeterminant
    putStrLn $ "Cholesky determinant success: " ++ show (isSuccess result1)

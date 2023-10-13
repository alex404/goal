{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Vectors and Matrices with statically typed dimensions based on storable vectors and using HMatrix where possible.
module Goal.Core.Vector.Storable.Linear where

--- Imports ---


-- Goal --


import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G
import GHC.TypeLits
import Goal.Core.Util (average,Triangular)


--- Types ---


data LinearType = Full | Symmetric | Diagonal | Scale | Identity

data Linear t n m where
    FullLinear :: S.Vector (m*n) Double -> Linear Full m n
    SymmetricLinear :: S.Vector (Triangular n) Double -> Linear Symmetric n n
    DiagonalLinear :: S.Vector n Double -> Linear Diagonal n n
    ScaleLinear :: Double -> Linear Scale n n
    IdentityLinear :: Linear Identity n n

deriving instance Show (Linear t m n)  -- for debugging


--- Construction and Manipulation ---


toMatrix :: (KnownNat m, KnownNat n) => Linear t m n -> S.Matrix m n Double
{-# INLINE toMatrix #-}
toMatrix (FullLinear m) = G.Matrix m
toMatrix (SymmetricLinear m) = S.triangularToSymmetric m
toMatrix (DiagonalLinear m) = S.diagonalMatrix m
toMatrix (ScaleLinear s) = S.diagonalMatrix (S.replicate s)
toMatrix IdentityLinear = S.matrixIdentity

type family LinearSize (t :: LinearType) (n :: Nat) (m :: Nat) :: Nat where
    LinearSize Full n m = n * m 
    LinearSize Symmetric n n = Triangular n
    LinearSize Diagonal n n = n
    LinearSize Scale n n = 1
    LinearSize Identity n n = 0

toVector :: (KnownNat m, KnownNat n) => Linear t m n -> S.Vector (LinearSize t m n) Double
{-# INLINE toVector #-}
toVector (FullLinear m) = m
toVector (SymmetricLinear m) = m
toVector (DiagonalLinear m) = m
toVector (ScaleLinear s) = S.replicate s
toVector IdentityLinear = S.empty

transpose :: (KnownNat m, KnownNat n) => Linear t m n -> Linear t n m
{-# INLINE transpose #-}
transpose m@(FullLinear _) = FullLinear . G.toVector . S.transpose $ toMatrix m
transpose m@(SymmetricLinear _) = m
transpose m@(DiagonalLinear _) = m
transpose m@(ScaleLinear _) = m
transpose IdentityLinear = IdentityLinear

class OuterProductable t m n where
    outerProduct :: (KnownNat m, KnownNat n) => S.Vector m Double -> S.Vector n Double -> Linear t m n
    averageOuterProduct :: (KnownNat m, KnownNat n) => [S.Vector m Double] -> [S.Vector n Double] -> Linear t m n

-- Full matrix instance: standard outer product
instance OuterProductable Full n m where
    outerProduct v1 v2 = FullLinear . G.toVector $ S.outerProduct v1 v2
    averageOuterProduct v1s v2s = FullLinear . G.toVector . S.averageOuterProduct $ zip v1s v2s

instance OuterProductable Symmetric n n where
    outerProduct v1 v2 = SymmetricLinear . S.lowerTriangular $ S.outerProduct v1 v2
    averageOuterProduct v1s v2s = SymmetricLinear . S.lowerTriangular . S.averageOuterProduct $ zip v1s v2s

-- Diagonal instance: elementwise product (should have same dimension)
instance OuterProductable Diagonal n n where
    outerProduct v1 v2 = DiagonalLinear $ v1 * v2
    averageOuterProduct v1s v2s = DiagonalLinear . average $ zipWith (*) v1s v2s

-- Diagonal instance: elementwise product (should have same dimension)
instance OuterProductable Scale n n where
    outerProduct v1 v2 = ScaleLinear . S.average $ v1 * v2
    averageOuterProduct v1s v2s = ScaleLinear . average $ zipWith (\v1 v2 -> S.average $ v1 * v2) v1s v2s

-- Diagonal instance: elementwise product (should have same dimension)
instance OuterProductable Identity n n where
    outerProduct _ _ = IdentityLinear
    averageOuterProduct _ _ = IdentityLinear

-- Other instances can be added similarly

--- Operatiors ---


-- | Inversion for general linear operators.
inverse :: KnownNat n => Linear t n n -> Linear t n n
{-# INLINE inverse #-}
inverse m@(FullLinear _) = FullLinear . G.toVector . S.inverse $ toMatrix m
inverse m@(SymmetricLinear _) = SymmetricLinear . S.lowerTriangular . S.inverse $ toMatrix m
inverse (DiagonalLinear m) = DiagonalLinear $ recip m
inverse (ScaleLinear s) = ScaleLinear (recip s)
inverse IdentityLinear = IdentityLinear

-- | Cholesky decomposition for symmetric linear operators (note: does not check positive
-- definiteness).
choleskyDecomposition :: KnownNat n => Linear Symmetric n n -> Linear Full n n
{-# INLINE choleskyDecomposition #-}
choleskyDecomposition m = FullLinear . G.toVector . S.transpose . S.unsafeCholesky $ toMatrix m

-- | Matrix inversion based on Cholesky decomposition for symmetric linear operators
-- (note: does not check positive definiteness).
choleskyInversion :: KnownNat n => Linear Symmetric n n -> Linear Symmetric n n
{-# INLINE choleskyInversion #-}
choleskyInversion m = SymmetricLinear . S.lowerTriangular . S.unsafeCholeskyInversion $ toMatrix m


--- Vector multiplication ---


matrixVectorMultiply :: (KnownNat n, KnownNat m) => Linear t m n -> S.Vector n Double -> S.Vector m Double
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply (FullLinear m) v = S.matrixVectorMultiply (G.Matrix m) v
matrixVectorMultiply (SymmetricLinear m) v = S.matrixVectorMultiply (S.triangularToSymmetric m) v
matrixVectorMultiply (DiagonalLinear m) v = m * v
matrixVectorMultiply (ScaleLinear s) v = S.scale s v
matrixVectorMultiply IdentityLinear v = v

matrixMap :: (KnownNat n, KnownNat m) => Linear t m n -> [S.Vector n Double] -> [S.Vector m Double]
{-# INLINE matrixMap #-}
matrixMap (FullLinear m) vs = S.matrixMap (G.Matrix m) vs
matrixMap (SymmetricLinear m) vs = S.matrixMap (S.triangularToSymmetric m) vs
matrixMap (DiagonalLinear m) vs = S.diagonalMatrixMap m vs
matrixMap (ScaleLinear s) vs = map (S.scale s) vs
matrixMap IdentityLinear vs = vs


--- Matrix multiplication ---


type family LinearMultiply (s :: LinearType) (t :: LinearType) :: LinearType where
    LinearMultiply Identity t = t
    LinearMultiply s Identity = s
    LinearMultiply Scale t = t
    LinearMultiply s Scale = s
    LinearMultiply Diagonal Diagonal = Diagonal
    LinearMultiply _ _ = Full


matrixMatrixMultiply 
    :: forall n m o s t
    . (KnownNat n, KnownNat m, KnownNat o) 
    => Linear s m n -> Linear t n o -> Linear (LinearMultiply s t) m o
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply m IdentityLinear = m
matrixMatrixMultiply IdentityLinear m = m
matrixMatrixMultiply (ScaleLinear s) (SymmetricLinear m) = SymmetricLinear $ S.scale s m
matrixMatrixMultiply (SymmetricLinear m) (ScaleLinear s) = SymmetricLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s) (DiagonalLinear m) = DiagonalLinear $ S.scale s m
matrixMatrixMultiply (DiagonalLinear m) (ScaleLinear s) = DiagonalLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s1) (ScaleLinear s2) = ScaleLinear $ s1 * s2
matrixMatrixMultiply (FullLinear m) (ScaleLinear s) = FullLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s) (FullLinear m) = FullLinear $ S.scale s m
matrixMatrixMultiply (DiagonalLinear m1) (DiagonalLinear m2) = DiagonalLinear $ m1 * m2
matrixMatrixMultiply (DiagonalLinear m1) m2@(FullLinear _) =
   FullLinear . G.toVector . S.diagonalMatrixMatrixMultiply m1 $ toMatrix m2
matrixMatrixMultiply  m2@(FullLinear _) (DiagonalLinear m1) =
   FullLinear . G.toVector . S.transpose . S.diagonalMatrixMatrixMultiply m1 . S.transpose $ toMatrix m2
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(FullLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(FullLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(FullLinear _) m2@(FullLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(DiagonalLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(DiagonalLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2




--- Internal ---


fullMultiply 
    :: forall n m o s t
    . (KnownNat n, KnownNat m, KnownNat o) 
    => Linear s m n -> Linear t n o -> Linear Full m o
fullMultiply m' m'' = 
  let mtx1 :: S.Matrix m n Double
      mtx1 = toMatrix m'
      mtx2 :: S.Matrix n o Double
      mtx2 = toMatrix m''
  in FullLinear . G.toVector $ S.matrixMatrixMultiply mtx1 mtx2

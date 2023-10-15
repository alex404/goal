{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | A generalized representation of a linear operator with efficient implementations for simpler
cases.
-}
module Goal.Core.Vector.Storable.Linear (
    -- * Types
    LinearRep (..),
    Linear (..),
    LinearParameters,

    -- * Construction and Manipulation
    toVector,
    toMatrix,
    transpose,
    SquareConstruct (..),
    LinearConstruct (..),

    -- * Operations
    inverse,
    choleskyDecomposition,
    determinant,
    inverseLogDeterminant,

    -- * Vector multiplication
    matrixVectorMultiply,
    matrixMap,

    -- * Matrix multiplication
    LinearMultiply,
    matrixMatrixMultiply,
) where

--- Imports ---

-- Goal --

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S

import Data.Proxy (Proxy (Proxy))
import GHC.TypeNats (KnownNat, Nat, natVal, type (*))

import Goal.Core.Util (Triangular, average, natValInt, square)

--- Types ---

-- | The internal representation of a linear operator.
data LinearRep = Full | Symmetric | PositiveDefinite | Diagonal | Scale | Identity

-- | A linear operator.
data Linear t n m where
    -- | A full matrix.
    FullLinear :: S.Vector (m * n) Double -> Linear Full m n
    -- | A symmetric matrix.
    SymmetricLinear :: S.Vector (Triangular n) Double -> Linear Symmetric n n
    -- | A positive definite matrix (note: this is not checked).
    PositiveDefiniteLinear :: S.Vector (Triangular n) Double -> Linear PositiveDefinite n n
    -- | A diagonal matrix.
    DiagonalLinear :: S.Vector n Double -> Linear Diagonal n n
    -- | A scalar matrix.
    ScaleLinear :: Double -> Linear Scale n n
    -- | The identity matrix.
    IdentityLinear :: Linear Identity n n

deriving instance Show (Linear t m n) -- for debugging

--- Construction and Manipulation ---

-- | Type-level representation of the number of parameters in a linear operator.
type family LinearParameters (t :: LinearRep) (n :: Nat) (m :: Nat) :: Nat where
    LinearParameters Full n m = n * m
    LinearParameters Symmetric n n = Triangular n
    LinearParameters PositiveDefinite n n = Triangular n
    LinearParameters Diagonal n n = n
    LinearParameters Scale n n = 1
    LinearParameters Identity n n = 0

-- | A vector of the parameters in a linear operator.
toVector :: (KnownNat m, KnownNat n) => Linear t m n -> S.Vector (LinearParameters t m n) Double
{-# INLINE toVector #-}
toVector (FullLinear m) = m
toVector (SymmetricLinear m) = m
toVector (PositiveDefiniteLinear m) = m
toVector (DiagonalLinear m) = m
toVector (ScaleLinear s) = S.replicate s
toVector IdentityLinear = S.empty

-- | Convert a linear operator to a storable 'Matrix'.
toMatrix :: (KnownNat m, KnownNat n) => Linear t m n -> S.Matrix m n Double
{-# INLINE toMatrix #-}
toMatrix (FullLinear m) = G.Matrix m
toMatrix (SymmetricLinear m) = S.triangularToSymmetric m
toMatrix (PositiveDefiniteLinear m) = S.triangularToSymmetric m
toMatrix (DiagonalLinear m) = S.diagonalMatrix m
toMatrix (ScaleLinear s) = S.diagonalMatrix (S.replicate s)
toMatrix IdentityLinear = S.matrixIdentity

-- | Transpose of the linear operator.
transpose :: (KnownNat m, KnownNat n) => Linear t m n -> Linear t n m
{-# INLINE transpose #-}
transpose m@(FullLinear _) = FullLinear . G.toVector . S.transpose $ toMatrix m
transpose m@(SymmetricLinear _) = m
transpose m@(PositiveDefiniteLinear _) = m
transpose m@(DiagonalLinear _) = m
transpose m@(ScaleLinear _) = m
transpose IdentityLinear = IdentityLinear

-- | Construction of linear operators from the outer product.
class SquareConstruct t n where
    identity :: (KnownNat n) => Linear t n n

instance SquareConstruct Full n where
    identity =
        let idnt :: S.Matrix n n Double
            idnt = S.matrixIdentity
         in FullLinear $ G.toVector idnt

instance SquareConstruct Symmetric n where
    identity =
        let idnt :: S.Matrix n n Double
            idnt = S.matrixIdentity
         in SymmetricLinear $ S.lowerTriangular idnt

instance SquareConstruct PositiveDefinite n where
    identity =
        let idnt :: S.Matrix n n Double
            idnt = S.matrixIdentity
         in PositiveDefiniteLinear $ S.lowerTriangular idnt

instance SquareConstruct Diagonal n where
    identity = DiagonalLinear $ S.replicate 1

instance SquareConstruct Scale n where
    identity = ScaleLinear 1

instance SquareConstruct Identity n where
    identity = IdentityLinear

-- | Construction of linear operators from the outer product.
class (SquareConstruct t m, KnownNat m, KnownNat n) => LinearConstruct t m n where
    outerProduct :: (KnownNat m, KnownNat n) => S.Vector m Double -> S.Vector n Double -> Linear t m n
    averageOuterProduct :: (KnownNat m, KnownNat n) => [S.Vector m Double] -> [S.Vector n Double] -> Linear t m n
    fromMatrix :: (KnownNat m, KnownNat n) => S.Matrix m n Double -> Linear t m n

-- Full matrix instance: standard outer product
instance (KnownNat n, KnownNat m) => LinearConstruct Full n m where
    outerProduct v1 v2 = FullLinear . G.toVector $ S.outerProduct v1 v2
    averageOuterProduct v1s v2s = FullLinear . G.toVector . S.averageOuterProduct $ zip v1s v2s
    fromMatrix = FullLinear . G.toVector

instance (KnownNat n) => LinearConstruct Symmetric n n where
    outerProduct v1 v2 = SymmetricLinear . S.lowerTriangular $ S.outerProduct v1 v2
    averageOuterProduct v1s v2s = SymmetricLinear . S.lowerTriangular . S.averageOuterProduct $ zip v1s v2s
    fromMatrix = SymmetricLinear . S.lowerTriangular

instance (KnownNat n) => LinearConstruct PositiveDefinite n n where
    outerProduct v1 v2 = PositiveDefiniteLinear . S.lowerTriangular $ S.outerProduct v1 v2
    averageOuterProduct v1s v2s = PositiveDefiniteLinear . S.lowerTriangular . S.averageOuterProduct $ zip v1s v2s
    fromMatrix = PositiveDefiniteLinear . S.lowerTriangular

-- Diagonal instance: elementwise product (should have same dimension)
instance (KnownNat n) => LinearConstruct Diagonal n n where
    outerProduct v1 v2 = DiagonalLinear $ v1 * v2
    averageOuterProduct v1s v2s = DiagonalLinear . average $ zipWith (*) v1s v2s
    fromMatrix = DiagonalLinear . S.takeDiagonal

-- Diagonal instance: elementwise product (should have same dimension)
instance (KnownNat n) => LinearConstruct Scale n n where
    outerProduct v1 v2 = ScaleLinear . S.average $ v1 * v2
    averageOuterProduct v1s v2s = ScaleLinear . average $ zipWith (\v1 v2 -> S.average $ v1 * v2) v1s v2s
    fromMatrix = ScaleLinear . S.average . S.takeDiagonal

-- Diagonal instance: elementwise product (should have same dimension)
instance (KnownNat n) => LinearConstruct Identity n n where
    outerProduct _ _ = IdentityLinear
    averageOuterProduct _ _ = IdentityLinear
    fromMatrix _ = IdentityLinear

-- Other instances can be added similarly

--- Operators ---

{- | Cholesky decomposition for symmetric linear operators (note: does not check positive
definiteness).
-}
choleskyDecomposition :: (KnownNat n) => Linear PositiveDefinite n n -> Linear Full n n
{-# INLINE choleskyDecomposition #-}
choleskyDecomposition m = FullLinear . G.toVector . S.transpose . S.unsafeCholesky $ toMatrix m

{- | Matrix inversion based on Cholesky decomposition for symmetric linear operators
(note: does not check positive definiteness).
-}
choleskyInversion :: (KnownNat n) => Linear PositiveDefinite n n -> Linear PositiveDefinite n n
{-# INLINE choleskyInversion #-}
choleskyInversion m = PositiveDefiniteLinear . S.lowerTriangular . S.unsafeCholeskyInversion $ toMatrix m

choleskyDeterminant :: (KnownNat n) => Linear PositiveDefinite n n -> Double
{-# INLINE choleskyDeterminant #-}
choleskyDeterminant = S.product . square . S.takeDiagonal . toMatrix . choleskyDecomposition

logCholeskyDeterminant :: (KnownNat n) => Linear PositiveDefinite n n -> Double
{-# INLINE logCholeskyDeterminant #-}
logCholeskyDeterminant = (* 2) . S.sum . log . S.takeDiagonal . toMatrix . choleskyDecomposition

-- | Inversion for general linear operators.
inverse :: (KnownNat n) => Linear t n n -> Linear t n n
{-# INLINE inverse #-}
inverse m@(FullLinear _) = FullLinear . G.toVector . S.inverse $ toMatrix m
inverse m@(SymmetricLinear _) = SymmetricLinear . S.lowerTriangular . S.inverse $ toMatrix m
inverse m@(PositiveDefiniteLinear _) = choleskyInversion m
inverse (DiagonalLinear m) = DiagonalLinear $ recip m
inverse (ScaleLinear s) = ScaleLinear (recip s)
inverse IdentityLinear = IdentityLinear

determinant :: forall t n. (KnownNat n) => Linear t n n -> Double
{-# INLINE determinant #-}
determinant m@(FullLinear _) = S.determinant $ toMatrix m
determinant m@(SymmetricLinear _) = S.determinant $ toMatrix m
determinant m@(PositiveDefiniteLinear _) = choleskyDeterminant m
determinant (DiagonalLinear m) = S.product m
determinant (ScaleLinear s) = s ^ natValInt (Proxy @n)
determinant IdentityLinear = 1

inverseLogDeterminant :: forall t n. (LinearConstruct t n n) => Linear t n n -> (Linear t n n, Double, Double)
{-# INLINE inverseLogDeterminant #-}
inverseLogDeterminant (ScaleLinear s) =
    let lndet = (fromIntegral . natVal $ Proxy @n) * log (abs s)
     in (inverse (ScaleLinear s), lndet, signum s)
inverseLogDeterminant (DiagonalLinear d) =
    let (lndet, sgn) =
            if S.all (> 0) $ signum d
                then (S.sum $ log d, 1)
                else
                    let prd = S.product d
                        sgn0 = signum prd
                     in (log $ abs prd, sgn0)
     in (inverse (DiagonalLinear d), lndet, sgn)
inverseLogDeterminant f@(PositiveDefiniteLinear _) =
    let chol = S.unsafeCholesky $ toMatrix f
        inv = S.unsafeCholeskyInversion0 chol
        lndet = (* 2) . S.sum . log $ S.takeDiagonal chol
     in (fromMatrix inv, lndet, 1)
inverseLogDeterminant f =
    let (imtx, lndet, sgn) = S.inverseLogDeterminant $ toMatrix f
     in (fromMatrix imtx, lndet, sgn)

--- Vector multiplication ---

-- | Multiply a vector by a linear operator.
matrixVectorMultiply :: (KnownNat n, KnownNat m) => Linear t m n -> S.Vector n Double -> S.Vector m Double
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply (FullLinear m) v = S.matrixVectorMultiply (G.Matrix m) v
matrixVectorMultiply (SymmetricLinear m) v = S.matrixVectorMultiply (S.triangularToSymmetric m) v
matrixVectorMultiply (PositiveDefiniteLinear m) v = S.matrixVectorMultiply (S.triangularToSymmetric m) v
matrixVectorMultiply (DiagonalLinear m) v = m * v
matrixVectorMultiply (ScaleLinear s) v = S.scale s v
matrixVectorMultiply IdentityLinear v = v

-- | Efficiently multiply a list of vectors by a linear operator.
matrixMap :: (KnownNat n, KnownNat m) => Linear t m n -> [S.Vector n Double] -> [S.Vector m Double]
{-# INLINE matrixMap #-}
matrixMap (FullLinear m) vs = S.matrixMap (G.Matrix m) vs
matrixMap (SymmetricLinear m) vs = S.matrixMap (S.triangularToSymmetric m) vs
matrixMap (PositiveDefiniteLinear m) vs = S.matrixMap (S.triangularToSymmetric m) vs
matrixMap (DiagonalLinear m) vs = S.diagonalMatrixMap m vs
matrixMap (ScaleLinear s) vs = map (S.scale s) vs
matrixMap IdentityLinear vs = vs

--- Matrix multiplication ---

-- | A representation of the simplest type yielded by a matrix-matrix multiplication.
type family LinearMultiply (s :: LinearRep) (t :: LinearRep) :: LinearRep where
    LinearMultiply Identity t = t
    LinearMultiply s Identity = s
    LinearMultiply Scale t = t
    LinearMultiply s Scale = s
    LinearMultiply Diagonal Diagonal = Diagonal
    LinearMultiply _ _ = Full

-- | Multiply two linear operators.
matrixMatrixMultiply ::
    forall n m o s t.
    (KnownNat n, KnownNat m, KnownNat o) =>
    Linear s m n ->
    Linear t n o ->
    Linear (LinearMultiply s t) m o
{-# INLINE matrixMatrixMultiply #-}
-- Identity
matrixMatrixMultiply m IdentityLinear = m
matrixMatrixMultiply IdentityLinear m = m
-- Scale
matrixMatrixMultiply (ScaleLinear s) (PositiveDefiniteLinear m) = PositiveDefiniteLinear $ S.scale s m
matrixMatrixMultiply (PositiveDefiniteLinear m) (ScaleLinear s) = PositiveDefiniteLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s) (SymmetricLinear m) = SymmetricLinear $ S.scale s m
matrixMatrixMultiply (SymmetricLinear m) (ScaleLinear s) = SymmetricLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s) (DiagonalLinear m) = DiagonalLinear $ S.scale s m
matrixMatrixMultiply (DiagonalLinear m) (ScaleLinear s) = DiagonalLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s1) (ScaleLinear s2) = ScaleLinear $ s1 * s2
matrixMatrixMultiply (FullLinear m) (ScaleLinear s) = FullLinear $ S.scale s m
matrixMatrixMultiply (ScaleLinear s) (FullLinear m) = FullLinear $ S.scale s m
-- Diagonal
matrixMatrixMultiply (DiagonalLinear d1) (DiagonalLinear d2) = DiagonalLinear $ d1 * d2
matrixMatrixMultiply d@(DiagonalLinear _) m@(FullLinear _) = diagonalMultiply d m
matrixMatrixMultiply m@(FullLinear _) d@(DiagonalLinear _) = transposeDiagonalMultiply d m
matrixMatrixMultiply m@(SymmetricLinear _) d@(DiagonalLinear _) = diagonalMultiply d m
matrixMatrixMultiply d@(DiagonalLinear _) m@(SymmetricLinear _) = transposeDiagonalMultiply d m
matrixMatrixMultiply m@(PositiveDefiniteLinear _) d@(DiagonalLinear _) = diagonalMultiply d m
matrixMatrixMultiply d@(DiagonalLinear _) m@(PositiveDefiniteLinear _) = transposeDiagonalMultiply d m
-- Symmetric
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(FullLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(FullLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(PositiveDefiniteLinear _) m2@(SymmetricLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(SymmetricLinear _) m2@(PositiveDefiniteLinear _) = fullMultiply m1 m2
-- PositiveDefinite
matrixMatrixMultiply m1@(PositiveDefiniteLinear _) m2@(PositiveDefiniteLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(PositiveDefiniteLinear _) m2@(FullLinear _) = fullMultiply m1 m2
matrixMatrixMultiply m1@(FullLinear _) m2@(PositiveDefiniteLinear _) = fullMultiply m1 m2
-- Full
matrixMatrixMultiply m1@(FullLinear _) m2@(FullLinear _) = fullMultiply m1 m2

--- Internal ---

fullMultiply ::
    forall n m o s t.
    (KnownNat n, KnownNat m, KnownNat o) =>
    Linear s m n ->
    Linear t n o ->
    Linear Full m o
fullMultiply m' m'' =
    let mtx1 :: S.Matrix m n Double
        mtx1 = toMatrix m'
        mtx2 :: S.Matrix n o Double
        mtx2 = toMatrix m''
     in FullLinear . G.toVector $ S.matrixMatrixMultiply mtx1 mtx2

diagonalMultiply ::
    forall n m t.
    (KnownNat n, KnownNat m) =>
    Linear Diagonal m m ->
    Linear t m n ->
    Linear Full m n
diagonalMultiply (DiagonalLinear d) m = FullLinear . G.toVector . S.diagonalMatrixMatrixMultiply d $ toMatrix m

transposeDiagonalMultiply ::
    forall n m t.
    (KnownNat n, KnownNat m) =>
    Linear Diagonal m m ->
    Linear t n m ->
    Linear Full n m
transposeDiagonalMultiply (DiagonalLinear d) m = FullLinear . G.toVector . S.transpose . S.diagonalMatrixMatrixMultiply d . S.transpose $ toMatrix m

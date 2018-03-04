{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLVM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector.Storable
    ( -- * Vector
      module Data.Vector.Storable.Sized
    -- * Matrix
    , Matrix (Matrix,toVector)
    -- ** Construction
    , fromRows
    , fromColumns
    , matrixIdentity
    , outerProduct
    -- ** Deconstruction
    , toRows
    , toColumns
    , nRows
    , nColumns
    -- ** Manipulation
    , columnVector
    , rowVector
    , diagonalConcat
    -- ** BLAS
    , dotProduct
    , matrixVectorMultiply
    , matrixMatrixMultiply
    , inverse
    , transpose
    ) where


--- Imports ---


import GHC.TypeLits
import Data.Proxy
import Control.DeepSeq
import Foreign.Storable
import Goal.Core.Vector.TypeLits
import Unsafe.Coerce
import Data.Vector.Storable.Sized
import Numeric.LinearAlgebra (Field,Numeric)

-- Qualified Imports --

import qualified Data.Vector.Storable as V
import qualified Numeric.LinearAlgebra as H

import Prelude hiding (concatMap,replicate,map,(++))


-- | Matrices with static dimensions.
newtype Matrix (m :: Nat) (n :: Nat) a = Matrix { toVector :: Vector (m*n) a }
    deriving (Eq,Show,NFData)

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector n a -> Matrix n 1 a
{-# INLINE columnVector #-}
columnVector = Matrix

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector n a -> Matrix 1 n a
{-# INLINE rowVector #-}
rowVector = Matrix

toHMatrix :: forall m n x . (KnownNat n, Storable x) => Matrix m n x -> H.Matrix x
toHMatrix (Matrix mtx) =
    H.reshape (natValInt (Proxy :: Proxy n)) $ fromSized mtx

fromHMatrix :: Numeric x => H.Matrix x -> Matrix m n x
fromHMatrix = unsafeCoerce . H.flatten

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: (KnownNat n, Storable x) => Vector m (Vector n x) -> Matrix m n x
{-# INLINE fromRows #-}
fromRows vs = Matrix $ concatMap id vs

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
fromColumns :: (KnownNat m, KnownNat n, Numeric x) => Vector n (Vector m x) -> Matrix m n x
{-# INLINE fromColumns #-}
fromColumns = transpose . fromRows

breakEvery :: forall n k a . (KnownNat n, KnownNat k, Storable a) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> unsafeCoerce $ V.unsafeSlice (i*k) k v)


--- BLAS ---


-- | The dot product of two numerical 'Vector's.
dotProduct :: Numeric x => Vector n x -> Vector n x -> x
{-# INLINE dotProduct #-}
dotProduct v1 v2 = H.dot (fromSized v1) (fromSized v2)

-- | Transpose a 'Matrix'.
transpose :: forall m n x . (KnownNat m, KnownNat n, Numeric x) => Matrix m n x -> Matrix n m x
{-# INLINE transpose #-}
transpose (Matrix mtx) =
    Matrix $ withVectorUnsafe (H.flatten . H.tr . H.reshape (natValInt (Proxy :: Proxy n))) mtx

-- | The number of rows in the 'Matrix'.
nRows :: forall m n a . KnownNat m => Matrix m n a -> Int
{-# INLINE nRows #-}
nRows _ = natValInt (Proxy :: Proxy m)

-- | The columns of rows in the 'Matrix'.
nColumns :: forall m n a . KnownNat n => Matrix m n a -> Int
{-# INLINE nColumns #-}
nColumns _ = natValInt (Proxy :: Proxy n)

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (KnownNat m, KnownNat n, Storable x) => Matrix m n x -> Vector m (Vector n x)
{-# INLINE toRows #-}
toRows (Matrix v) = breakEvery v

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns :: (KnownNat m, KnownNat n, Numeric x) => Matrix m n x -> Vector n (Vector m x)
{-# INLINE toColumns #-}
toColumns = toRows . transpose

-- | Diagonally concatenate two matrices, padding the gaps with zeroes.
diagonalConcat
    :: (KnownNat n, KnownNat m, KnownNat o, KnownNat p, Numeric x)
    => Matrix n m x -> Matrix o p x -> Matrix (n+o) (m+p) x
{-# INLINE diagonalConcat #-}
diagonalConcat mtx10 mtx20 =
    let mtx1 = toHMatrix mtx10
        mtx2 = toHMatrix mtx20
     in fromHMatrix $ H.diagBlock [mtx1,mtx2]

-- | Invert a 'Matrix'.
inverse :: forall n x . (KnownNat n, Field x) => Matrix n n x -> Matrix n n x
{-# INLINE inverse #-}
inverse (Matrix mtx) =
    Matrix $ withVectorUnsafe (H.flatten . H.inv . H.reshape (natValInt (Proxy :: Proxy n))) mtx

-- | The outer product of two 'Vector's.
outerProduct :: (KnownNat m, KnownNat n, Numeric x) => Vector m x -> Vector n x -> Matrix m n x
{-# INLINE outerProduct #-}
outerProduct v1 v2 =
    fromHMatrix $ H.outer (fromSized v1) (fromSized v2)

-- | The identity 'Matrix'.
matrixIdentity :: forall n x . (KnownNat n, Numeric x, Num x) => Matrix n n x
{-# INLINE matrixIdentity #-}
matrixIdentity =
    fromHMatrix . H.ident $ natValInt (Proxy :: Proxy n)

-- | Apply a linear transformation to a 'Vector'.
matrixVectorMultiply :: (KnownNat m, KnownNat n, Numeric x)
                     => Matrix m n x -> Vector n x -> Vector m x
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply mtx v =
    unsafeCoerce $ toHMatrix mtx H.#> fromSized v

-- | Multiply a 'Matrix' with a second 'Matrix'.
matrixMatrixMultiply
    :: (KnownNat m, KnownNat n, KnownNat o, Numeric x)
    => Matrix m n x
    -> Matrix n o x
    -> Matrix m o x
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply mtx1 mtx2 = fromHMatrix $ toHMatrix mtx1 H.<> toHMatrix mtx2

-- | Prety print the values of a 'Matrix' (for extremely simple values of pretty).
prettyPrintMatrix :: (KnownNat m, KnownNat n, Numeric a, Show a) => Matrix m n a -> IO ()
prettyPrintMatrix = print . toHMatrix

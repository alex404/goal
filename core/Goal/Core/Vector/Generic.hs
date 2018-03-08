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
module Goal.Core.Vector.Generic
    ( -- * Vector
      module Data.Vector.Generic.Sized
    , GVector
    , concat
    , doubleton
    , breakEvery
    -- * Matrix
    , Matrix (Matrix,toVector)
    -- ** Construction
    , fromRows
    , fromColumns
    -- ** Deconstruction
    , toRows
    , toColumns
    , nRows
    , nColumns
    -- ** Manipulation
    , columnVector
    , rowVector
    ) where


--- Imports ---


import GHC.TypeLits
import Data.Proxy
import Control.DeepSeq
import Goal.Core.Vector.TypeLits
import Data.Vector.Generic.Sized

import qualified Data.Vector.Generic as G

import Unsafe.Coerce

import Prelude hiding (concatMap,concat)


--- Vector ---

type GVector = G.Vector

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
concat :: (KnownNat n, GVector v x, GVector v (Vector v n x)) => Vector v m (Vector v n x) -> Vector v (m*n) x
{-# INLINE concat #-}
concat = concatMap id

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
doubleton :: GVector v x => x -> x -> Vector v 2 x
{-# INLINE doubleton #-}
doubleton x1 x2 = cons x1 $ singleton x2



-- | Matrices with static dimensions.
newtype Matrix v (m :: Nat) (n :: Nat) a = Matrix { toVector :: Vector v (m*n) a }
    deriving (Eq,Show,NFData)

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector v n a -> Matrix v n 1 a
{-# INLINE columnVector #-}
columnVector = Matrix

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector v n a -> Matrix v 1 n a
{-# INLINE rowVector #-}
rowVector = Matrix

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: (GVector v x, GVector v (Vector v n x), KnownNat n) => Vector v m (Vector v n x) -> Matrix v m n x
{-# INLINE fromRows #-}
fromRows = Matrix . concat

transpose
    :: forall v m n a . (KnownNat m, KnownNat n, GVector v Int, GVector v a)
    => Matrix v m n a -> Matrix v n m a
{-# INLINE transpose #-}
transpose (Matrix v0) =
    let v = fromSized v0
        m = natValInt (Proxy :: Proxy m)
        n = natValInt (Proxy :: Proxy n)
        vi = G.concatMap (\i -> G.generate m (\j -> i + j*n)) $ G.generate n id
     in Matrix . unsafeCoerce $ G.unsafeBackpermute v vi


-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
fromColumns
    :: (GVector v x, GVector v (Vector v m x), GVector v Int, KnownNat n, KnownNat m)
    => Vector v n (Vector v m x) -> Matrix v m n x
{-# INLINE fromColumns #-}
fromColumns = transpose . fromRows

breakEvery
    :: forall v n k a . (GVector v (Vector v k a), GVector v a, KnownNat n, KnownNat k)
    => Vector v (n*k) a -> Vector v n (Vector v k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> unsafeCoerce $ G.unsafeSlice (i*k) k v)


--- BLAS ---


-- | The number of rows in the 'Matrix'.
nRows :: forall v m n a . KnownNat m => Matrix v m n a -> Int
{-# INLINE nRows #-}
nRows _ = natValInt (Proxy :: Proxy m)

-- | The columns of rows in the 'Matrix'.
nColumns :: forall v m n a . KnownNat n => Matrix v m n a -> Int
{-# INLINE nColumns #-}
nColumns _ = natValInt (Proxy :: Proxy n)

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (GVector v (Vector v n a), GVector v a, KnownNat n, KnownNat m)
       => Matrix v m n a -> Vector v m (Vector v n a)
{-# INLINE toRows #-}
toRows (Matrix v) = breakEvery v

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns
    :: (GVector v (Vector v m a), GVector v a, KnownNat m, KnownNat n, GVector v Int)
    => Matrix v m n a -> Vector v n (Vector v m a)
{-# INLINE toColumns #-}
toColumns = toRows . transpose

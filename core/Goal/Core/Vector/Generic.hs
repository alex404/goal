{-# LANGUAGE StandaloneDeriving,GeneralizedNewtypeDeriving #-}
 {-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
-- | Vectors and Matrices with statically typed dimensions.
module Goal.Core.Vector.Generic
    ( -- * Vector
      module Data.Vector.Generic.Sized
    , VectorClass
    -- * Construction
    , doubleton
    , range
    , breakEvery
    -- * Deconstruction
    , concat
    -- * Matrix
    , Matrix (Matrix,toVector)
    , nRows
    , nColumns
    -- ** Construction
    , fromRows
    , fromColumns
    -- ** Deconstruction
    , toPair
    , toRows
    , toColumns
    -- ** Manipulation
    , columnVector
    , rowVector
    -- ** BLAS
    , transpose
    , dotProduct
    , weakDotProduct
    , outerProduct
    , matrixVectorMultiply
    , matrixMatrixMultiply
    ) where


--- Imports ---


-- Goal --

import Goal.Core.Util hiding (breakEvery,range)

-- Unqualified --

import GHC.TypeNats
import Data.Proxy
import Control.DeepSeq
import Data.Vector.Generic.Sized
import Data.Vector.Generic.Sized.Internal
import Foreign.Storable
import Prelude hiding (concatMap,concat,map,sum,replicate)

-- Qualified --

import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S

import Numeric.LinearAlgebra (Numeric)

--- Vector ---


type VectorClass = G.Vector

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
concat :: (KnownNat n, G.Vector v x, G.Vector v (Vector v n x)) => Vector v m (Vector v n x) -> Vector v (m*n) x
{-# INLINE concat #-}
concat = concatMap id

-- | Collect two values into a length 2 'Vector'.
doubleton :: G.Vector v x => x -> x -> Vector v 2 x
{-# INLINE doubleton #-}
doubleton x1 x2 = cons x1 $ singleton x2

-- | Breaks a 'Vector' into a Vector of Vectors.
breakEvery
    :: forall v n k a . (G.Vector v a, G.Vector v (Vector v k a), KnownNat n, KnownNat k)
    => Vector v (n*k) a -> Vector v n (Vector v k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> Vector $ G.unsafeSlice (finiteInt i*k) k v)

-- | Reshapes a length 2 'Vector' into a pair of values.
toPair :: G.Vector v a => Vector v 2 a -> (a,a)
{-# INLINE toPair #-}
toPair v = (unsafeIndex v 0, unsafeIndex v 1)

-- | Uniform partition of an interval into a 'Vector'.
range
    :: forall v n x. (G.Vector v x, KnownNat n, Fractional x)
    => x -> x -> Vector v n x
{-# INLINE range #-}
range mn mx =
    let n = natValInt (Proxy :: Proxy n)
        stp = (mx - mn)/fromIntegral (n-1)
     in enumFromStepN mn stp


--- Matrix ---


-- | Matrices with static dimensions.
newtype Matrix v (m :: Nat) (n :: Nat) a = Matrix { toVector :: Vector v (m*n) a }
    deriving (Eq,Show,NFData)

deriving instance (KnownNat m, KnownNat n, Storable x) => Storable (Matrix S.Vector m n x)
deriving instance (KnownNat m, KnownNat n, Numeric x, Num x)
  => Num (Matrix S.Vector m n x)
deriving instance (KnownNat m, KnownNat n, Numeric x, Fractional x)
  => Fractional (Matrix S.Vector m n x)
deriving instance (KnownNat m, KnownNat n, Numeric x, Floating x)
  => Floating (Matrix S.Vector m n x)

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector v n a -> Matrix v n 1 a
{-# INLINE columnVector #-}
columnVector = Matrix

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector v n a -> Matrix v 1 n a
{-# INLINE rowVector #-}
rowVector = Matrix

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: (G.Vector v x, G.Vector v (Vector v n x), KnownNat n) => Vector v m (Vector v n x) -> Matrix v m n x
{-# INLINE fromRows #-}
fromRows = Matrix . concat

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
fromColumns
    :: (G.Vector v x, G.Vector v Int, G.Vector v (Vector v n x), G.Vector v (Vector v m x), KnownNat n, KnownNat m)
    => Vector v n (Vector v m x) -> Matrix v m n x
{-# INLINE fromColumns #-}
fromColumns = transpose . fromRows

-- | The number of rows in the 'Matrix'.
nRows :: forall v m n a . KnownNat m => Matrix v m n a -> Int
{-# INLINE nRows #-}
nRows _ = natValInt (Proxy :: Proxy m)

-- | The number of columns in the 'Matrix'.
nColumns :: forall v m n a . KnownNat n => Matrix v m n a -> Int
{-# INLINE nColumns #-}
nColumns _ = natValInt (Proxy :: Proxy n)

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (G.Vector v a, G.Vector v (Vector v n a), KnownNat n, KnownNat m)
       => Matrix v m n a -> Vector v m (Vector v n a)
{-# INLINE toRows #-}
toRows (Matrix v) = breakEvery v

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns
    :: (G.Vector v a, G.Vector v (Vector v m a), KnownNat m, KnownNat n, G.Vector v Int)
    => Matrix v m n a -> Vector v n (Vector v m a)
{-# INLINE toColumns #-}
toColumns = toRows . transpose


--- BLAS ---


-- | Pure implementation of 'Matrix' transposition.
transpose
    :: forall v m n a . (KnownNat m, KnownNat n, G.Vector v Int, G.Vector v a, G.Vector v (Vector v m a))
    => Matrix v m n a -> Matrix v n m a
{-# INLINE transpose #-}
transpose (Matrix v) =
    let n = natValInt (Proxy :: Proxy n)
     in fromRows $ generate (\j -> generate (\i -> unsafeIndex v $ finiteInt j + finiteInt i*n) :: Vector v m a)

-- | Pure implementation of the dot product.
dotProduct :: (G.Vector v x, Num x) => Vector v n x -> Vector v n x -> x
{-# INLINE dotProduct #-}
dotProduct v1 v2 = weakDotProduct (fromSized v1) (fromSized v2)

-- | Pure implementation of the outer product.
outerProduct
    :: ( KnownNat m, KnownNat n, Num x
       , G.Vector v Int, G.Vector v x, G.Vector v (Vector v n x), G.Vector v (Vector v m x), G.Vector v (Vector v 1 x) )
     => Vector v n x -> Vector v m x -> Matrix v n m x
{-# INLINE outerProduct #-}
outerProduct v1 v2 = matrixMatrixMultiply (columnVector v1) (rowVector v2)

-- | Pure implementation of the dot product on standard vectors.
weakDotProduct :: (G.Vector v x, Num x) => v x -> v x -> x
{-# INLINE weakDotProduct #-}
weakDotProduct v1 v2 = G.foldl foldFun 0 (G.enumFromN 0 (G.length v1) :: S.Vector Int)
    where foldFun d i = d + G.unsafeIndex v1 i * G.unsafeIndex v2 i

-- | Pure 'Matrix' x 'Vector' multiplication.
matrixVectorMultiply
    :: (KnownNat m, KnownNat n, G.Vector v x, G.Vector v (Vector v n x), Num x)
    => Matrix v m n x
    -> Vector v n x
    -> Vector v m x
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply mtx v =
    map (dotProduct v) $ toRows mtx

-- | Pure 'Matrix' x 'Matrix' multiplication.
matrixMatrixMultiply
    :: ( KnownNat m, KnownNat n, KnownNat o, Num x
       , G.Vector v Int, G.Vector v x, G.Vector v (Vector v m x), G.Vector v (Vector v n x), G.Vector v (Vector v o x) )
    => Matrix v m n x
    -> Matrix v n o x
    -> Matrix v m o x
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply mtx1 mtx2 =
    fromColumns . map (matrixVectorMultiply mtx1) $ toColumns mtx2


--- Numeric Classes ---


--instance (Storable x, Numeric x, KnownNat n, KnownNat m)
--  => Num (Matrix S.Vector n m x) where
--    {-# INLINE (+) #-}
--    (+) (Matrix (Vector v1)) (Matrix (Vector v2)) = Matrix $ Vector (H.add v1 v2)
--    {-# INLINE (*) #-}
--    (*) (Matrix xs) (Matrix xs') = Matrix $ xs * xs'
--    {-# INLINE negate #-}
--    negate (Matrix (Vector v)) = Matrix $ Vector (H.scale (-1) v)
--    {-# INLINE abs #-}
--    abs (Matrix xs) = Matrix $ abs xs
--    {-# INLINE signum #-}
--    signum (Matrix xs) = Matrix $ signum xs
--    {-# INLINE fromInteger #-}
--    fromInteger x = Matrix . replicate $ fromInteger x

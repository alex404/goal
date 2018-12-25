{-# LANGUAGE
   DataKinds,
   TypeOperators,
   RankNTypes,
   GADTs,
   ScopedTypeVariables
   #-}
-- | Vectors and Matrices with statically typed dimensions based on boxed vectors.

module Goal.Core.Vector.Boxed
    ( -- * Vector
      module Data.Vector.Sized
      -- ** Blas
    , concat
    , doubleton
    , breakEvery
    , range
    , breakStream
    , toPair
    -- * Matrix
    , Matrix
    -- ** Construction
    , fromRows
    , fromColumns
    , matrixIdentity
    , outerProduct
    , diagonalConcat
    , generateP
    , generatePM
    -- ** Deconstruction
    , toRows
    , toColumns
    , nRows
    , nColumns
    -- ** Manipulation
    , foldr1
    , columnVector
    , rowVector
    -- , diagonalConcat
    -- ** BLAS
    , dotProduct
    -- , determinant
    , matrixVectorMultiply
    , matrixMatrixMultiply
    , inverse
    , transpose
    , convertMatrix
    ) where

--- Imports ---


-- Goal --

import Goal.Core.Util hiding (breakEvery,range)
import qualified Goal.Core.Util (breakEvery)

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S

-- Unqualified --

import Prelude hiding (concat,zipWith,(++),replicate,foldr1)
import GHC.TypeNats
import Data.Proxy

-- Qualified --

import qualified Data.Vector as B
import qualified Data.Vector.Mutable as BM
import qualified Control.Monad.ST as ST
import Data.Vector.Sized hiding (foldr1)
import qualified Data.Vector.Generic.Sized.Internal as I

-- Qualified Imports --

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
concat :: KnownNat n => Vector m (Vector n x) -> Vector (m*n) x
{-# INLINE concat #-}
concat = G.concat

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
doubleton :: x -> x -> Vector 2 x
{-# INLINE doubleton #-}
doubleton = G.doubleton

-- | The number of rows in the 'Matrix'.
nRows :: forall m n a . KnownNat m => Matrix m n a -> Int
{-# INLINE nRows #-}
nRows = G.nRows

-- | The number of columns in the 'Matrix'.
nColumns :: forall m n a . KnownNat n => Matrix m n a -> Int
{-# INLINE nColumns #-}
nColumns = G.nColumns

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (KnownNat m, KnownNat n) => Matrix m n x -> Vector m (Vector n x)
{-# INLINE toRows #-}
toRows = G.toRows

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns :: (KnownNat m, KnownNat n) => Matrix m n x -> Vector n (Vector m x)
{-# INLINE toColumns #-}
toColumns = G.toColumns

-- | Partition of an interval.
range :: (KnownNat n, Fractional x) => x -> x -> Vector n x
{-# INLINE range #-}
range = G.range

-- | Cycles a list of elements and breaks it up into a list of 'Vector's.
breakStream :: forall n a. KnownNat n => [a] -> [Vector n a]
{-# INLINE breakStream #-}
breakStream as =
    I.Vector . B.fromList <$> Goal.Core.Util.breakEvery (natValInt (Proxy :: Proxy n)) (cycle as)
-- | Converts a length two 'Vector' into a pair of elements.
toPair :: Vector 2 x -> (x,x)
{-# INLINE toPair #-}
toPair = G.toPair

-- | foldr1 with a more manageable type signature.
foldr1 :: (KnownNat n, 1 <= n) => (x -> x -> x) -> Vector n x -> x
{-# INLINE foldr1 #-}
foldr1 f v = B.foldr1 f $ fromSized v

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector n a -> Matrix n 1 a
{-# INLINE columnVector #-}
columnVector = G.columnVector

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector n a -> Matrix 1 n a
{-# INLINE rowVector #-}
rowVector = G.rowVector

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: KnownNat n => Vector m (Vector n x) -> Matrix m n x
{-# INLINE fromRows #-}
fromRows = G.fromRows

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromColumns :: (KnownNat n, KnownNat m) => Vector n (Vector m x) -> Matrix m n x
{-# INLINE fromColumns #-}
fromColumns = G.fromColumns

-- | Matrices with static dimensions (boxed).
type Matrix = G.Matrix B.Vector

-- | Diagonally concatenate two matrices, padding the gaps with zeroes (pure implementation).
diagonalConcat
    :: (KnownNat n, KnownNat m, KnownNat o, KnownNat p, Num a)
    => Matrix n m a -> Matrix o p a -> Matrix (n+o) (m+p) a
{-# INLINE diagonalConcat #-}
diagonalConcat mtx1 mtx2 =
    let rws1 = (++ replicate 0) <$> toRows mtx1
        rws2 = (replicate 0 ++) <$> toRows mtx2
     in fromRows $ rws1 ++ rws2

-- | Breaks a 'Vector' into a Vector of Vectors.
breakEvery :: (KnownNat n, KnownNat k) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEvery #-}
breakEvery = G.breakEvery

-- | Pure implementation of the dot product.
dotProduct :: Num x => Vector n x -> Vector n x -> x
{-# INLINE dotProduct #-}
dotProduct = G.dotProduct

-- | Pure implementation of the outer product.
outerProduct
    :: (KnownNat m, KnownNat n, Num x)
    => Vector m x -> Vector n x -> Matrix m n x
{-# INLINE outerProduct #-}
outerProduct = G.outerProduct

-- | Pure implementation of 'Matrix' transposition.
transpose
    :: (KnownNat m, KnownNat n, Num x)
    => Matrix m n x -> Matrix n m x
{-# INLINE transpose #-}
transpose = G.transpose

-- | Pure 'Matrix' x 'Vector' multiplication.
matrixVectorMultiply
    :: (KnownNat m, KnownNat n, Num x)
    => Matrix m n x -> Vector n x -> Vector m x
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply mtx = G.toVector . matrixMatrixMultiply mtx . columnVector

-- | Convert a boxed-based 'Matrix' to a storable-based matrix.
convertMatrix
    :: (KnownNat m, KnownNat n)
    => Matrix m n Double -> S.Matrix m n Double
convertMatrix (G.Matrix v) = G.Matrix $ G.convert v

-- | The identity 'Matrix'.
matrixIdentity :: (KnownNat n, Num a) => Matrix n n a
{-# INLINE matrixIdentity #-}
matrixIdentity =
    fromRows $ generate (\i -> generate (\j -> if finiteInt i == finiteInt j then 1 else 0))

-- | Pure implementation of matrix inversion.
inverse :: forall a n. (Fractional a, Ord a, KnownNat n) => Matrix n n a -> Maybe (Matrix n n a)
{-# INLINE inverse #-}
inverse mtx =
    let rws = fromSized $ fromSized <$> zipWith (++) (toRows mtx) (toRows matrixIdentity)
        n = natValInt (Proxy :: Proxy n)
        rws' = B.foldM' eliminateRow rws $ B.generate n id
     in G.Matrix . I.Vector . B.concatMap (B.drop n) <$> rws'

eliminateRow :: (Ord a, Fractional a) => B.Vector (B.Vector a) -> Int -> Maybe (B.Vector (B.Vector a))
{-# INLINE eliminateRow #-}
eliminateRow mtx k = do
    mtx' <- pivotRow k mtx
    return . nullifyRows k $ normalizePivot k mtx'

pivotRow :: (Fractional a, Ord a) => Int -> B.Vector (B.Vector a) -> Maybe (B.Vector (B.Vector a))
{-# INLINE pivotRow #-}
pivotRow k rws =
    let l = (+k) . B.maxIndex $ abs . flip B.unsafeIndex k . B.take (B.length rws) <$> B.drop k rws
        ak = B.unsafeIndex rws k B.! l
     in if abs ak < 1e-10 then Nothing
                  else ST.runST $ do
                           mrws <- B.thaw rws
                           BM.unsafeSwap mrws k l
                           Just <$> B.freeze mrws

normalizePivot :: Fractional a => Int -> B.Vector (B.Vector a) -> B.Vector (B.Vector a)
{-# INLINE normalizePivot #-}
normalizePivot k rws = ST.runST $ do
    let ak = recip . flip B.unsafeIndex k $ B.unsafeIndex rws k
    mrws <- B.thaw rws
    BM.modify mrws ((*ak) <$>) k
    B.freeze mrws

nullifyRows :: Fractional a => Int -> B.Vector (B.Vector a) -> B.Vector (B.Vector a)
{-# INLINE nullifyRows #-}
nullifyRows k rws =
    let rwk = B.unsafeIndex rws k
        ak = B.unsafeIndex rwk k
        generator i = if i == k then 0 else B.unsafeIndex (B.unsafeIndex rws i) k / ak
        as = B.generate (B.length rws) generator
     in B.zipWith (B.zipWith (-)) rws $ (\a -> (*a) <$> rwk) <$> as

-- | Pure 'Matrix' x 'Matrix' multiplication.
matrixMatrixMultiply
    :: forall m n o a. (KnownNat m, KnownNat n, KnownNat o, Num a)
    => Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply (G.Matrix (I.Vector v)) wm =
    let n = natValInt (Proxy :: Proxy n)
        o = natValInt (Proxy :: Proxy o)
        (G.Matrix (I.Vector w')) = G.transpose wm
        f k = let (i,j) = divMod (finiteInt k) o
                  slc1 = B.unsafeSlice (i*n) n v
                  slc2 = B.unsafeSlice (j*n) n w'
               in G.weakDotProduct slc1 slc2
     in G.Matrix $ G.generate f


-- | Vector generation given based on Proxied Nats. Incorporates a size
-- constraint into the generating function to allow functions which are bounded
-- by the size of the vector.
generateP
    :: forall n x . KnownNat n
    => (forall i j . (KnownNat i, KnownNat j, (i + j + 1) ~ n) => Proxy i -> x)
    -> Vector n x
generateP = G.generateP

-- | Vector generation given based on Proxied Nats (Monadic Version).
generatePM
    :: forall n m x . (KnownNat n, Monad m)
    => (forall i j . (KnownNat i, KnownNat j, (i + j + 1) ~ n) => Proxy i -> m x)
    -> m (Vector n x)
generatePM = G.generatePM



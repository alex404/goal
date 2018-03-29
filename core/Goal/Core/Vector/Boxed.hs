-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLBM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector.Boxed
    ( -- * Vector
      module Data.Vector.Sized
    , BaseVector
      -- ** Blas
    , concat
    , doubleton
    , breakEvery
    , range
    , toPair
    -- * Matrix
    , Matrix
    -- ** Construction
    , fromRows
    , fromColumns
    , matrixIdentity
    , outerProduct
    , diagonalConcat
    -- ** Deconstruction
    , toRows
    , toColumns
    , nRows
    , nColumns
    -- ** Manipulation
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

import Goal.Core.Vector.TypeLits
import qualified Data.Vector as B
import qualified Data.Vector.Mutable as BM
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S

import qualified Control.Monad.ST as ST
import Data.Vector.Sized
import GHC.TypeLits
import Data.Proxy
import qualified Data.Vector.Generic.Sized.Internal as I

import Prelude hiding (concat,zipWith,(++),replicate)

-- Qualified Imports --

type BaseVector = B.Vector

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

-- | The columns of rows in the 'Matrix'.
nColumns :: forall m n a . KnownNat n => Matrix m n a -> Int
{-# INLINE nColumns #-}
nColumns = G.nColumns

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (KnownNat m, KnownNat n) => Matrix m n x -> Vector m (Vector n x)
{-# INLINE toRows #-}
toRows = G.toRows

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toColumns :: (KnownNat m, KnownNat n) => Matrix m n x -> Vector n (Vector m x)
{-# INLINE toColumns #-}
toColumns = G.toColumns

-- | Range
range :: (KnownNat n, Fractional x) => x -> x -> Vector n x
{-# INLINE range #-}
range = G.range

-- | Range
toPair :: Vector 2 x -> (x,x)
{-# INLINE toPair #-}
toPair = G.toPair

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

type Matrix = G.Matrix B.Vector

-- | Diagonally concatenate two matrices, padding the gaps with zeroes.
diagonalConcat
    :: (KnownNat n, KnownNat m, KnownNat o, KnownNat p, Num a)
    => Matrix n m a -> Matrix o p a -> Matrix (n+o) (m+p) a
{-# INLINE diagonalConcat #-}
diagonalConcat mtx1 mtx2 =
    let rws1 = (++ replicate 0) <$> toRows mtx1
        rws2 = (replicate 0 ++) <$> toRows mtx2
     in fromRows $ rws1 ++ rws2

breakEvery :: (KnownNat n, KnownNat k) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEvery #-}
breakEvery = G.breakEvery

dotProduct :: Num x => Vector n x -> Vector n x -> x
{-# INLINE dotProduct #-}
dotProduct = G.dotProduct

outerProduct
    :: (KnownNat m, KnownNat n, Num x)
    => Vector m x -> Vector n x -> Matrix m n x
{-# INLINE outerProduct #-}
outerProduct = G.outerProduct

transpose
    :: (KnownNat m, KnownNat n, Num x)
    => Matrix m n x -> Matrix n m x
{-# INLINE transpose #-}
transpose = G.transpose

matrixVectorMultiply
    :: (KnownNat m, KnownNat n, Num x)
    => Matrix m n x -> Vector n x -> Vector m x
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply mtx = G.toVector . matrixMatrixMultiply mtx . columnVector

convertMatrix
    :: (KnownNat m, KnownNat n)
    => Matrix m n Double -> S.Matrix m n Double
convertMatrix (G.Matrix v) = G.Matrix $ G.convert v

-- | The identity 'Matrix'.
matrixIdentity :: (KnownNat n, Num a) => Matrix n n a
{-# INLINE matrixIdentity #-}
matrixIdentity =
    fromRows $ generate (\i -> generate (\j -> if finiteInt i == finiteInt j then 1 else 0))

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


--deleteRow :: (KnownNat n, KnownNat m, KnownNat k, k <= n-1) => Proxy k -> Matrix n m a -> Matrix (n-1) m a
--deleteRow prxyk mtx =
--    let rws = toRows mtx
--        (hrws,trws) = splitV0 prxyk rws
--     in fromRows . joinV hrws $ tailV trws
--
--deleteColumn :: (KnownNat n, KnownNat m, KnownNat k, k <= m-1, 1 <= m) => Proxy k -> Matrix n m a -> Matrix n (m-1) a
--deleteColumn prxyk mtx =
--    let rws = toColumns mtx
--        (hrws,trws) = splitV0 prxyk rws
--     in fromColumns . joinV hrws $ tailV trws
--
--minorMatrix
--    :: (Num a, KnownNat n, KnownNat i, KnownNat j, i <= n-1, 1 <= i, j <= n-1, 1 <= j, 1 <= n)
--    => Proxy i -> Proxy j ->  Matrix n n a -> Matrix (n-1) (n-1) a
--minorMatrix prxyi prxyj = deleteRow prxyi . deleteColumn prxyj
--
--cofactor
--    :: (Num x, KnownNat n, KnownNat i, KnownNat j, i <= n-1, 1 <= i, j <= n-1, 1 <= j, 1 <= n-1, 1 <= n, 3 <= n)
--    => Proxy i -> Proxy j ->  Matrix n n x -> x
--cofactor prxyi prxyj mtx =
--    let i = natValInt prxyi
--        j = natValInt prxyj
--        mtx' = minorMatrix prxyi prxyj mtx
--     in fromIntegral ((-1)^(i+j)) * determinantV mtx'
--
--determinantV2 :: Num x => B.Vector x -> x
--determinantV2 v =
--    let [a,b,c,d] = B.toList v
--     in a*d - b*c
--
--determinantVN :: Num x => Matrix n n x -> x
--determinantVN v = undefined
--
--determinantV0 :: (KnownNat n, Num x, 1 <= n) => Proxy n -> Matrix n n x -> x
--determinantV0 prxyn mtx =
--    case natValInt prxyn of
--      1 -> B.head . weakVector $ toVector mtx
--      2 -> determinantV2 (weakVector $ toVector mtx)
--      _ -> determinantVN mtx
--
--determinantV :: (KnownNat n, Num x, 1 <= n) => Matrix n n x -> x
--determinantV = determinantV0 Proxy
--
--
---- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
--fromRows :: Vector m (Vector n a) -> Matrix m n a
--{-# INLINE fromRows #-}
--fromRows (Vector vs) = Matrix . Vector $ B.concatMap weakVector vs
--
---- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
--fromColumns :: (KnownNat m, KnownNat n) => Vector n (Vector m a) -> Matrix m n a
--{-# INLINE fromColumns #-}
--fromColumns = matrixTranspose . fromRows
--
---- | The number of rows in the 'Matrix'.
--nRows :: (KnownNat m, KnownNat n) => Matrix m n a -> Int
--{-# INLINE nRows #-}
--nRows = nRows0 Proxy
--
--nRows0 :: KnownNat m => Proxy m -> Matrix m n a -> Int
--{-# INLINE nRows0 #-}
--nRows0 prxyn _ = natValInt prxyn
--
---- | The columns of rows in the 'Matrix'.
--nColumns :: KnownNat n => Matrix m n a -> Int
--{-# INLINE nColumns #-}
--nColumns = nColumns0 Proxy
--
--nColumns0 :: KnownNat n => Proxy n -> Matrix m n a -> Int
--{-# INLINE nColumns0 #-}
--nColumns0 prxym _ = natValInt prxym
--
---- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
--toRows :: (KnownNat m, KnownNat n) => Matrix m n a -> Vector m (Vector n a)
--{-# INLINE toRows #-}
--toRows (Matrix v) = breakEveryV v
--
---- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
--toColumns :: (KnownNat m, KnownNat n) => Matrix m n a -> Vector n (Vector m a)
--{-# INLINE toColumns #-}
--toColumns = toRows . matrixTranspose
--
---- | Transpose a 'Matrix'.
--matrixTranspose :: (KnownNat m, KnownNat n) => Matrix m n a -> Matrix n m a
--{-# INLINE matrixTranspose #-}
--matrixTranspose = matrixTranspose0 Proxy Proxy
--
---- | Invert a 'Matrix' with Gaussian elimination. This is not terribly
---- efficient, but works.
--matrixInverse :: (Fractional a, Ord a, KnownNat n) => Matrix n n a -> Maybe (Matrix n n a)
--{-# INLINE matrixInverse #-}
--matrixInverse = matrixInverse0 Proxy
--
---- | The outer product of two 'Vector's.
--outerProduct :: (KnownNat m, KnownNat n, Num a) => Vector m a -> Vector n a -> Matrix m n a
--{-# INLINE outerProduct #-}
--outerProduct v1 v2 =
--    matrixMatrixMultiply (columnVector v1) (rowVector v2)
--
----- Internal ---
--
--splitV0 :: (KnownNat k) => Proxy k -> Vector n a -> (Vector k a, Vector (n-k) a)
--{-# INLINE splitV0 #-}
--splitV0 prxy (Vector v) =
--    let k = natValInt prxy
--        (v1,v2) = B.splitAt k v
--     in (Vector v1,Vector v2)
--
--replicate0 :: KnownNat n => Proxy n -> a -> Vector n a
--{-# INLINE replicate0 #-}
--replicate0 prxy a = Vector $ B.replicate (natValInt prxy) a
--
--replicateM0 :: (Monad m, KnownNat n) => Proxy n -> m a -> m (Vector n a)
--{-# INLINE replicateM0 #-}
--replicateM0 prxy ma = Vector <$> B.replicateM (natValInt prxy) ma
--
--generateV0 :: KnownNat n => Proxy n -> (Int -> a) -> Vector n a
--{-# INLINE generateV0 #-}
--generateV0 prxy f =
--    let n = natValInt prxy
--     in Vector $ B.generate n f
--
--rangeV0 :: (KnownNat n, 2 <= n, Fractional x) => Proxy n -> x -> x -> Vector n x
--{-# INLINE rangeV0 #-}
--rangeV0 prxy mn mx =
--    let n = natValInt prxy
--        stp = (mx - mn) / fromIntegral (n-1)
--     in Vector $ B.enumFromStepN mn stp n
--
--generateMV0 :: (Monad m, KnownNat n) => Proxy n -> (Int -> m a) -> m (Vector n a)
--{-# INLINE generateMV0 #-}
--generateMV0 prxy f =
--    let n = natValInt prxy
--     in Vector <$> B.generateM n f
--
--
--breakStream :: KnownNat n => [a] -> [Vector n a]
--{-# INLINE breakStream #-}
--breakStream = breakStream0 Proxy
--
--breakStream0 :: KnownNat n => Proxy n -> [a] -> [Vector n a]
--{-# INLINE breakStream0 #-}
--breakStream0 prxyn as =
--    Vector . B.fromList <$> breakEvery (natValInt prxyn) (cycle as)
--
--breakEveryV :: (KnownNat n, KnownNat k) => Vector (n*k) a -> Vector n (Vector k a)
--{-# INLINE breakEveryV #-}
--breakEveryV = breakEveryV0 Proxy Proxy
--
--breakEveryV0 :: (KnownNat n, KnownNat k) => Proxy n -> Proxy k -> Vector (n*k) a -> Vector n (Vector k a)
--{-# INLINE breakEveryV0 #-}
--breakEveryV0 prxyn prxyk (Vector v) =
--    let n = natValInt prxyn
--        k = natValInt prxyk
--     in Vector $ Vector <$> B.generate n (\i -> B.unsafeSlice (i*k) k v)
--


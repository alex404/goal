{-# LANGUAGE UndecidableInstances,GeneralizedNewtypeDeriving,DeriveTraversable #-}

-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLVM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector.Internal where


--- Imports ---


import Control.DeepSeq
import Data.Ratio
import Data.Complex
import Goal.Core.Util (breakEvery)


-- Qualified Imports --

import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Mutable as VM
import qualified Control.Monad.ST as ST
import qualified Numeric.FFT.Vector.Invertible as F


deleteRow :: (KnownNat n, KnownNat m, KnownNat k, k <= n-1) => Proxy k -> Matrix n m a -> Matrix (n-1) m a
deleteRow prxyk mtx =
    let rws = toRows mtx
        (hrws,trws) = splitV0 prxyk rws
     in fromRows . joinV hrws $ tailV trws

deleteColumn :: (KnownNat n, KnownNat m, KnownNat k, k <= m-1, 1 <= m) => Proxy k -> Matrix n m a -> Matrix n (m-1) a
deleteColumn prxyk mtx =
    let rws = toColumns mtx
        (hrws,trws) = splitV0 prxyk rws
     in fromColumns . joinV hrws $ tailV trws

minorMatrix
    :: (Num a, KnownNat n, KnownNat i, KnownNat j, i <= n-1, 1 <= i, j <= n-1, 1 <= j, 1 <= n)
    => Proxy i -> Proxy j ->  Matrix n n a -> Matrix (n-1) (n-1) a
minorMatrix prxyi prxyj = deleteRow prxyi . deleteColumn prxyj

cofactor
    :: (Num x, KnownNat n, KnownNat i, KnownNat j, i <= n-1, 1 <= i, j <= n-1, 1 <= j, 1 <= n-1, 1 <= n, 3 <= n)
    => Proxy i -> Proxy j ->  Matrix n n x -> x
cofactor prxyi prxyj mtx =
    let i = natValInt prxyi
        j = natValInt prxyj
        mtx' = minorMatrix prxyi prxyj mtx
     in fromIntegral ((-1)^(i+j)) * determinantV mtx'

determinantV2 :: Num x => V.Vector x -> x
determinantV2 v =
    let [a,b,c,d] = V.toList v
     in a*d - b*c

determinantVN :: Num x => Matrix n n x -> x
determinantVN v = undefined

determinantV0 :: (KnownNat n, Num x, 1 <= n) => Proxy n -> Matrix n n x -> x
determinantV0 prxyn mtx =
    case natValInt prxyn of
      1 -> V.head . weakVector $ toVector mtx
      2 -> determinantV2 (weakVector $ toVector mtx)
      _ -> determinantVN mtx

determinantV :: (KnownNat n, Num x, 1 <= n) => Matrix n n x -> x
determinantV = determinantV0 Proxy


-- | Diagonally concatenate two matrices, padding the gaps with zeroes.
diagonalConcat
    :: (KnownNat n, KnownNat m, KnownNat o, KnownNat p, Num a)
    => Matrix n m a -> Matrix o p a -> Matrix (n+o) (m+p) a
{-# INLINE diagonalConcat #-}
diagonalConcat mtx1 mtx2 =
    let rws1 = flip joinV (replicateV 0) <$> toRows mtx1
        rws2 = joinV (replicateV 0) <$> toRows mtx2
     in fromRows $ joinV rws1 rws2

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: Vector m (Vector n a) -> Matrix m n a
{-# INLINE fromRows #-}
fromRows (Vector vs) = Matrix . Vector $ V.concatMap weakVector vs

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
fromColumns :: (KnownNat m, KnownNat n) => Vector n (Vector m a) -> Matrix m n a
{-# INLINE fromColumns #-}
fromColumns = matrixTranspose . fromRows

-- | The number of rows in the 'Matrix'.
nRows :: (KnownNat m, KnownNat n) => Matrix m n a -> Int
{-# INLINE nRows #-}
nRows = nRows0 Proxy

nRows0 :: KnownNat m => Proxy m -> Matrix m n a -> Int
{-# INLINE nRows0 #-}
nRows0 prxyn _ = natValInt prxyn

-- | The columns of rows in the 'Matrix'.
nColumns :: KnownNat n => Matrix m n a -> Int
{-# INLINE nColumns #-}
nColumns = nColumns0 Proxy

nColumns0 :: KnownNat n => Proxy n -> Matrix m n a -> Int
{-# INLINE nColumns0 #-}
nColumns0 prxym _ = natValInt prxym

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of rows.
toRows :: (KnownNat m, KnownNat n) => Matrix m n a -> Vector m (Vector n a)
{-# INLINE toRows #-}
toRows (Matrix v) = breakEveryV v

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns :: (KnownNat m, KnownNat n) => Matrix m n a -> Vector n (Vector m a)
{-# INLINE toColumns #-}
toColumns = toRows . matrixTranspose

-- | Transpose a 'Matrix'.
matrixTranspose :: (KnownNat m, KnownNat n) => Matrix m n a -> Matrix n m a
{-# INLINE matrixTranspose #-}
matrixTranspose = matrixTranspose0 Proxy Proxy

-- | Invert a 'Matrix' with Gaussian elimination. This is not terribly
-- efficient, but works.
matrixInverse :: (Fractional a, Ord a, KnownNat n) => Matrix n n a -> Maybe (Matrix n n a)
{-# INLINE matrixInverse #-}
matrixInverse = matrixInverse0 Proxy

-- | The identity 'Matrix'.
matrixIdentity :: (KnownNat n, Num a) => Matrix n n a
{-# INLINE matrixIdentity #-}
matrixIdentity =
    fromRows $ generateV (\i -> generateV (\j -> if i == j then 1 else 0))

-- | The outer product of two 'Vector's.
outerProduct :: (KnownNat m, KnownNat n, Num a) => Vector m a -> Vector n a -> Matrix m n a
{-# INLINE outerProduct #-}
outerProduct v1 v2 =
    matrixMatrixMultiply (columnVector v1) (rowVector v2)

--- Internal ---

ratVal0 :: (KnownNat n, KnownNat d) => Proxy n -> Proxy d -> Proxy (n / d) -> Rational
ratVal0 prxyn prxyd _ = natVal prxyn % natVal prxyd

matrixInverse0 :: (Fractional a, Ord a, KnownNat n) => Proxy n -> Matrix n n a -> Maybe (Matrix n n a)
{-# INLINE matrixInverse0 #-}
matrixInverse0 prxyn mtx =
    let rws = weakVector $ weakVector <$> zipWithV joinV (toRows mtx) (toRows matrixIdentity)
        n = natValInt prxyn
        rws' = V.foldM' eliminateRow rws $ V.generate n id
     in Matrix . Vector . V.concatMap (V.drop n) <$> rws'

eliminateRow :: (Ord a, Fractional a) => V.Vector (V.Vector a) -> Int -> Maybe (V.Vector (V.Vector a))
{-# INLINE eliminateRow #-}
eliminateRow mtx k = do
    mtx' <- pivotRow k mtx
    return . nullifyRows k $ normalizePivot k mtx'

pivotRow :: (Fractional a, Ord a) => Int -> V.Vector (V.Vector a) -> Maybe (V.Vector (V.Vector a))
{-# INLINE pivotRow #-}
pivotRow k rws =
    let l = (+k) . V.maxIndex $ abs . flip V.unsafeIndex k . V.take (V.length rws) <$> V.drop k rws
        ak = V.unsafeIndex rws k V.! l
     in if abs ak < 1e-10 then Nothing
                  else ST.runST $ do
                           mrws <- V.thaw rws
                           VM.unsafeSwap mrws k l
                           Just <$> V.freeze mrws

normalizePivot :: Fractional a => Int -> V.Vector (V.Vector a) -> V.Vector (V.Vector a)
{-# INLINE normalizePivot #-}
normalizePivot k rws = ST.runST $ do
    let ak = recip . flip V.unsafeIndex k $ V.unsafeIndex rws k
    mrws <- V.thaw rws
    VM.modify mrws ((*ak) <$>) k
    V.freeze mrws

nullifyRows :: Fractional a => Int -> V.Vector (V.Vector a) -> V.Vector (V.Vector a)
{-# INLINE nullifyRows #-}
nullifyRows k rws =
    let rwk = V.unsafeIndex rws k
        ak = V.unsafeIndex rwk k
        generator i = if i == k then 0 else V.unsafeIndex (V.unsafeIndex rws i) k / ak
        as = V.generate (V.length rws) generator
     in V.zipWith (V.zipWith (-)) rws $ (\a -> (*a) <$> rwk) <$> as

splitV0 :: (KnownNat k) => Proxy k -> Vector n a -> (Vector k a, Vector (n-k) a)
{-# INLINE splitV0 #-}
splitV0 prxy (Vector v) =
    let k = natValInt prxy
        (v1,v2) = V.splitAt k v
     in (Vector v1,Vector v2)

replicate0 :: KnownNat n => Proxy n -> a -> Vector n a
{-# INLINE replicate0 #-}
replicate0 prxy a = Vector $ V.replicate (natValInt prxy) a

replicateM0 :: (Monad m, KnownNat n) => Proxy n -> m a -> m (Vector n a)
{-# INLINE replicateM0 #-}
replicateM0 prxy ma = Vector <$> V.replicateM (natValInt prxy) ma

generateV0 :: KnownNat n => Proxy n -> (Int -> a) -> Vector n a
{-# INLINE generateV0 #-}
generateV0 prxy f =
    let n = natValInt prxy
     in Vector $ V.generate n f

rangeV0 :: (KnownNat n, 2 <= n, Fractional x) => Proxy n -> x -> x -> Vector n x
{-# INLINE rangeV0 #-}
rangeV0 prxy mn mx =
    let n = natValInt prxy
        stp = (mx - mn) / fromIntegral (n-1)
     in Vector $ V.enumFromStepN mn stp n

generateMV0 :: (Monad m, KnownNat n) => Proxy n -> (Int -> m a) -> m (Vector n a)
{-# INLINE generateMV0 #-}
generateMV0 prxy f =
    let n = natValInt prxy
     in Vector <$> V.generateM n f

weakDotProduct :: Num a => V.Vector a -> V.Vector a -> a
{-# INLINE weakDotProduct #-}
weakDotProduct v1 v2 = V.foldl foldFun 0 (V.enumFromN 0 (V.length v1))
    where foldFun d i = d + V.unsafeIndex v1 i * V.unsafeIndex v2 i

breakStream :: KnownNat n => [a] -> [Vector n a]
{-# INLINE breakStream #-}
breakStream = breakStream0 Proxy

breakStream0 :: KnownNat n => Proxy n -> [a] -> [Vector n a]
{-# INLINE breakStream0 #-}
breakStream0 prxyn as =
    Vector . V.fromList <$> breakEvery (natValInt prxyn) (cycle as)

breakEveryV :: (KnownNat n, KnownNat k) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEveryV #-}
breakEveryV = breakEveryV0 Proxy Proxy

breakEveryV0 :: (KnownNat n, KnownNat k) => Proxy n -> Proxy k -> Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEveryV0 #-}
breakEveryV0 prxyn prxyk (Vector v) =
    let n = natValInt prxyn
        k = natValInt prxyk
     in Vector $ Vector <$> V.generate n (\i -> V.unsafeSlice (i*k) k v)

matrixTranspose0 :: (KnownNat m, KnownNat n) => Proxy m -> Proxy n -> Matrix m n a -> Matrix n m a
{-# INLINE matrixTranspose0 #-}
matrixTranspose0 prxym prxyn (Matrix (Vector v)) =
    let m = natValInt prxym
        n = natValInt prxyn
        vi = V.concatMap (\i -> V.generate m (\j -> i + j*n)) $ V.generate n id
     in Matrix . Vector $ V.unsafeBackpermute v vi

matrixMatrixMultiply0
    :: (KnownNat m, KnownNat n, KnownNat o, Num a)
    => Proxy n -> Proxy o -> Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply0 #-}
matrixMatrixMultiply0 prxyn prxyo (Matrix (Vector v)) wm =
    let n = natValInt prxyn
        o = natValInt prxyo
        (Matrix (Vector w')) = matrixTranspose wm
        f k = let (i,j) = divMod k o
                  slc1 = V.unsafeSlice (i*n) n v
                  slc2 = V.unsafeSlice (j*n) n w'
               in weakDotProduct slc1 slc2
     in Matrix $ generateV f

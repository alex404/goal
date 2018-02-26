{-# LANGUAGE BangPatterns,GeneralizedNewtypeDeriving,RankNTypes,TypeOperators,DeriveTraversable,FlexibleContexts,DataKinds,TypeFamilies #-}

-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLVM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector
    ( -- * Vector
    Vector
      -- ** Construction
    , strongVector
    , strongVectorFromList
    , empty
    , singleton
    , doubleton
    , replicateV
    , generateV
    , rangeV
    -- *** Monadic
    , replicateMV
    , generateMV
      -- ** Deconstruction
    , headV
    , tailV
    , headTail
    , toPair
    , splitV
    , lookupV
    -- ** Concatenation
    , (&)
    , snoc
    , joinV
    , flattenV
    -- ** Manipulation
    , reverseV
    , zipV
    , unzipV
    , breakEveryV
    -- ** Computation
    , maxIndexV
    , zipWithV
    , foldlV'
    , foldl1V'
    , scanl1V'
    , dotProduct
    , deconvolve
    , deconvolve'
    -- ** Streaming
    , breakStream
    -- * Matrix
    , Matrix (Matrix,toVector)
    , matrixTranspose
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
    -- ** Computation
    , matrixVectorMultiply
    , matrixMatrixMultiply
    , matrixInverse
      -- * Miscellaneous
    , natValInt
    , prettyPrintMatrix
    -- * Type Rationals
    , Rat
    , type (/)
    , ratVal
    ) where


--- Imports ---


import GHC.TypeLits
import Data.Proxy
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


lookupV :: Vector n a -> Int -> Maybe a
{-# INLINE lookupV #-}
lookupV (Vector v) k = v V.!? k


--- BLAS ---


-- | The dot product of two numerical 'Vector's.
dotProduct :: Num a => Vector n a -> Vector n a -> a
{-# INLINE dotProduct #-}
dotProduct (Vector v1) (Vector v2) = weakDotProduct v1 v2
--
-- | Apply a linear transformation to a 'Vector'.
matrixVectorMultiply :: (KnownNat m, KnownNat n, Num a)
                     => Matrix m n a -> Vector n a -> Vector m a
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply m v =
    let Matrix v' = matrixMatrixMultiply m $ columnVector v
     in v'

-- | Multiply a 'Matrix' with a second 'Matrix'.
matrixMatrixMultiply :: (KnownNat m, KnownNat n, KnownNat o, Num a) => Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply = matrixMatrixMultiply0 Proxy Proxy

-- | The dot product of two numerical 'Vector's.
--dotProduct :: Num a => Vector n a -> Vector n a -> a
--{-# INLINE dotProduct #-}
--dotProduct v1 v2 = sum $ zipWithV (*) v1 v2
--
---- | The dot product of two numerical 'Vector's.
--matrixVectorMultiply
--    :: (KnownNat m, KnownNat n, Num a)
--    => Matrix m n a
--    -> Vector n a
--    -> Vector m a
--{-# INLINE matrixVectorMultiply #-}
--matrixVectorMultiply mtx v = dotProduct v <$> toRows mtx
--
---- | The dot product of two numerical 'Vector's.
--matrixMatrixMultiply
--    :: (KnownNat m, KnownNat n, KnownNat o, Num a)
--    => Matrix m n a
--    -> Matrix n o a
--    -> Matrix m o a
--{-# INLINE matrixMatrixMultiply #-}
--matrixMatrixMultiply mtx1 mtx2 = fromColumns $ matrixVectorMultiply mtx1 <$> toColumns mtx2

--- Type Rations ---


-- | Type level rational numbers. This implementation does not currently permit negative numbers.
data Rat (n :: Nat) (d :: Nat)

-- | Infix 'Rat'.
type (/) n d = Rat n d

-- | Recover a rational value from a 'Proxy'.
ratVal :: (KnownNat n, KnownNat d) => Proxy (n / d) -> Rational
ratVal = ratVal0 Proxy Proxy


--- Datatypes ---

-- | A newtype wrapper around Data.Vector with static length checks.
newtype Vector (n :: Nat) a = Vector {weakVector :: V.Vector a} deriving (Eq,Show,Functor,Foldable,Traversable,NFData)

vectorize0 :: (G.Vector v x, KnownNat n) => Proxy n -> v x -> Vector n x
vectorize0 prxyn v =
    let n = natValInt prxyn
     in if n == G.length v
           then Vector $ V.convert v
           else error "Vector Length Mismatch"

strongVector :: (G.Vector v x, KnownNat n) => v x -> Vector n x
strongVector = vectorize0 Proxy

strongVectorFromList :: KnownNat n => [x] -> Vector n x
strongVectorFromList = vectorize0 Proxy . V.fromList


{-
readPrec0 :: KnownNat n => Proxy n -> T.ReadPrec (V.Vector x) -> T.ReadPrec (Vector n x)
readPrec0 prxyn readPrec' = do
    let n = natValInt prxyn
    n' <- V.length <$> readPrec'
    if n' == n
       then Vector <$> readPrec'
       else fail "Vector Length Mismatch"

instance (KnownNat n, Read (V.Vector x)) => Read (Vector n x) where
    readPrec = readPrec0 Proxy readPrec
    -}

-- | The empty vector.
empty :: Vector 0 a
{-# INLINE empty #-}
empty = Vector V.empty

-- | A vector singleton.
singleton :: a -> Vector 1 a
{-# INLINE singleton #-}
singleton = Vector . V.singleton

-- | A vector doubleton.
doubleton :: a -> a -> Vector 2 a
{-# INLINE doubleton #-}
doubleton a1 a2 = Vector $ V.fromList [a1,a2]

-- | Recover an 'Int' (as opposed to 'Integer') from a Proxy.
natValInt :: forall n proxy . KnownNat n => proxy n -> Int
{-# INLINE natValInt #-}
natValInt = fromIntegral . natVal

-- | Infix cons for 'Vector's.
(&) :: a -> Vector n a -> Vector (n+1) a
{-# INLINE (&) #-}
(&) a (Vector v) = Vector $ V.cons a v
infixr 1 &

-- | Append for 'Vector's.
snoc :: Vector n a -> a -> Vector (n+1) a
{-# INLINE snoc #-}
snoc (Vector v) a = Vector $ V.snoc v a

-- | Concatenate two 'Vector's.
joinV :: Vector m a -> Vector n a -> Vector (n+m) a
{-# INLINE joinV #-}
joinV (Vector v) (Vector w) = Vector $ v V.++ w

-- | Flatten a 'Vector' of 'Vector's.
flattenV :: Vector m (Vector n a) -> Vector (n * m) a
{-# INLINE flattenV #-}
flattenV (Vector vs) = Vector . V.concatMap id $ weakVector <$> vs

-- | Return the first element of a 'Vector'.
headV :: Vector n a -> a
{-# INLINE headV #-}
headV (Vector v) = V.unsafeIndex v 0

-- | Return all but the first element of a 'Vector'.
tailV :: (1 <= n) => Vector n a ->  Vector (n-1) a
{-# INLINE tailV #-}
tailV (Vector v) = Vector $ V.tail v

-- | Return the head and tail of a 'Vector'.
headTail :: (1 <= n) => Vector n a -> (a, Vector (n-1) a)
{-# INLINE headTail #-}
headTail (Vector v) =
    let (a,v') = V.splitAt 1 v
     in (V.head a,Vector v')

-- | Convert a length 2 'Vector' into a pair of elements.
toPair :: Vector 2 a -> (a,a)
{-# INLINE toPair #-}
toPair (Vector v) = (V.unsafeIndex v 0, V.unsafeIndex v 1)

-- | Split a 'Vector' into two 'Vector's.
splitV :: (KnownNat k, KnownNat (n-k)) => Vector n a -> (Vector k a, Vector (n-k) a)
{-# INLINE splitV #-}
splitV = splitV0 Proxy

-- | Replicate a single element into a 'Vector'.
replicateV :: KnownNat n => a -> Vector n a
{-# INLINE replicateV #-}
replicateV = replicate0 Proxy

-- | Repeat an action to generate a 'Vector'.
replicateMV :: (Monad m, KnownNat n) => m a -> m (Vector n a)
{-# INLINE replicateMV #-}
replicateMV = replicateM0 Proxy

-- | Generate a 'Vector' by stepping uniformly from the first input to the
-- second input.
rangeV :: (KnownNat n, 2 <= n, Fractional x) => x -> x -> Vector n x
{-# INLINE rangeV  #-}
rangeV = rangeV0 Proxy

-- | Generate a 'Vector' as a function of the element index.
generateV :: KnownNat n => (Int -> a) -> Vector n a
{-# INLINE generateV  #-}
generateV = generateV0 Proxy

-- | Generate a 'Vector' with actions that are functions of the element index.
generateMV :: (Monad m, KnownNat n) => (Int -> m a) -> m (Vector n a)
{-# INLINE generateMV  #-}
generateMV = generateMV0 Proxy

-- | Reverse a 'Vector'.
reverseV :: Vector n a -> Vector n a
{-# INLINE reverseV  #-}
reverseV (Vector v) = Vector $ V.reverse v

-- | Zip two 'Vector's together.
zipV :: Vector n a -> Vector n b -> Vector n (a,b)
{-# INLINE zipV #-}
zipV (Vector va) (Vector vb) = Vector $ V.zip va vb

-- | Unzip two 'Vector's.
unzipV :: Vector n (a,b) -> (Vector n a,Vector n b)
{-# INLINE unzipV #-}
unzipV (Vector vab) =
    let (va,vb) = V.unzip vab
     in (Vector va, Vector vb)

-- | Index of the maximum element of the 'Vector'.
maxIndexV :: (Ord a, 1 <= n) => Vector n a -> Int
{-# INLINE maxIndexV #-}
maxIndexV (Vector v) = V.maxIndex v

-- | Zip two 'Vector's together with a function.
zipWithV :: (a -> b -> c) -> Vector n a -> Vector n b -> Vector n c
{-# INLINE zipWithV #-}
zipWithV f (Vector va) (Vector vb) = Vector $ V.zipWith f va vb

-- | Strictly left-fold the elements of the 'Vector' into a value.
foldlV' :: (a -> b -> a) -> a -> Vector n b -> a
{-# INLINE foldlV' #-}
foldlV' f a (Vector v) = V.foldl' f a v

-- | Strictly left-fold the elements of the 'Vector' into a value, where the
-- initial element is given by the head of the 'Vector'.
foldl1V' :: (a -> a -> a) -> Vector n a -> a
{-# INLINE foldl1V' #-}
foldl1V' f (Vector v) = V.foldl1' f v

-- | Strictly left-scan the elements of the 'Vector' into a 'Vector'.
scanl1V' :: (a -> a -> a) -> Vector n a -> Vector n a
{-# INLINE scanl1V' #-}
scanl1V' f (Vector v) = Vector $ V.scanl1' f v

-- | Deconvolves the signal with the convolutional kernel.
deconvolve
    :: Vector n Double -- ^ Signal
    -> Vector n Double -- ^ Kernel
    -> Vector n Double -- ^ Deconvolved Signal
{-# INLINE deconvolve #-}
deconvolve (Vector v1) (Vector v2) =
    let v1' = F.run F.dft $ (:+ 0) <$> v1
        v2' = F.run F.dft $ (:+ 0) <$> v2
     in Vector . fmap realPart . F.run F.idft $ V.zipWith (/) v1' v2'

-- | Deconvolves the signal with the convolutional kernel.
deconvolve'
    :: Vector n Double -- ^ Signal
    -> Vector n Double -- ^ Kernel
    -> Vector n Double -- ^ Deconvolved Signal
{-# INLINE deconvolve' #-}
deconvolve' (Vector v1) (Vector v2) =
    let v1' = F.run F.dftR2C v1
        v2' = F.run F.dftR2C v2
     in Vector . F.run F.dftC2R $ V.zipWith (/) v1' v2'


--- Matrices ---


-- | Prety print the values of a 'Matrix' (for extremely simple values of pretty).
prettyPrintMatrix :: (KnownNat m, KnownNat n, Show a) => Matrix m n a -> IO ()
prettyPrintMatrix mtx = do
    let rws = toRows mtx
    sequence_ $ print . V.toList . weakVector <$> rws

-- | Matrices with static dimensions.
newtype Matrix (m :: Nat) (n :: Nat) a = Matrix { toVector :: Vector (m*n) a } deriving (Eq,Show,Functor,Foldable,Traversable,NFData)

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector n a -> Matrix n 1 a
{-# INLINE columnVector #-}
columnVector = Matrix

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector n a -> Matrix 1 n a
{-# INLINE rowVector #-}
rowVector = Matrix

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
eliminateRow !mtx k = do
    mtx' <- pivotRow k mtx
    return . nullifyRows k $ normalizePivot k mtx'

pivotRow :: (Fractional a, Ord a) => Int -> V.Vector (V.Vector a) -> Maybe (V.Vector (V.Vector a))
{-# INLINE pivotRow #-}
pivotRow k !rws =
    let l = (+k) . V.maxIndex $ abs . flip V.unsafeIndex k . V.take (V.length rws) <$> V.drop k rws
        ak = V.unsafeIndex rws k V.! l
     in if abs ak < 1e-10 then Nothing
                  else ST.runST $ do
                           mrws <- V.thaw rws
                           VM.unsafeSwap mrws k l
                           Just <$> V.freeze mrws

normalizePivot :: Fractional a => Int -> V.Vector (V.Vector a) -> V.Vector (V.Vector a)
{-# INLINE normalizePivot #-}
normalizePivot k !rws = ST.runST $ do
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

{-
matrixMatrixMultiply0' :: (KnownNat m, KnownNat n, KnownNat o, Num a)
                       => Proxy m -> Proxy n -> Proxy o -> Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply0' #-}
matrixMatrixMultiply0' prxym prxyn prxyo (Matrix (Vector v)) wm =
    let m = natValInt prxym
        n = natValInt prxyn
        o = natValInt prxyo
        (Matrix (Vector w')) = matrixTranspose wm
        vi = V.concatMap (\i -> V.generate o (\j -> (i,j))) $ V.generate m id
        g (i,j) = let slc1 = V.unsafeSlice (i*n) n v
                      slc2 = V.unsafeSlice (j*n) n w'
                   in weakDotProduct slc1 slc2
     in Matrix . Vector $ V.map g vi

matrixMatrixMultiply0' :: (KnownNat n, KnownNat o, KnownNat (o*n), KnownNat (m*o), Num a)
                      => Proxy n -> Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply0' #-}
matrixMatrixMultiply0' prxyn (Matrix (Vector v)) wm =
    let n = natValInt prxyn
        (Matrix (Vector w)) = matrixTranspose wm
        vs = breakEveryV n v
        ws = breakEveryV n w
        f vr = sum . V.zipWith (*) vr <$> ws
     in Matrix (Vector $ V.concatMap f vs)

matrixMatrixMultiply' :: (KnownNat n, KnownNat o, KnownNat (o*n), KnownNat (m*o), Num a)
                     => Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply' #-}
matrixMatrixMultiply' = matrixMatrixMultiply0' Proxy

weakDotProduct :: Num a => V.Vector a -> V.Vector a -> a
{-# INLINE weakDotProduct #-}
weakDotProduct v1 v2 = numLoopFold 0 (V.length v1 - 1) 0 $
  \r i -> V.unsafeIndex v1 i * V.unsafeIndex v2 i + r
-}
{-
matrixMatrixMultiply0' :: (KnownNat n, KnownNat o, KnownNat (o*n), KnownNat (m*o), Num a)
                      => Proxy n -> Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply0' #-}
matrixMatrixMultiply0' prxyn (Matrix (Vector v)) wm =
    let n = natValInt prxyn
        (Matrix (Vector w)) = matrixTranspose wm
        vs = breakEveryV n v
        ws = breakEveryV n w
        f vr = sum . V.zipWith (*) vr <$> ws
     in Matrix (Vector $ V.concatMap f vs)

matrixMatrixMultiply' :: (KnownNat n, KnownNat o, KnownNat (o*n), KnownNat (m*o), Num a)
                     => Matrix m n a -> Matrix n o a -> Matrix m o a
{-# INLINE matrixMatrixMultiply' #-}
matrixMatrixMultiply' = matrixMatrixMultiply0' Proxy

weakDotProduct :: Num a => V.Vector a -> V.Vector a -> a
{-# INLINE weakDotProduct #-}
weakDotProduct v1 v2 = numLoopFold 0 (V.length v1 - 1) 0 $
  \r i -> V.unsafeIndex v1 i * V.unsafeIndex v2 i + r

weakDotProduct :: Num a => V.Vector a -> V.Vector a -> a
{-# INLINE weakDotProduct #-}
weakDotProduct v1 v2 = dotProductHelper (V.length v1 - 1) 0
    where dotProductHelper 0 !d = d + V.unsafeIndex v1 0 * V.unsafeIndex v2 0
          dotProductHelper !k !d =
              let d0 = V.unsafeIndex v1 k * V.unsafeIndex v2 k
               in dotProductHelper (k-1) (d0 + d)


-}

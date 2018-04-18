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
    , BaseVector
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
    , fromHMatrix
    -- ** Deconstruction
    , toRows
    , toColumns
    , nRows
    , nColumns
    , toHMatrix
    -- ** Manipulation
    , columnVector
    , rowVector
    , diagonalConcat
    , foldr1
    , zipFold
    -- ** BLAS
    , add
    , scale
    , average
    , dotProduct
    , determinant
    , matrixVectorMultiply
    , matrixMatrixMultiply
    , inverse
    , transpose
    -- ** Convolutions
    , crossCorrelate2d'
    , crossCorrelate2d
    , convolve2d
    , kernelDifferential
    -- * Miscellaneous
    , prettyPrintMatrix
    ) where


--- Imports ---


import GHC.TypeLits
import Data.Proxy
import Foreign.Storable
import Goal.Core.Vector.TypeLits
import Data.Vector.Storable.Sized hiding (foldr1)
import Numeric.LinearAlgebra (Field,Numeric)

-- Qualified Imports --

import qualified Data.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G
import qualified Data.Vector.Generic.Sized.Internal as G
import qualified Numeric.LinearAlgebra as H
import qualified Data.List as L

import Prelude hiding (concat,foldr1,concatMap,replicate,(++),reverse,map,length,zipWith)


--- Generic ---

type BaseVector = S.Vector

-- | Matrices with static dimensions.
type Matrix = G.Matrix S.Vector

zipFold :: (KnownNat n, Storable x, Storable y) => (z -> x -> y -> z) -> z -> Vector n x -> Vector n y -> z
{-# INLINE zipFold #-}
zipFold f z0 xs ys =
    let n = length xs
        foldfun z i = f z (unsafeIndex xs i) (unsafeIndex ys i)
     in L.foldl' foldfun z0 [0..n-1]

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
foldr1 :: (KnownNat n, 1 <= n, Storable x) => (x -> x -> x) -> Vector n x -> x
{-# INLINE foldr1 #-}
foldr1 f (G.Vector xs) = S.foldr1 f xs

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
concat :: (KnownNat n, Storable x) => Vector m (Vector n x) -> Vector (m*n) x
{-# INLINE concat #-}
concat = G.concat

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
doubleton :: Storable x => x -> x -> Vector 2 x
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
toRows :: (KnownNat m, KnownNat n, Storable x) => Matrix m n x -> Vector m (Vector n x)
{-# INLINE toRows #-}
toRows = G.toRows

-- | Turn a 'Vector' into a single column 'Matrix'.
columnVector :: Vector n a -> Matrix n 1 a
{-# INLINE columnVector #-}
columnVector = G.columnVector

-- | Turn a 'Vector' into a single row 'Matrix'.
rowVector :: Vector n a -> Matrix 1 n a
{-# INLINE rowVector #-}
rowVector = G.rowVector

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the rows.
fromRows :: (KnownNat n, Storable x) => Vector m (Vector n x) -> Matrix m n x
{-# INLINE fromRows #-}
fromRows = G.fromRows

-- | Range
range :: (KnownNat n, Fractional x, Storable x) => x -> x -> Vector n x
{-# INLINE range #-}
range = G.range

-- | Range
toPair :: Storable x => Vector 2 x -> (x,x)
{-# INLINE toPair #-}
toPair = G.toPair


--- HMatrix ---


toHMatrix :: forall m n x . (KnownNat n, Storable x) => Matrix m n x -> H.Matrix x
{-# INLINE toHMatrix #-}
toHMatrix (G.Matrix mtx) =
    H.reshape (natValInt (Proxy :: Proxy n)) $ fromSized mtx

fromHMatrix :: Numeric x => H.Matrix x -> Matrix m n x
{-# INLINE fromHMatrix #-}
fromHMatrix = G.Matrix . G.Vector . H.flatten

-- | Convert a 'Matrix' into a 'Vector' of 'Vector's of columns.
toColumns :: (KnownNat m, KnownNat n, Numeric x) => Matrix m n x -> Vector n (Vector m x)
{-# INLINE toColumns #-}
toColumns = toRows . transpose

-- | Create a 'Matrix' from a 'Vector' of 'Vector's which represent the columns.
fromColumns :: (KnownNat m, KnownNat n, Numeric x) => Vector n (Vector m x) -> Matrix m n x
{-# INLINE fromColumns #-}
fromColumns = transpose . fromRows

breakEvery :: forall n k a . (KnownNat n, KnownNat k, Storable a) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> G.Vector $ S.unsafeSlice (finiteInt i*k) k v)


--- BLAS ---


-- | The dot product of two numerical 'Vector's.
add :: Numeric x => Vector n x -> Vector n x -> Vector n x
{-# INLINE add #-}
add (G.Vector v1) (G.Vector v2) = G.Vector (H.add v1 v2)

-- | The dot product of two numerical 'Vector's.
scale :: Numeric x => x -> Vector n x -> Vector n x
{-# INLINE scale #-}
scale x (G.Vector v) = G.Vector (H.scale x v)

-- | The dot product of two numerical 'Vector's.
average :: (Numeric x, Fractional x) => Vector n x -> x
{-# INLINE average #-}
average (G.Vector v) = H.sumElements v / fromIntegral (S.length v)

-- | The dot product of two numerical 'Vector's.
dotProduct :: Numeric x => Vector n x -> Vector n x -> x
{-# INLINE dotProduct #-}
dotProduct v1 v2 = H.dot (fromSized v1) (fromSized v2)

-- | The dot product of two numerical 'Vector's.
determinant :: (KnownNat n, Field x) => Matrix n n x -> x
{-# INLINE determinant #-}
determinant = H.det . toHMatrix

-- | Transpose a 'Matrix'.
transpose :: forall m n x . (KnownNat m, KnownNat n, Numeric x) => Matrix m n x -> Matrix n m x
{-# INLINE transpose #-}
transpose (G.Matrix mtx) =
    G.Matrix $ withVectorUnsafe (H.flatten . H.tr . H.reshape (natValInt (Proxy :: Proxy n))) mtx

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
inverse (G.Matrix mtx) =
    G.Matrix $ withVectorUnsafe (H.flatten . H.inv . H.reshape (natValInt (Proxy :: Proxy n))) mtx

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
    G.Vector $ toHMatrix mtx H.#> fromSized v

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


--- Convolutions ---


--convolve2d
--    :: (KnownNat kr, KnownNat kc, KnownNat mr, KnownNat mc, Num)
--    => Matrix kr kc x
--    -> Matrix mr mc x
--    -> Matrix mr mc x
--{-# INLINE convolve2d #-}
--convolve2d krn mtx = undefined
--    let generator idx =
--            let (r,c) = divMod idx oc
--             in matrixVectorMultiply kmtx $! cutWindow cnv Proxy Proxy Proxy Proxy img r c
--     in flattenV $! generateV generator

--toKernelMatrix
--    :: (KnownNat kr, KnownNat kc, KnownNat mr, KnownNat mc)
--    => Matrix kr kc
--    => Matrix mr mc
--    => Point c (Convolutional rd r c ih oh om im) x -> Matrix nk ((2*rd+1) * (2*rd+1) * nd) x
--{-# INLINE toKernelMatrix #-}
--toKernelMatrix prms = Matrix $ coordinates prms

--matrixIndex
--    :: forall m n x . (KnownNat m, KnownNat n, Storable x)
--    => Matrix m n x -> Integer -> Integer -> Maybe x
--matrixIndex (G.Matrix v) r c0 = do
--    let n = natVal (Proxy :: Proxy n)
--    c <- if c0 < n then Just c0 else Nothing
--    idx <- F.packFinite $ r*n + c
--    return $ index v idx

--windowIndices
--    :: forall rdkr rdkc mr mc x . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, Num x, Storable x)
--    => Matrix (2*rdkr+1) (2*rdkc+1) x
--    -> Matrix mr mc x
--    -> Int
--    -> Int
--    -> Vector ((2*rdkr+1)*(2*rdkc+1)) Int
--{-# INLINE windowIndices #-}
--windowIndices _ _ r c =
--    let rdkr = natValInt (Proxy :: Proxy rdkr)
--        rdkc = natValInt (Proxy :: Proxy rdkc)
--        mc = natValInt (Proxy :: Proxy mc)
--        kc = (2*rdkr + 1)
--        reIndex idx =
--            let (ir,ic) = divMod (fromIntegral idx) kc
--             in (ir + r - rdkr) * mc + (ic + c - rdkc)
--     in generate reIndex

kernelIndices'
    :: forall rdkc mr mc md . (KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md)
    => Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Proxy md
    -> Int
    -> Int
    -> Vector (mr*mc*md) Int
{-# INLINE kernelIndices' #-}
kernelIndices' prdkc pmc pmr pmd kr kc =
    let rdkc = natValInt prdkc
        mr = natValInt pmr
        mc = natValInt pmc
        md = natValInt pmd
        mc' = mc + 2*rdkc
        mrc = mr*mc*md
        reIndex idx =
            let (idx',idd) = divMod idx md
                (ir,ic) = divMod idx' mc
             in ((ir + kr) * mc' + (ic + kc)) * md + idd
     in G.Vector $ S.generate mrc reIndex

padMatrix'
    :: forall rdkr rdkc mr mc md x
    . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md, Num x, Storable x)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Proxy md
    -> Vector (mr*mc*md) x
    -> Vector ((mr + 2*rdkr)*(mc + 2*rdkc)*md) x
{-# INLINE padMatrix' #-}
padMatrix' _ _ _ _ _ v =
    let mtx :: Matrix mr mc (Vector md x)
        mtx = G.Matrix $ breakEvery v
        pdrs :: Vector rdkr (Vector mc (Vector md x))
        pdrs = replicate . replicate $ replicate 0
        mtx' = fromRows $ pdrs ++ toRows mtx ++ pdrs
        pdcs :: Vector rdkc (Vector (mr + 2*rdkr) (Vector md x))
        pdcs = replicate . replicate $ replicate 0
     in concat . G.toVector . G.fromColumns $ pdcs ++ G.toColumns mtx' ++ pdcs

im2colIndices'
    :: (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Proxy md
    -> Vector ((2*rdkr+1)*(2*rdkc+1)*mr*mc*md) Int
{-# INLINE im2colIndices' #-}
im2colIndices' prdkr prdkc pmr pmc pmd =
    let rdkr = natValInt prdkr
        rdkc = natValInt prdkc
        dms = ((2*rdkr+1)*(2*rdkc+1))
        kc = (2*rdkc + 1)
        reWindow idx =
            let (ir,ic) = divMod idx kc
             in kernelIndices' prdkc pmr pmc pmd ir ic
          in concatMap reWindow . G.Vector $ S.generate dms id

im2col'
    :: forall rdkr rdkc mr mc md x
    . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md, Num x, Storable x)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Proxy md
    -> Vector (mr*mc*md) x
    -> Matrix ((2*rdkr+1)*(2*rdkc+1)*md) (mr*mc) x
{-# INLINE im2col' #-}
im2col' prdkr prdkc pmr pmc pmd mtx =
    let idxs = im2colIndices' prdkr prdkc pmr pmc pmd
        mtx' = padMatrix' prdkr prdkc pmr pmc pmd mtx
     in G.Matrix $ backpermute mtx' idxs

crossCorrelate2d'
    :: forall nk rdkr rdkc mr mc md x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md
      , KnownNat nk, Numeric x, Storable x )
      => Proxy rdkr
      -> Proxy rdkc
      -> Proxy mr
      -> Proxy mc
      -> Proxy md
      -> Matrix nk ((2*rdkr+1)*(2*rdkc+1)*md) x
      -> Vector (mr*mc*md) x
      -> Matrix nk (mr*mc) x
{-# INLINE crossCorrelate2d' #-}
crossCorrelate2d' prdkr prdkc pmr pmc pmd krns v =
    let mtx = im2col' prdkr prdkc pmr pmc pmd v
     in matrixMatrixMultiply krns mtx


kernelIndices
    :: forall rdkc mr mc . (KnownNat rdkc, KnownNat mr, KnownNat mc)
    => Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Int
    -> Int
    -> Vector (mr*mc) Int
{-# INLINE kernelIndices #-}
kernelIndices prdkc pmc pmr kr kc =
    let rdkc = natValInt prdkc
        mr = natValInt pmr
        mc = natValInt pmc
        mc' = mc + 2*rdkc
        mrc = mr*mc
        reIndex idx =
            let (ir,ic) = divMod idx mc
             in (ir + kr) * mc' + (ic + kc)
     in G.Vector $ S.generate mrc reIndex

padMatrix
    :: forall rdkr rdkc mr mc x . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, Num x, Storable x)
    => Proxy rdkr
    -> Proxy rdkc
    -> Matrix mr mc x
    -> Matrix (mr + 2*rdkr) (mc + 2*rdkc) x
{-# INLINE padMatrix #-}
padMatrix _ _ mtx =
    let pdrs :: Vector rdkr (Vector mc x)
        pdrs = replicate (replicate 0)
        mtx' = fromRows $ pdrs ++ toRows mtx ++ pdrs
        pdcs :: Vector rdkc (Vector (mr + 2*rdkr) x)
        pdcs = replicate (replicate 0)
     in G.fromColumns $ pdcs ++ G.toColumns mtx' ++ pdcs

im2colIndices
    :: (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Vector ((2*rdkr+1)*(2*rdkc+1)*mr*mc) Int
{-# INLINE im2colIndices #-}
im2colIndices prdkc prdkr pmr pmc =
    let rdkc = natValInt prdkc
        rdkr = natValInt prdkr
        dms = ((2*rdkr+1)*(2*rdkc+1))
        kc = (2*rdkc + 1)
        reWindow idx =
            let (ir,ic) = divMod idx kc
             in kernelIndices prdkc pmr pmc ir ic
          in concatMap reWindow . G.Vector $ S.generate dms id

im2col
    :: forall rdkr rdkc mr mc x . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, Num x, Storable x)
    => Vector ((2*rdkr+1)*(2*rdkc+1)*mr*mc) Int
    -> Matrix mr mc x
    -> Matrix ((2*rdkr+1)*(2*rdkc+1)) (mr*mc) x
{-# INLINE im2col #-}
im2col idxs mtx =
    let mtx' = padMatrix (Proxy :: Proxy rdkc) (Proxy :: Proxy rdkr) mtx
     in G.Matrix $ backpermute (G.toVector mtx') idxs

-- | The dot product of two numerical 'Vector's.
addMatrix :: Numeric x => Matrix m n x -> Matrix m n x -> Matrix m n x
{-# INLINE addMatrix #-}
addMatrix mtx1 mtx2 = G.Matrix $ add (G.toVector mtx1) (G.toVector mtx2)

matrixMap :: (Storable x, Storable y, KnownNat m, KnownNat n)
          => (x -> y) -> Matrix m n x -> Matrix m n y
{-# INLINE matrixMap #-}
matrixMap f (G.Matrix v) = G.Matrix $ map f v

crossCorrelate2d0
    :: forall nk rdkr rdkc mr mc x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat nk, Numeric x, Storable x )
      => Vector ((2*rdkr+1)*(2*rdkc+1)*mr*mc) Int
      -> Vector nk (Matrix (2*rdkr+1) (2*rdkc+1) x)
      -> Matrix mr mc x
      -> Vector nk (Matrix mr mc x)
{-# INLINE crossCorrelate2d0 #-}
crossCorrelate2d0 idxs krns mtx =
    let mtx' = im2col idxs mtx
        krn' = fromRows . map G.toVector $ krns
     in map G.Matrix . toRows $ matrixMatrixMultiply krn' mtx'

crossCorrelate2d
    :: forall nk rdkr rdkc mr mc d x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat nk, KnownNat d, Numeric x, Storable x )
      => Vector d (Vector nk (Matrix (2*rdkr+1) (2*rdkc+1) x))
      -> Vector d (Matrix mr mc x)
      -> Vector nk (Matrix mr mc x)
{-# INLINE crossCorrelate2d #-}
crossCorrelate2d =
    let prdkr = Proxy :: Proxy rdkr
        prdkc = Proxy :: Proxy rdkc
        pmr = Proxy :: Proxy mr
        pmc = Proxy ::Proxy mc
        idxs = im2colIndices prdkr prdkc pmr pmc
     in zipFold (\zs krns mtx -> zipWith addMatrix zs $ crossCorrelate2d0 idxs krns mtx) (replicate . G.Matrix $ replicate 0)

rotateKernel :: (KnownNat m, KnownNat n, Storable x) => Matrix m n x -> Matrix m n x
{-# INLINE rotateKernel #-}
rotateKernel = fromRows . reverse . map reverse . toRows

convolve2d
    :: forall nk rdkr rdkc mr mc d x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat d, KnownNat nk, Numeric x, Storable x )
      => Vector d (Vector nk (Matrix (2*rdkr+1) (2*rdkc+1) x))
      -> Vector nk (Matrix mr mc x)
      -> Vector d (Matrix mr mc x)
{-# INLINE convolve2d #-}
convolve2d krns0 mtxs =
    let krns = toRows . matrixMap rotateKernel . G.transpose $ fromRows krns0
     in crossCorrelate2d krns mtxs

kernelDifferential
    :: forall nk rdkr rdkc mr mc d x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat d, KnownNat nk, Numeric x, Storable x )
      => Vector d (Vector nk (Matrix (2*rdkr+1) (2*rdkc+1) x))
      -> Vector nk (Matrix mr mc x)
      -> Vector d (Matrix mr mc x)
      -> Vector d (Vector nk (Matrix (2*rdkr+1) (2*rdkc+1) x))
{-# INLINE kernelDifferential #-}
kernelDifferential = undefined

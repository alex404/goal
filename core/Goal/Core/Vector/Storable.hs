 {-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
-- | Vectors and Matrices with statically typed dimensions based on storable vectors and using HMatrix where possible.
module Goal.Core.Vector.Storable
    ( -- * Vector
      module Data.Vector.Storable.Sized
      -- ** Construction
    , doubleton
    , range
      -- ** Deconstruction
    , concat
    , breakEvery
    , toPair
    -- ** Computation
    , average
    , zipFold
    -- * Matrix
    , Matrix
    , nRows
    , nColumns
    -- ** Construction
    , fromRows
    , fromColumns
    , matrixIdentity
    , outerProduct
    , sumOuterProduct
    , averageOuterProduct
    , weightedAverageOuterProduct
    , diagonalMatrix
    , fromLowerTriangular
    -- ** Deconstruction
    , toRows
    , toColumns
    , lowerTriangular
    , takeDiagonal
    -- ** Manipulation
    , columnVector
    , rowVector
    , combineTriangles
    , diagonalConcat
    -- ** Computation
    , trace
    , withMatrix
    -- *** BLAS
    , scale
    , add
    , dotProduct
    , dotMap
    , matrixVectorMultiply
    , matrixMatrixMultiply
    , matrixMap
    , eigens
    , isSemiPositiveDefinite
    , determinant
    , inverseLogDeterminant
    , inverse
    , pseudoInverse
    , matrixRoot
    , transpose
    -- *** Least Squares
    , linearLeastSquares
    , meanSquaredError
    , rSquared
    , l2Norm
    , unsafeCholesky
    -- *** Convolutions
    , crossCorrelate2d
    , convolve2d
    , kernelOuterProduct
    , kernelTranspose
    -- ** Miscellaneous
    , prettyPrintMatrix
    ) where


--- Imports ---


-- Goal --

import Goal.Core.Util hiding (average,breakEvery,range)

-- Unqualified --

import Data.Proxy
import Data.Complex
import Foreign.Storable
import Data.Vector.Storable.Sized
import Numeric.LinearAlgebra (Field,Numeric)
import GHC.TypeNats
import Prelude hiding (concat,foldr1,concatMap,replicate,(++),length,map,sum,zip,and)

-- Qualified --

import qualified Prelude
import qualified Data.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G
import qualified Data.Vector.Generic.Sized.Internal as G
import qualified Numeric.LinearAlgebra as H
import qualified Data.List as L


--- Generic ---


-- | Matrices with static dimensions (storable).
type Matrix = G.Matrix S.Vector

-- | A fold over pairs of elements of 'Vector's of equal length.
zipFold :: (KnownNat n, Storable x, Storable y) => (z -> x -> y -> z) -> z -> Vector n x -> Vector n y -> z
{-# INLINE zipFold #-}
zipFold f z0 xs ys =
    let n = length xs
        foldfun z i = f z (unsafeIndex xs i) (unsafeIndex ys i)
     in L.foldl' foldfun z0 [0..n-1]

-- | Concatenates a vector of vectors into a single vector.
concat :: (KnownNat n, Storable x) => Vector m (Vector n x) -> Vector (m*n) x
{-# INLINE concat #-}
concat = G.concat

-- | Collect two values into a length 2 'Vector'.
doubleton :: Storable x => x -> x -> Vector 2 x
{-# INLINE doubleton #-}
doubleton = G.doubleton

-- | The number of rows in the 'Matrix'.
nRows :: forall m n a . KnownNat m => Matrix m n a -> Int
{-# INLINE nRows #-}
nRows = G.nRows

-- | The columns of columns in the 'Matrix'.
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

-- | Uniform partition of an interval into a 'Vector'.
range :: (KnownNat n, Fractional x, Storable x) => x -> x -> Vector n x
{-# INLINE range #-}
range = G.range

-- | Reshapes a length 2 'Vector' into a pair of values.
toPair :: Storable x => Vector 2 x -> (x,x)
{-# INLINE toPair #-}
toPair = G.toPair


--- HMatrix ---


-- | Converts a pure, Storable-based 'Matrix' into an HMatrix matrix.
toHMatrix :: forall m n x . (KnownNat n, KnownNat m, H.Element x, Storable x) => Matrix m n x -> H.Matrix x
{-# INLINE toHMatrix #-}
toHMatrix (G.Matrix mtx) =
    let n = natValInt (Proxy :: Proxy n)
        m = natValInt (Proxy :: Proxy m)
     in if n == 0
           then H.fromRows $ Prelude.replicate m S.empty
           else H.reshape n $ fromSized mtx

-- | Converts an HMatrix matrix into a pure, Storable-based 'Matrix'.
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

-- | Breaks a 'Vector' into a Vector of Vectors.
breakEvery :: forall n k a . (KnownNat n, KnownNat k, Storable a) => Vector (n*k) a -> Vector n (Vector k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> G.Vector $ S.unsafeSlice (finiteInt i*k) k v)


--- BLAS ---


-- | The sum of two 'Vector's.
add :: Numeric x => Vector n x -> Vector n x -> Vector n x
{-# INLINE add #-}
add (G.Vector v1) (G.Vector v2) = G.Vector (H.add v1 v2)

-- | Scalar multiplication of a 'Vector'.
scale :: Numeric x => x -> Vector n x -> Vector n x
{-# INLINE scale #-}
scale x (G.Vector v) = G.Vector (H.scale x v)

-- | Apply a 'Vector' operation to a 'Matrix'.
withMatrix :: (Vector (n*m) x -> Vector (n*m) x) -> Matrix n m x -> Matrix n m x
{-# INLINE withMatrix #-}
withMatrix f (G.Matrix v) = G.Matrix $ f v

-- | Returns the lower triangular part of a square matrix.
lowerTriangular :: forall n x . (Storable x, H.Element x, KnownNat n) => Matrix n n x -> Vector (Triangular n) x
{-# INLINE lowerTriangular #-}
lowerTriangular mtx =
    let hmtx = toHMatrix mtx
        rws = H.toRows hmtx
        rws' = Prelude.zipWith S.take [1..] rws
     in G.Vector $ S.concat rws'
--    let n = natValInt (Proxy :: Proxy n)
--        idxs = G.Vector . S.fromList
--            $ Prelude.concat [ from2Index n <$> Prelude.zip (repeat k) [0..k] | k <- [0..n-1] ]
--     in backpermute xs idxs

toTriangularIndex :: (Int,Int) -> Int
{-# INLINE toTriangularIndex #-}
toTriangularIndex (i,j)
    | i >= j = triangularNumber i + j
    | otherwise = toTriangularIndex (j,i)

-- | Constructs a `Matrix` from a lower triangular part.
fromLowerTriangular :: forall n x . (Storable x, KnownNat n) => Vector (Triangular n) x -> Matrix n n x
{-# INLINE fromLowerTriangular #-}
fromLowerTriangular xs =
    let n = natValInt (Proxy :: Proxy n)
        idxs = generate (toTriangularIndex . to2Index n . finiteInt)
     in G.Matrix $ backpermute xs idxs

-- | Build a matrix with the given diagonal, lower triangular part given by the
-- first matrix, and upper triangular part given by the second matrix.
combineTriangles
    :: (KnownNat k, Storable x)
    => Vector k x -- ^ Diagonal
    -> Matrix k k x -- ^ Lower triangular source
    -> Matrix k k x -- ^ Upper triangular source
    -> Matrix k k x
{-# INLINE combineTriangles #-}
combineTriangles (G.Vector diag) crs1 crs2 =
    fromRows $ generate (generator (toRows crs1) (toRows crs2))
        where
            generator rws1 rws2 fnt =
                let (G.Vector rw1) = index rws1 fnt
                    (G.Vector rw2) = index rws2 fnt
                    i = fromIntegral fnt
                 in G.Vector $ S.take i rw1 S.++ S.cons (diag S.! i) (S.drop (i+1) rw2)




to2Index :: Int -> Int -> (Int,Int)
{-# INLINE to2Index #-}
to2Index nj ij = divMod ij nj

--from2Index :: Int -> (Int,Int) -> Int
--{-# INLINE from2Index #-}
--from2Index nj (i,j) = i*nj + j

-- | The average of a 'Vector' of elements.
average :: (Numeric x, Fractional x) => Vector n x -> x
{-# INLINE average #-}
average (G.Vector v) = H.sumElements v / fromIntegral (S.length v)

-- | The dot product of two numerical 'Vector's.
dotProduct :: Numeric x => Vector n x -> Vector n x -> x
{-# INLINE dotProduct #-}
dotProduct v1 v2 = H.dot (fromSized v1) (fromSized v2)

-- | The determinant of a 'Matrix'.
diagonalMatrix :: forall n x . (KnownNat n, Field x) => Vector n x -> Matrix n n x
{-# INLINE diagonalMatrix #-}
diagonalMatrix v =
    let n = natValInt (Proxy :: Proxy n)
     in fromHMatrix $ H.diagRect 0 (fromSized v) n n

-- | The determinant of a 'Matrix'.
takeDiagonal :: (KnownNat n, Field x) => Matrix n n x -> Vector n x
{-# INLINE takeDiagonal #-}
takeDiagonal = G.Vector . H.takeDiag . toHMatrix

-- | The determinant of a 'Matrix'.
trace :: (KnownNat n, Field x) => Matrix n n x -> x
{-# INLINE trace #-}
trace = S.sum . H.takeDiag . toHMatrix

-- | Returns the eigenvalues and eigenvectors 'Matrix'.
eigens :: (KnownNat n, Field x) => Matrix n n x -> (Vector n (Complex Double), Vector n (Vector n (Complex Double)))
{-# INLINE eigens #-}
eigens mtx =
    let (exs,evs) = H.eig $ toHMatrix mtx
     in (G.Vector exs, G.Vector . S.fromList $ G.Vector <$> H.toColumns evs)

-- | Test if the matrix is semi-positive definite.
isSemiPositiveDefinite :: (KnownNat n, Field x) => Matrix n n x -> Bool
{-# INLINE isSemiPositiveDefinite #-}
isSemiPositiveDefinite =
    and . map ((0 <=) . realPart) . fst . eigens

-- | Returns the inverse, the logarithm of the absolute value of the
-- determinant, and the sign of the determinant of a given matrix.
inverseLogDeterminant :: (KnownNat n, Field x) => Matrix n n x -> (Matrix n n x, x, x)
{-# INLINE inverseLogDeterminant #-}
inverseLogDeterminant mtx =
    let (imtx,(ldet,sgn)) = H.invlndet $ toHMatrix mtx
     in (fromHMatrix imtx, ldet, sgn)

-- | The determinant of a 'Matrix'.
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

-- | Invert a 'Matrix'.
pseudoInverse :: forall n x . (KnownNat n, Field x) => Matrix n n x -> Matrix n n x
{-# INLINE pseudoInverse #-}
pseudoInverse (G.Matrix mtx) =
    G.Matrix $ withVectorUnsafe (H.flatten . H.pinv . H.reshape (natValInt (Proxy :: Proxy n))) mtx

-- | Square root of a 'Matrix'.
matrixRoot :: forall n x . (KnownNat n, Field x) => Matrix n n x -> Matrix n n x
{-# INLINE matrixRoot #-}
matrixRoot (G.Matrix mtx) =
    G.Matrix $ withVectorUnsafe (H.flatten . H.sqrtm . H.reshape (natValInt (Proxy :: Proxy n))) mtx

-- | The outer product of two 'Vector's.
outerProduct :: (KnownNat m, KnownNat n, Numeric x) => Vector m x -> Vector n x -> Matrix m n x
{-# INLINE outerProduct #-}
outerProduct v1 v2 =
    fromHMatrix $ H.outer (fromSized v1) (fromSized v2)

-- | The average outer product of two lists of 'Vector's.
sumOuterProduct :: (KnownNat m, KnownNat n, Fractional x, Numeric x) => [(Vector m x,Vector n x)] -> Matrix m n x
{-# INLINE sumOuterProduct #-}
sumOuterProduct v12s =
    let (v1s,v2s) = L.unzip v12s
        mtx1 = H.fromColumns $ fromSized <$> v1s
        mtx2 = H.fromRows $ fromSized <$> v2s
     in fromHMatrix (mtx1 H.<> mtx2)

-- | The average outer product of two lists of 'Vector's.
averageOuterProduct :: (KnownNat m, KnownNat n, Fractional x, Numeric x) => [(Vector m x,Vector n x)] -> Matrix m n x
{-# INLINE averageOuterProduct #-}
averageOuterProduct v12s =
    let (v1s,v2s) = L.unzip v12s
        mtx1 = H.fromColumns $ fromSized <$> v1s
        (_,n) = H.size mtx1
        mtx2 = H.scale (1/fromIntegral n) . H.fromRows $ fromSized <$> v2s
     in fromHMatrix (mtx1 H.<> mtx2)

-- | The average outer product of two lists of 'Vector's.
weightedAverageOuterProduct
    :: ( KnownNat m, KnownNat n, Fractional x, Numeric x )
    => [(x,Vector m x,Vector n x)]
    -> Matrix m n x
{-# INLINE weightedAverageOuterProduct #-}
weightedAverageOuterProduct wv12s =
    let (ws,v1s,v2s) = L.unzip3 wv12s
        v1s' = L.zipWith H.scale ws $ fromSized <$> v1s
        mtx1 = H.fromColumns v1s'
        mtx2 = H.fromRows $ fromSized <$> v2s
     in fromHMatrix (mtx1 H.<> mtx2)

-- | The identity 'Matrix'.
matrixIdentity :: forall n x . (KnownNat n, Numeric x, Num x) => Matrix n n x
{-# INLINE matrixIdentity #-}
matrixIdentity =
    fromHMatrix . H.ident $ natValInt (Proxy :: Proxy n)

-- | Apply a linear transformation to a 'Vector'.
dotMap :: (KnownNat n, Numeric x) => Vector n x -> [Vector n x] -> [x]
{-# INLINE dotMap #-}
dotMap v vs =
    let mtx' = H.fromRows $ fromSized <$> vs
     in H.toList $ mtx' H.#> fromSized v
--     in if S.null w
--           then replicate 0
--           else fmap G.Vector . H.toColumns $ toHMatrix mtx H.<> mtx'

-- | Apply a linear transformation to a 'Vector'.
matrixMap :: (KnownNat m, KnownNat n, Numeric x)
                     => Matrix m n x -> [Vector n x] -> [Vector m x]
{-# INLINE matrixMap #-}
matrixMap mtx vs =
    let mtx' = H.fromColumns $ fromSized <$> vs
     in fmap G.Vector . H.toColumns $ toHMatrix mtx H.<> mtx'
--     in if S.null w
--           then replicate 0
--           else fmap G.Vector . H.toColumns $ toHMatrix mtx H.<> mtx'


-- | Map a linear transformation over a list of 'Vector's.
matrixVectorMultiply :: (KnownNat m, KnownNat n, Numeric x)
                     => Matrix m n x -> Vector n x -> Vector m x
{-# INLINE matrixVectorMultiply #-}
matrixVectorMultiply mtx v =
    G.Vector $ toHMatrix mtx H.#> fromSized v
--    let w = toHMatrix mtx H.#> fromSized v
--     in if S.null w
--           then replicate 0
--           else G.Vector w

-- | Multiply a 'Matrix' with a second 'Matrix'.
matrixMatrixMultiply
    :: (KnownNat m, KnownNat n, KnownNat o, Numeric x)
    => Matrix m n x
    -> Matrix n o x
    -> Matrix m o x
{-# INLINE matrixMatrixMultiply #-}
matrixMatrixMultiply mtx1 mtx2 = fromHMatrix $ toHMatrix mtx1 H.<> toHMatrix mtx2

-- | Pretty print the values of a 'Matrix' (for extremely simple values of pretty).
prettyPrintMatrix :: (KnownNat m, KnownNat n, Numeric a, Show a) => Matrix m n a -> IO ()
prettyPrintMatrix = print . toHMatrix

-- | The Mean Squared difference between two vectors.
meanSquaredError
    :: KnownNat k
    => Vector k Double
    -> Vector k Double
    -> Double
{-# INLINE meanSquaredError #-}
meanSquaredError ys yhts = average $ map square (ys - yhts)

-- | L2 length of a vector.
l2Norm
    :: KnownNat k
    => Vector k Double
    -> Double
{-# INLINE l2Norm #-}
l2Norm (G.Vector xs) = H.norm_2 xs

-- | Computes the coefficient of determintation for the given outputs and model
-- predictions.
rSquared
    :: KnownNat k
    => Vector k Double -- ^ Dependent variable observations
    -> Vector k Double -- ^ Predicted Values
    -> Double -- ^ R-squared
{-# INLINE rSquared #-}
rSquared ys yhts =
    let ybr = average ys
        ssres = sum $ map square (ys - yhts)
        sstot = sum $ map (square . subtract ybr) ys
     in 1 - (ssres/sstot)

-- | Solves the linear least squares problem.
linearLeastSquares
    :: KnownNat l
    => [Vector l Double] -- ^ Independent variable observations
    -> [Double] -- ^ Dependent variable observations
    -> Vector l Double -- ^ Parameter estimates
{-# INLINE linearLeastSquares #-}
linearLeastSquares as xs =
    G.Vector $ H.fromRows (fromSized <$> as) H.<\> S.fromList xs

unsafeCholesky
    :: (KnownNat n, Field x, Storable x)
    => Matrix n n x
    -> Matrix n n x
unsafeCholesky =
    transpose . fromHMatrix . H.chol . H.trustSym . toHMatrix


--- Convolutions ---


to3Index :: Int -> Int -> Int -> (Int,Int,Int)
{-# INLINE to3Index #-}
to3Index nj nk ijk =
    let nj' = nj*nk
        (i,jk) = divMod ijk nj'
        (j,k) = divMod jk nk
     in (i,j,k)

from3Index :: Int -> Int -> (Int,Int,Int) -> Int
{-# INLINE from3Index #-}
from3Index nj nk (i,j,k) =
    let nj' = nj*nk
     in i*nj' + j*nk + k


windowIndices
    :: forall rdkr rdkc mr mc . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy mr
    -> Proxy mc
    -> Int
    -> Int
    -> Int
    -> Vector (mr*mc) Int
{-# INLINE windowIndices #-}
windowIndices prdkr prdkc pmr pmc kd kr kc =
    let rdkr = natValInt prdkr
        rdkc = natValInt prdkc
        mr = natValInt pmr
        mc = natValInt pmc
        mrc = mr*mc
        nj' = mr + 2*rdkr
        nk' = mc + 2*rdkc
        reIndex idx =
            let (j,k) = divMod idx mc
             in from3Index nj' nk' (kd,j+kr,k+kc)
     in G.Vector $ S.generate mrc reIndex

padMatrix
    :: forall rdkr rdkc mr mc md x
    . (KnownNat rdkr, KnownNat rdkc, KnownNat md, KnownNat mr, KnownNat mc, Num x, Storable x)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy md
    -> Proxy mr
    -> Proxy mc
    -> Vector (md*mr*mc) x
    -> Vector (md*(mr + 2*rdkr)*(mc + 2*rdkc)) x
{-# INLINE padMatrix #-}
padMatrix _ _ _ _ _ v =
    let mtxs :: Vector md (Matrix mr mc x)
        mtxs = map G.Matrix $ breakEvery v
        pdrs :: Vector rdkr (Vector mc x)
        pdrs = replicate $ replicate 0
        mtxs' = map (\mtx -> fromRows $ pdrs ++ toRows mtx ++ pdrs) mtxs
        pdcs :: Vector rdkc (Vector (mr + 2*rdkr) x)
        pdcs = replicate $ replicate 0
     in concatMap G.toVector $ map (\mtx' -> G.fromColumns $ pdcs ++ G.toColumns mtx' ++ pdcs) mtxs'

im2colIndices
    :: forall rdkr rdkc mr mc md
     . (KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc, KnownNat md)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy md
    -> Proxy mr
    -> Proxy mc
    -> Vector (((2*rdkr+1)*(2*rdkc+1)*md)*(mr*mc)) Int
{-# INLINE im2colIndices #-}
im2colIndices prdkr prdkc _ pmr pmc =
    let rdkr = natValInt prdkr
        rdkc = natValInt prdkc
        nj = (2*rdkr + 1)
        nk = (2*rdkc + 1)
        reWindow idx =
            let (i,j,k) = to3Index nj nk idx
             in windowIndices prdkr prdkc pmr pmc i j k
          in (concatMap reWindow :: Vector ((2*rdkr+1)*(2*rdkc+1)*md) Int -> Vector (((2*rdkr+1)*(2*rdkc+1)*md)*(mr*mc)) Int) $ generate finiteInt

im2col
    :: forall rdkr rdkc md mr mc x
    . (KnownNat rdkr, KnownNat rdkc, KnownNat mc, KnownNat md, KnownNat mr, Num x, Storable x)
    => Proxy rdkr
    -> Proxy rdkc
    -> Proxy md
    -> Proxy mr
    -> Proxy mc
    -> Vector (md*mr*mc) x
    -> Matrix (md*(2*rdkr+1)*(2*rdkc+1)) (mr*mc) x
{-# INLINE im2col #-}
im2col prdkr prdkc pmd pmr pmc mtx =
    let idxs = im2colIndices prdkr prdkc pmd pmr pmc
        mtx' = padMatrix prdkr prdkc pmd pmr pmc mtx
     in G.Matrix $ backpermute mtx' idxs

-- | 2d cross-correlation of a kernel over a matrix of values.
crossCorrelate2d
    :: forall nk rdkr rdkc mr mc md x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat md, KnownNat mr, KnownNat mc
      , KnownNat nk, Numeric x, Storable x )
      => Proxy rdkr -- ^ Number of Kernel rows
      -> Proxy rdkc -- ^ Number of Kernel columns
      -> Proxy mr -- ^ Number of Matrix/Image rows
      -> Proxy mc -- ^ Number of Kernel/Image columns
      -> Matrix nk (md*(2*rdkr+1)*(2*rdkc+1)) x -- ^ Kernels (nk is their number)
      -> Matrix md (mr*mc) x -- ^ Image (md is the depth)
      -> Matrix nk (mr*mc) x -- ^ Cross-correlated image
{-# INLINE crossCorrelate2d #-}
crossCorrelate2d prdkr prdkc pmr pmc krns (G.Matrix v) =
    let pmd = Proxy :: Proxy md
        mtx = im2col prdkr prdkc pmd pmr pmc v
     in matrixMatrixMultiply krns mtx

to4Index :: Int -> Int -> Int -> Int -> (Int,Int,Int,Int)
{-# INLINE to4Index #-}
to4Index nj nk nl ijkl =
    let nk' = nl*nk
        nj' = nj*nk'
        (i,jkl) = divMod ijkl nj'
        (j,kl) = divMod jkl nk'
        (k,l) = divMod kl nl
     in (i,j,k,l)

from4Index :: Int -> Int -> Int -> (Int,Int,Int,Int) -> Int
{-# INLINE from4Index #-}
from4Index nj nk nl (i,j,k,l) =
    let nk' = nl*nk
        nj' = nj*nk'
     in i*nj' + j*nk' + k*nl + l

kernelTransposeIndices
    :: (KnownNat nk, KnownNat md, KnownNat rdkr, KnownNat rdkc)
    => Proxy nk
    -> Proxy md
    -> Proxy rdkr
    -> Proxy rdkc
    -> Vector (nk*md*(2*rdkr+1)*(2*rdkc+1)) Int
{-# INLINE kernelTransposeIndices #-}
kernelTransposeIndices pnk pmd prdkr prdkc =
    let nkrn = natValInt pnk
        md = natValInt pmd
        rdkr = natValInt prdkr
        rdkc = natValInt prdkc
        dmkr = 2*rdkr+1
        dmkc = 2*rdkc+1
        nl = dmkc
        nk = dmkr
        nj = nkrn
        nl' = dmkc
        nk' = dmkr
        nj' = md
        reIndex idx =
            let (i,j,k,l) = to4Index nj nk nl idx
             in from4Index nj' nk' nl' (j,i,nk-1-k,nl-1-l)
     in generate (reIndex . fromIntegral)

-- | The transpose of a convolutional kernel.
kernelTranspose
    :: (KnownNat nk, KnownNat md, KnownNat rdkr, KnownNat rdkc, Numeric x, Storable x)
    => Proxy nk
    -> Proxy md
    -> Proxy rdkr
    -> Proxy rdkc
    -> Matrix nk (md*(2*rdkr+1)*(2*rdkc+1)) x -- ^ Kernels (nk is their number)
    -> Matrix md (nk*(2*rdkr+1)*(2*rdkc+1)) x -- ^ Kernels (nk is their number)
{-# INLINE kernelTranspose #-}
kernelTranspose pnk pmd prdkr prdkc (G.Matrix kv) = G.Matrix . backpermute kv $ kernelTransposeIndices pnk pmd prdkr prdkc

-- | 2d convolution of a kernel over a matrix of values. This is the adjoint of crossCorrelate2d.
convolve2d
    :: forall nk rdkr rdkc md mr mc x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat md, KnownNat nk, Numeric x, Storable x )
      => Proxy rdkr -- ^ Number of Kernel rows
      -> Proxy rdkc -- ^ Number of Kernel columns
      -> Proxy mr -- ^ Number of Matrix/Image rows
      -> Proxy mc -- ^ Number of Kernel/Image columns
      -> Matrix nk (md*(2*rdkr+1)*(2*rdkc+1)) x -- ^ Kernels (nk is their number)
      -> Matrix nk (mr*mc) x -- ^ Dual image (nk is its depth)
      -> Matrix md (mr*mc) x -- ^ Convolved image
{-# INLINE convolve2d #-}
convolve2d prdkr prdkc pmr pmc krn mtxs =
    let pnk = Proxy :: Proxy nk
        pmd = Proxy :: Proxy md
        krn' = kernelTranspose pnk pmd prdkr prdkc krn
     in crossCorrelate2d prdkr prdkc pmr pmc krn' mtxs

-- | The outer product of an image and a dual image to produce a convolutional kernel.
kernelOuterProduct
    :: forall nk rdkr rdkc md mr mc x
    . ( KnownNat rdkr, KnownNat rdkc, KnownNat mr, KnownNat mc
      , KnownNat md, KnownNat nk, Numeric x, Storable x )
      => Proxy rdkr -- ^ Number of Kernel rows
      -> Proxy rdkc -- ^ Number of Kernel columns
      -> Proxy mr -- ^ Number of Matrix/Image rows
      -> Proxy mc -- ^ Number of Kernel/Image columns
      -> Matrix nk (mr*mc) x -- ^ Dual image (nk is its depth)
      -> Matrix md (mr*mc) x -- ^ Image (md is the depth)
      -> Matrix nk (md*(2*rdkr+1)*(2*rdkc+1)) x -- ^ Kernels
{-# INLINE kernelOuterProduct #-}
kernelOuterProduct prdkr prdkc pmr pmc omtx (G.Matrix v) =
    let pmd = Proxy :: Proxy md
        imtx = im2col prdkr prdkc pmd pmr pmc v
     in matrixMatrixMultiply omtx $ transpose imtx

{-# LANGUAGE
   DataKinds,
   FlexibleInstances,
   FlexibleContexts,
   KindSignatures,
   ConstraintKinds,
   TypeOperators,
   TypeApplications,
   ScopedTypeVariables,
   RankNTypes,
   StandaloneDeriving,
   GeneralizedNewtypeDeriving,
   GADTs
   #-}

-- | Vectors and Matrices with statically typed dimensions.
module Goal.Core.Vector.Generic
    ( -- * Vector
      module Data.Vector.Generic.Sized
    , VectorClass
    , concat
    , doubleton
    , breakEvery
    , range
    , generateP
    , generatePM
    -- * Matrix
    , Matrix (Matrix,toVector)
    -- ** Construction
    , fromRows
    , fromColumns
    -- ** Deconstruction
    , toPair
    , toRows
    , toColumns
    , nRows
    , nColumns
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

import Data.Functor.Identity
import GHC.TypeNats
import Data.Proxy
import Control.DeepSeq
import Data.Vector.Generic.Sized
import Data.Vector.Generic.Sized.Internal
import Foreign.Storable
import Data.Type.Equality
import Numeric.Natural
import Prelude hiding (concatMap,concat,map,sum)

-- Qualified --

import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S


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

-- | Matrices with static dimensions.
newtype Matrix v (m :: Nat) (n :: Nat) a = Matrix { toVector :: Vector v (m*n) a }
    deriving (Eq,Show,NFData)

deriving instance (KnownNat m, KnownNat n, Storable x) => Storable (Matrix S.Vector m n x)

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

-- | Breaks a 'Vector' into a Vector of Vectors.
breakEvery
    :: forall v n k a . (G.Vector v a, G.Vector v (Vector v k a), KnownNat n, KnownNat k)
    => Vector v (n*k) a -> Vector v n (Vector v k a)
{-# INLINE breakEvery #-}
breakEvery v0 =
    let k = natValInt (Proxy :: Proxy k)
        v = fromSized v0
     in generate (\i -> Vector $ G.unsafeSlice (finiteInt i*k) k v)

-- | The number of rows in the 'Matrix'.
nRows :: forall v m n a . KnownNat m => Matrix v m n a -> Int
{-# INLINE nRows #-}
nRows _ = natValInt (Proxy :: Proxy m)

-- | The number of columns in the 'Matrix'.
nColumns :: forall v m n a . KnownNat n => Matrix v m n a -> Int
{-# INLINE nColumns #-}
nColumns _ = natValInt (Proxy :: Proxy n)

-- | Reshapes a length 2 'Vector' into a pair of values.
toPair :: G.Vector v a => Vector v 2 a -> (a,a)
toPair v = (unsafeIndex v 0, unsafeIndex v 1)

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

-- | Uniform partition of an interval into a 'Vector'.
range
    :: forall v n x. (G.Vector v x, KnownNat n, Fractional x)
    => x -> x -> Vector v n x
{-# INLINE range #-}
range mn mx =
    let n = natValInt (Proxy :: Proxy n)
        stp = (mx - mn)/fromIntegral (n-1)
     in enumFromStepN mn stp


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


--- GenerateP ---


generateP0
    :: forall n m i x . (KnownNat n, KnownNat i, Monad m)
    => Proxy n
    -> Proxy i
    -> Natural
    -> (forall i' j . (KnownNat i', KnownNat j, (i' + j + 1) ~ n) => Proxy i' -> m x)
    -> m [x]
generateP0 prxn prxi 0 f =
    case sameNat (Proxy @ (i+1)) prxn of
        Just Refl -> (:[]) <$> f prxi
        Nothing -> error "misuse of generateP0 function"
generateP0 prxn prxi j f = case someNatVal j of
    SomeNat (_ :: Proxy j) -> case sameNat (Proxy @ (i+j+1)) prxn of
        Just Refl -> do
            xm <- f prxi
            (xm :) <$> generateP0 prxn (Proxy @ (i+1)) (j-1) f
        Nothing -> error "misuse of generateP0 function"
    --        case (same
    --generateP0 prxn kp f `S.snoc` f (Proxy :: Proxy k)

-- | Vector generation given based on Proxied Nats. Incorporates a size
-- constraint into the generating function to allow functions which are bounded
-- by the size of the vector.
generateP
    :: forall v n x . (G.Vector v x, KnownNat n)
    => (forall i j . (KnownNat i, KnownNat j, (i + j + 1) ~ n) => Proxy i -> x)
    -> Vector v n x
generateP f0 =
    let prxn = Proxy @ n
        f :: (KnownNat i, KnownNat j, (i + j + 1) ~ n) => Proxy i -> Identity x
        f = return . f0
     in Vector . G.fromList . runIdentity $ generateP0 prxn (Proxy @ 0) (natVal prxn - 1) f

-- | Vector generation given based on Proxied Nats (Monadic Version).
generatePM
    :: forall v n m x . (G.Vector v x, KnownNat n, Monad m)
    => (forall i j . (KnownNat i, KnownNat j, (i + j + 1) ~ n) => Proxy i -> m x)
    -> m (Vector v n x)
generatePM f =
    let prxn = Proxy @ n
     in Vector . G.fromList <$> generateP0 prxn (Proxy @ 0) (natVal prxn - 1) f

---- | Right now the evaluated values are 1..k, which is a bit unusual.
--generatePM0'
--    :: forall n k x m . (KnownNat n, k <= n, KnownNat k, Monad m, Storable x, NFData x)
--    => Proxy n
--    -> NatPeano k
--    -> (forall j . (KnownNat j, j <= n, 1 <= n) => Proxy j -> m x)
--    -> m (S.Vector x)
--{-# INLINE generatePM0' #-}
--generatePM0' _ PeanoZero _ = return S.empty
--generatePM0' prxn (PeanoSucc kp) f = do
--        x <- f (Proxy :: Proxy k)
--        deepseq x $ (`S.snoc` x) <$> generatePM0' prxn kp f
--
--generatePM'
--    :: forall n m x . (KnownNat n, Monad m, Storable x, NFData x)
--    => (forall j . (KnownNat j, j <= n, 1 <= n) => Proxy j -> m x)
--    -> m (Vector n x)
--generatePM' f = I.Vector <$> generatePM0' (Proxy :: Proxy n) (natSingleton :: NatPeano n) f
--
---- | Right now the evaluated values are 1..k, which is a bit unusual.
--generatePM0
--    :: forall n k x m . (KnownNat n, k <= n, KnownNat k, Monad m, Storable x)
--    => Proxy n
--    -> NatPeano k
--    -> (forall j . (KnownNat j, j <= n, 1 <= n) => Proxy j -> m x)
--    -> m (S.Vector x)
--{-# INLINE generatePM0 #-}
--generatePM0 _ PeanoZero _ = return S.empty
--generatePM0 prxn (PeanoSucc kp) f = do
--        x <- f (Proxy :: Proxy k)
--        (`S.snoc` x) <$> generatePM0 prxn kp f
--
--generatePM
--    :: forall n m x . (KnownNat n, Monad m, Storable x)
--    => (forall j . (KnownNat j, j <= n, 1 <= n) => Proxy j -> m x)
--    -> m (Vector n x)
--generatePM f = I.Vector <$> generatePM0 (Proxy :: Proxy n) (natSingleton :: NatPeano n) f
--
--
--

--- Least Squares ---

---- | Linear least squares estimation.
--linearLeastSquares
--    :: (KnownNat l, KnownNat k, 1 <= k)
--    => Vector k (Vector l Double) -- ^ Independent variable observations
--    -> Vector k Double -- ^ Dependent variable observations
--    -> Vector l Double -- ^ Parameter estimates
--{-# INLINE linearLeastSquares #-}
--linearLeastSquares xs ys =
--    let mtx = fromRows xs
--     in linearLeastSquares0 mtx ys
--
---- | Linear least squares estimation, where the design matrix is provided directly.
--linearLeastSquares0
--    :: (KnownNat l, KnownNat k)
--    => Matrix k l Double -- ^ Independent variable observations
--    -> Vector k Double -- ^ Dependent variable observations
--    -> Vector l Double -- ^ Parameter estimates
--{-# INLINE linearLeastSquares0 #-}
--linearLeastSquares0 mtx ys =
--    let tmtx = transpose mtx
--        prj = matrixMatrixMultiply (inverse $ matrixMatrixMultiply tmtx mtx) tmtx
--     in matrixVectorMultiply prj ys

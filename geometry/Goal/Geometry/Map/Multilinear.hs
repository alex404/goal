{-# LANGUAGE UndecidableInstances #-}
-- | This module provides tools for working with tensors, affine transformations, and general
-- multilinear objects.

module Goal.Geometry.Map.Multilinear
    ( -- * Bilinear Forms
    Bilinear ((<.<), (<$<))
    -- * Tensors
    , Product
    -- ** Matrix Construction
    , toMatrix
    , fromMatrix
    , (>.<)
    -- ** Computation
    , (<#>)
    , inverse
    , transpose
    , determinant
    -- * Affine Functions
    , Affine
    , type (<*)
    -- ** Reshape Affine functions
    , splitAffine
    , joinAffine
    ) where

--- Imports ---

-- Package --

import Goal.Core

import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Map

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G
--import qualified Goal.Core.Vector.Boxed as B


-- Bilinear Forms --

class (Dimension m ~ Dimension n) => Square v m n x where
-- | The inverse of a tensor.
    inverse
        :: Point v (Function c d) (Product m n) x
        -> Point v (Function d c) (Product n m) x
    determinant :: Point v (Function c d) (Product m n) x -> x


-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class Apply v c d f x => Bilinear v c d f x where
    (<.<)
        :: Point v (Dual d) (Codomain f) x
        -> Point v (Function c d) f x
        -> Point v (Dual c) (Domain f) x
    (<$<) :: KnownNat k
        => Vector v k (Point v (Dual d) (Codomain f) x)
        -> Point v (Function c d) f x
        -> Vector v k (Point v (Dual c) (Domain f) x)

-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Product m n

-- | Infix synonym for a 'Tensor'.
--type (m * n) = Product m n
--infixr 6 *

-- | Converts a point on a 'Product' manifold into a 'Matrix'.
toMatrix :: (Manifold m, Manifold n) => SPoint c (Product m n) x -> S.Matrix (Dimension m) (Dimension n) x
{-# INLINE toMatrix #-}
toMatrix (Point xs) = Matrix xs

-- | Converts a 'Matrix' into a 'Point' on a 'Product' 'Manifold'.
fromMatrix :: S.Matrix (Dimension m) (Dimension n) x -> SPoint c (Product m n) x
{-# INLINE fromMatrix #-}
fromMatrix (Matrix xs) = Point xs

-- | The transpose of a tensor.
transpose
    :: (Manifold m, Manifold n, Numeric x)
    => SPoint (Function c d) (Product m n) x
    -> SPoint (Function (Dual d) (Dual c)) (Product n m) x
{-# INLINE transpose #-}
transpose (Point xs) = fromMatrix . S.transpose $ Matrix xs

-- | Tensor x Tensor multiplication.
(<#>) :: (Manifold m, Manifold n, Manifold o, Numeric x)
      => SPoint (Function d e) (Product m n) x
      -> SPoint (Function c d) (Product n o) x
      -> SPoint (Function c e) (Product m o) x
{-# INLINE (<#>) #-}
(<#>) m1 m2 =
    fromMatrix $ S.matrixMatrixMultiply (toMatrix m1) (toMatrix m2)

-- | '>.<' denotes the outer product between two points. It provides a way of
-- constructing matrices of the 'Tensor' product space.
(>.<) :: (Manifold m, Manifold n, Numeric x)
      => SPoint d m x
      -> SPoint c n x
      -> SPoint (Function (Dual c) d) (Product m n) x
(>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs


--- Affine Functions ---

-- | An 'Affine' 'Manifold' represents linear transformations followed by a translation.
data Affine f

-- | Infix synonym for 'Affine'.
type (m <* n) = Affine (Product m n)
infixr 6 <*

-- | Split a 'Point' on an 'Affine' 'Manifold' into a 'Point' which represents the translation, and a tensor.
splitAffine
    :: (GVector v x, Map f)
    => Point v (Function c d) (Affine f) x
    -> (Point v d (Codomain f) x, Point v (Function c d) f x)
{-# INLINE splitAffine #-}
splitAffine (Point cppqs) =
    let (cps,cpqs) = G.splitAt cppqs
     in (Point cps, Point cpqs)

-- | Join a translation and a tensor into a 'Point' on an 'Affine' 'Manifold'.
joinAffine
    :: (GVector v x, Map f)
    => Point v d (Codomain f) x
    -> Point v (Function c d) f x
    -> Point v (Function c d) (Affine f) x
{-# INLINE joinAffine #-}
joinAffine (Point cps) (Point cpqs) = Point $ cps G.++ cpqs


--- Instances ---

-- Tensors --

instance (Manifold m, Manifold n) => Manifold (Product m n) where
    type Dimension (Product m n) = Dimension m * Dimension n

instance (Manifold m, Manifold n) => Map (Product m n) where
    type Domain (Product m n) = n
    type Codomain (Product m n) = m

instance (Manifold m, Manifold n, Numeric x) => Apply SVector c d (Product m n) x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs =
        S.map Point . S.toColumns . S.matrixMatrixMultiply (toMatrix pq) . S.fromColumns $ S.map coordinates qs

instance (Manifold m, Manifold n, Numeric x) => Bilinear SVector c d (Product m n) x where
    {-# INLINE (<.<) #-}
    (<.<) q pq = Point $ S.matrixVectorMultiply (toMatrix $ transpose pq) $ coordinates q
    {-# INLINE (<$<) #-}
    (<$<) qs pq =
        S.map Point . S.toColumns . S.matrixMatrixMultiply (toMatrix $ transpose pq) . S.fromColumns $ S.map coordinates qs

-- Affine Maps --

instance Map f => Manifold (Affine f) where
    type Dimension (Affine f) = Dimension (Codomain f) + Dimension f

instance Map f => Map (Affine f) where
    type Domain (Affine f) = Domain f
    type Codomain (Affine f) = Codomain f

instance (Bilinear v c d f x, Num x) => Apply v c d (Affine f) x where
    {-# INLINE (>.>) #-}
    (>.>) ppq q =
        let (p,pq) = splitAffine ppq
         in p <+> pq >.> q
    {-# INLINE (>$>) #-}
    (>$>) ppq qs =
        let (p,pq) = splitAffine ppq
         in G.map (p <+>) (pq >$> qs)

instance (Manifold m, Manifold n, Dimension m ~ Dimension n, KnownNat (Dimension n), Field x)
  => Square SVector m n x where
    {-# INLINE inverse #-}
    inverse p = fromMatrix . S.inverse $ toMatrix p
    {-# INLINE determinant #-}
    determinant = S.determinant . toMatrix




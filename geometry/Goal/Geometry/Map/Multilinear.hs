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

class (Dimension m ~ Dimension n) => Square m n v x where
-- | The inverse of a tensor.
    inverse
        :: Point (Function c d) (Product m n) v x
        -> Point (Function d c) (Product n m) v x
    determinant :: Point (Function c d) (Product m n) v x -> x


-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class Apply c d f v x => Bilinear c d f v x where
    (<.<)
        :: Point (Dual d) (Codomain f) v x
        -> Point (Function c d) f v x
        -> Point (Dual c) (Domain f) v x
    (<$<) :: KnownNat k
        => Vector v k (Point (Dual d) (Codomain f) v x)
        -> Point (Function c d) f v x
        -> Vector v k (Point (Dual c) (Domain f) v x)

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

-- | Tensor v x Tensor multiplication.
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
    => Point (Function c d) (Affine f) v x
    -> (Point d (Codomain f) v x, Point (Function c d) f v x)
{-# INLINE splitAffine #-}
splitAffine (Point cppqs) =
    let (cps,cpqs) = G.splitAt cppqs
     in (Point cps, Point cpqs)

-- | Join a translation and a tensor into a 'Point' on an 'Affine' 'Manifold'.
joinAffine
    :: (GVector v x, Map f)
    => Point d (Codomain f) v x
    -> Point (Function c d) f v x
    -> Point (Function c d) (Affine f) v x
{-# INLINE joinAffine #-}
joinAffine (Point cps) (Point cpqs) = Point $ cps G.++ cpqs


--- Instances ---

-- Tensors --

instance (Manifold m, Manifold n) => Manifold (Product m n) where
    type Dimension (Product m n) = Dimension m * Dimension n

instance (Manifold m, Manifold n) => Map (Product m n) where
    type Domain (Product m n) = n
    type Codomain (Product m n) = m

instance (Manifold m, Manifold n, Numeric x) => Apply c d (Product m n) SVector x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs =
        S.map Point . S.toColumns . S.matrixMatrixMultiply (toMatrix pq) . S.fromColumns $ S.map coordinates qs

instance (Manifold m, Manifold n, Numeric x) => Bilinear c d (Product m n) SVector  x where
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

instance (Bilinear c d f v x, Num x) => Apply c d (Affine f) v x where
    {-# INLINE (>.>) #-}
    (>.>) ppq q =
        let (p,pq) = splitAffine ppq
         in p <+> pq >.> q
    {-# INLINE (>$>) #-}
    (>$>) ppq qs =
        let (p,pq) = splitAffine ppq
         in G.map (p <+>) (pq >$> qs)

instance (Manifold m, Manifold n, Dimension m ~ Dimension n, KnownNat (Dimension n), Field x)
  => Square m n SVector x where
    {-# INLINE inverse #-}
    inverse p = fromMatrix . S.inverse $ toMatrix p
    {-# INLINE determinant #-}
    determinant = S.determinant . toMatrix




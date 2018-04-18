{-# LANGUAGE UndecidableInstances #-}
-- | This module provides tools for working with tensors, affine transformations, and general
-- multilinear objects.

module Goal.Geometry.Map.Multilinear
    ( -- * Bilinear Forms
    Bilinear ((<.<),(<$<))
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


-- Bilinear Forms --

-- | The inverse of a tensor.
inverse
    :: (Manifold m, Manifold n, Dimension m ~ Dimension n)
    => Point (Function c d) (Product m n)
    -> Point (Function d c) (Product n m)
{-# INLINE inverse #-}
inverse p = fromMatrix . S.inverse $ toMatrix p

determinant
    :: (Manifold m, Manifold n, Dimension m ~ Dimension n)
    => Point (Function c d) (Product m n)
    -> Double
{-# INLINE determinant #-}
determinant = S.determinant . toMatrix

-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class Apply c d f => Bilinear c d f where
    (<.<) :: Point (Dual d) (Codomain f)
          -> Point (Function c d) f
          -> Point (Dual c) (Domain f)
    (<$<) :: KnownNat k
          => Point (Dual d) (Replicated k (Codomain f))
          -> Point (Function c d) f
          -> Point (Dual c) (Replicated k (Domain f))
    (<$<) xs f = mapReplicatedPoint (<.< f) xs

-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Product m n

-- | Infix synonym for a 'Tensor'.
--type (m * n) = Product m n
--infixr 6 *

-- | Converts a point on a 'Product' manifold into a 'Matrix'.
toMatrix :: (Manifold m, Manifold n) => Point c (Product m n) -> S.Matrix (Dimension m) (Dimension n) Double
{-# INLINE toMatrix #-}
toMatrix (Point xs) = G.Matrix xs

-- | Converts a point on a 'Product' manifold into a 'Matrix'.
replicatedToMatrix :: (Manifold m, KnownNat k)
                   => Point c (Replicated k m)
                   -> S.Matrix k (Dimension m) Double
{-# INLINE replicatedToMatrix #-}
replicatedToMatrix (Point xs) = G.Matrix xs

-- | Converts a 'Matrix' into a 'Point' on a 'Product' 'Manifold'.
fromMatrix :: S.Matrix (Dimension m) (Dimension n) Double -> Point c (Product m n)
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs

-- | Converts a point on a 'Product' manifold into a 'Matrix'.
replicatedFromMatrix :: (Manifold m, KnownNat k)
                     => S.Matrix k (Dimension m) Double
                     -> Point c (Replicated k m)
{-# INLINE replicatedFromMatrix #-}
replicatedFromMatrix (G.Matrix xs) = Point xs

-- | The transpose of a tensor.
transpose
    :: (Manifold m, Manifold n)
    => Point (Function c d) (Product m n)
    -> Point (Function (Dual d) (Dual c)) (Product n m)
{-# INLINE transpose #-}
transpose (Point xs) = fromMatrix . S.transpose $ G.Matrix xs

-- | Tensor Tensor multiplication.
(<#>) :: (Manifold m, Manifold n, Manifold o)
      => Point (Function d e) (Product m n)
      -> Point (Function c d) (Product n o)
      -> Point (Function c e) (Product m o)
{-# INLINE (<#>) #-}
(<#>) m1 m2 =
    fromMatrix $ S.matrixMatrixMultiply (toMatrix m1) (toMatrix m2)

-- | '>.<' denotes the outer product between two points. It provides a way of
-- constructing matrices of the 'Tensor' product space.
(>.<) :: (Manifold m, Manifold n)
      => Point d m
      -> Point c n
      -> Point (Function (Dual c) d) (Product m n)
{-# INLINE (>.<) #-}
(>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs


--- Affine Functions ---

-- | An 'Affine' 'Manifold' represents linear transformations followed by a translation.
data Affine f

-- | Infix synonym for 'Affine'.
type (m <* n) = Affine (Product m n)
infixr 6 <*

-- | Split a 'Point' on an 'Affine' 'Manifold' into a 'Point' which represents the translation, and a tensor.
splitAffine
    :: Map f
    => Point (Function c d) (Affine f)
    -> (Point d (Codomain f), Point (Function c d) f)
{-# INLINE splitAffine #-}
splitAffine (Point cppqs) =
    let (cps,cpqs) = S.splitAt cppqs
     in (Point cps, Point cpqs)

-- | Join a translation and a tensor into a 'Point' on an 'Affine' 'Manifold'.
joinAffine
    :: Map f
    => Point d (Codomain f)
    -> Point (Function c d) f
    -> Point (Function c d) (Affine f)
{-# INLINE joinAffine #-}
joinAffine (Point cps) (Point cpqs) = Point $ cps S.++ cpqs


--- Instances ---

-- Tensors --

instance (Manifold m, Manifold n) => Manifold (Product m n) where
    type Dimension (Product m n) = Dimension m * Dimension n

instance (Manifold m, Manifold n) => Map (Product m n) where
    type Domain (Product m n) = n
    type Codomain (Product m n) = m

instance (Manifold m, Manifold n) => Apply c d (Product m n) where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs =
        replicatedFromMatrix . S.transpose . S.matrixMatrixMultiply (toMatrix pq) . S.transpose $ replicatedToMatrix qs

instance (Manifold m, Manifold n) => Bilinear c d (Product m n) where
    {-# INLINE (<.<) #-}
    (<.<) q pq = Point $ S.matrixVectorMultiply (toMatrix $ transpose pq) $ coordinates q
    {-# INLINE (<$<) #-}
    (<$<) qs pq =
        replicatedFromMatrix . S.matrixMatrixMultiply (replicatedToMatrix qs) $ toMatrix pq


-- Affine Maps --


instance Map f => Manifold (Affine f) where
    type Dimension (Affine f) = Dimension (Codomain f) + Dimension f

instance Map f => Map (Affine f) where
    type Domain (Affine f) = Domain f
    type Codomain (Affine f) = Codomain f

instance Bilinear c d f => Apply c d (Affine f) where
    {-# INLINE (>.>) #-}
    (>.>) ppq q =
        let (p,pq) = splitAffine ppq
         in p <+> pq >.> q
    {-# INLINE (>$>) #-}
    (>$>) ppq qs =
        let (p,pq) = splitAffine ppq
         in mapReplicatedPoint (p <+>) (pq >$> qs)

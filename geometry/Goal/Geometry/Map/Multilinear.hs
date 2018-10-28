{-# LANGUAGE UndecidableInstances #-}
-- | This module provides tools for working with tensors, affine transformations, and general
-- multilinear objects.

module Goal.Geometry.Map.Multilinear
    ( -- * Bilinear Forms
    Bilinear ((<.<),(<$<),(>.<),transpose)
    -- * Tensors
    , Tensor
    -- ** Matrix Construction
    , toMatrix
    , fromMatrix
    , toRows
    , fromRows
    -- ** Computation
    , (<#>)
    , inverse
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


-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class (Manifold m, Manifold n, Manifold (f m n)) => Bilinear f m n where
    -- | Transposed application.
    (<.<) :: (Manifold m, Manifold n)
          => Dual d # m
          -> Function c d # f m n
          -> Dual c # n
    -- | Mapped transposed application.
    (<$<) :: (Manifold m, Manifold n)
          => [Dual d # m]
          -> Function c d # f m n
          -> [Dual c # n]
    --(<$<) xs f = mapReplicatedPoint (<.< f) xs
    -- | Tensor outer product.
    (>.<) :: (Manifold m, Manifold n)
          => d # m
          -> c # n
          -> Function (Dual c) d # f m n
    -- | Tensor transpose.
    transpose :: (Manifold m, Manifold n)
              => c #> d # f m n
              -> Dual d #> Dual c # f n m


-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Tensor m n

-- | The inverse of a tensor.
inverse
    :: (Manifold m, Manifold n, Dimension m ~ Dimension n)
    => Point (Function c d) (Tensor m n)
    -> Point (Function d c) (Tensor n m)
{-# INLINE inverse #-}
inverse p = fromMatrix . S.inverse $ toMatrix p

-- | The determinant of a tensor.
determinant
    :: (Manifold m, Manifold n, Dimension m ~ Dimension n)
    => Point (Function c d) (Tensor m n)
    -> Double
{-# INLINE determinant #-}
determinant = S.determinant . toMatrix

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
toMatrix :: (Manifold m, Manifold n) => Point c (Tensor m n) -> S.Matrix (Dimension m) (Dimension n) Double
{-# INLINE toMatrix #-}
toMatrix (Point xs) = G.Matrix xs

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
toRows :: (Manifold m, Manifold n) => Dual c #> c # Tensor m n -> S.Vector (Dimension m) (c # n)
{-# INLINE toRows #-}
toRows tns = S.map Point . S.toRows $ toMatrix tns

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
fromRows :: (Manifold m, Manifold n) => S.Vector (Dimension m) (c # n) -> Dual c #> c # Tensor m n
{-# INLINE fromRows #-}
fromRows rws = fromMatrix . S.fromRows $ S.map coordinates rws

-- | Converts a 'Matrix' into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension m) (Dimension n) Double -> Point c (Tensor m n)
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
--replicatedToMatrix :: (Manifold m, KnownNat k)
--                   => Point c (Replicated k m)
--                   -> S.Matrix k (Dimension m) Double
--{-# INLINE replicatedToMatrix #-}
--replicatedToMatrix (Point xs) = G.Matrix xs
--
---- | Converts a point on a 'Tensor manifold into a 'Matrix'.
--replicatedFromMatrix :: (Manifold m, KnownNat k)
--                     => S.Matrix k (Dimension m) Double
--                     -> Point c (Replicated k m)
--{-# INLINE replicatedFromMatrix #-}
--replicatedFromMatrix (G.Matrix xs) = Point xs

-- | Tensor Tensor multiplication.
(<#>) :: (Manifold m, Manifold n, Manifold o)
      => Point (Function d e) (Tensor m n)
      -> Point (Function c d) (Tensor n o)
      -> Point (Function c e) (Tensor m o)
{-# INLINE (<#>) #-}
(<#>) m1 m2 =
    fromMatrix $ S.matrixMatrixMultiply (toMatrix m1) (toMatrix m2)


--- Affine Functions ---


-- | An 'Affine' 'Manifold' represents linear transformations followed by a translation.
data Affine (f :: * -> * -> *) m n

-- | Infix synonym for 'Affine'.
type (m <* n) = Affine Tensor m n
infixr 6 <*

-- | Split a 'Point' on an 'Affine' 'Manifold' into a 'Point' which represents the translation, and a tensor.
splitAffine
    :: (Manifold m, Manifold n, Manifold (f m n))
    => Function c d # Affine f m n
    -> (d # m, Function c d # f m n)
{-# INLINE splitAffine #-}
splitAffine (Point cppqs) =
    let (cps,cpqs) = S.splitAt cppqs
     in (Point cps, Point cpqs)

-- | Join a translation and a tensor into a 'Point' on an 'Affine' 'Manifold'.
joinAffine
    :: (Manifold m, Manifold n)
    => d # m
    -> Function c d # f m n
    -> Function c d # Affine f m n
{-# INLINE joinAffine #-}
joinAffine (Point cps) (Point cpqs) = Point $ cps S.++ cpqs


--- Instances ---

-- Tensors --

instance (Manifold m, Manifold n) => Manifold (Tensor m n) where
    type Dimension (Tensor m n) = Dimension m * Dimension n

instance (Manifold m, Manifold n) => Map c d Tensor m n where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = Point <$> S.matrixMap (toMatrix pq) (coordinates <$> qs)
    --{-# INLINE (>#>) #-}
    --(>#>) pq qs =
    --    replicatedFromMatrix . S.transpose . S.matrixMatrixMultiply (toMatrix pq) . S.transpose $ replicatedToMatrix qs

instance (Manifold m, Manifold n) => Bilinear Tensor m n where
    {-# INLINE (<.<) #-}
    (<.<) q pq = Point $ S.matrixVectorMultiply (toMatrix $ transpose pq) $ coordinates q
    {-# INLINE (<$<) #-}
    (<$<) qs pq = map Point . S.matrixMap (toMatrix $ transpose pq) $ coordinates <$> qs
    --{-# INLINE (<#<) #-}
    --(<#<) qs pq =
    --    replicatedFromMatrix . S.matrixMatrixMultiply (replicatedToMatrix qs) $ toMatrix pq
    {-# INLINE (>.<) #-}
    (>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs
    {-# INLINE transpose #-}
    transpose (Point xs) = fromMatrix . S.transpose $ G.Matrix xs


-- Affine Maps --


instance (Manifold m, Manifold n, Manifold (f m n)) => Manifold (Affine f m n) where
    type Dimension (Affine f m n) = Dimension m + Dimension (f m n)

instance Map c d f m n => Map c d (Affine f) m n where
    {-# INLINE (>.>) #-}
    (>.>) ppq q =
        let (p,pq) = splitAffine ppq
         in p <+> pq >.> q
    {-# INLINE (>$>) #-}
    (>$>) ppq qs =
        let (p,pq) = splitAffine ppq
         in (p <+>) <$> (pq >$> qs)
    --{-# INLINE (>#>) #-}
    --(>#>) ppq qs =
    --    let (p,pq) = splitAffine ppq
    --     in mapReplicatedPoint (p <+>) (pq >#> qs)

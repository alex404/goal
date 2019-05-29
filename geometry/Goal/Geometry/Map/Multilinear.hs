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
class (Manifold x, Manifold y, Manifold (f y x)) => Bilinear f y x where
    -- | Transposed application.
    (<.<) :: (Manifold x, Manifold y) => d #* y -> c #> d # f y x -> c #* x
    -- | Mapped transposed application.
    (<$<) :: (Manifold x, Manifold y) => [d #* y] -> c #> d # f y x -> [c #* x]
    -- | Tensor outer product.
    (>.<) :: (Manifold x, Manifold y) => d # y -> c # x -> Dual c #> d # f y x
    -- | Tensor transpose.
    transpose :: (Manifold x, Manifold y) => c #> d # f y x -> d #> c #* f x y


-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Tensor y x

-- | The inverse of a tensor.
inverse
    :: (Manifold x, Manifold y, Dimension x ~ Dimension y)
    => c #> d # Tensor y x -> d #> c # Tensor x y
{-# INLINE inverse #-}
inverse p = fromMatrix . S.inverse $ toMatrix p

-- | The determinant of a tensor.
determinant
    :: (Manifold x, Manifold y, Dimension x ~ Dimension y)
    => c #> d # Tensor y x
    -> Double
{-# INLINE determinant #-}
determinant = S.determinant . toMatrix

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
toMatrix :: (Manifold x, Manifold y) => c # Tensor y x -> S.Matrix (Dimension y) (Dimension x) Double
{-# INLINE toMatrix #-}
toMatrix (Point xs) = G.Matrix xs

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
toRows :: (Manifold x, Manifold y) => Dual c #> c # Tensor y x -> S.Vector (Dimension y) (c # x)
{-# INLINE toRows #-}
toRows tns = S.map Point . S.toRows $ toMatrix tns

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
fromRows :: (Manifold x, Manifold y) => S.Vector (Dimension y) (c # x) -> Dual c #> c # Tensor y x
{-# INLINE fromRows #-}
fromRows rws = fromMatrix . S.fromRows $ S.map coordinates rws

-- | Converts a 'Matrix' into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension y) (Dimension x) Double -> c # Tensor y x
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
--replicatedToMatrix :: (Manifold x, KnownNat k)
--                   => Point c (Replicated k m)
--                   -> S.Matrix k (Dimension x) Double
--{-# INLINE replicatedToMatrix #-}
--replicatedToMatrix (Point xs) = G.Matrix xs
--
---- | Converts a point on a 'Tensor manifold into a 'Matrix'.
--replicatedFromMatrix :: (Manifold x, KnownNat k)
--                     => S.Matrix k (Dimension x) Double
--                     -> Point c (Replicated k m)
--{-# INLINE replicatedFromMatrix #-}
--replicatedFromMatrix (G.Matrix xs) = Point xs

-- | Tensor Tensor multiplication.
(<#>) :: (Manifold x, Manifold y, Manifold z)
      => d #> e # Tensor z y
      -> c #> d # Tensor y x
      -> c #> e # Tensor z x
{-# INLINE (<#>) #-}
(<#>) m1 m2 =
    fromMatrix $ S.matrixMatrixMultiply (toMatrix m1) (toMatrix m2)


--- Affine Functions ---


-- | An 'Affine' 'Manifold' represents linear transformations followed by a translation.
data Affine (f :: Type -> Type -> Type) y x

-- | Infix synonym for 'Affine'.
type (y <* x) = Affine Tensor y x
infixr 6 <*

-- | Split a 'Point' on an 'Affine' 'Manifold' into a 'Point' which represents the translation, and a tensor.
splitAffine
    :: (Manifold x, Manifold y, Manifold (f y x))
    => c #> d # Affine f y x
    -> (d # y, c #> d # f y x)
{-# INLINE splitAffine #-}
splitAffine (Point cppqs) =
    let (cps,cpqs) = S.splitAt cppqs
     in (Point cps, Point cpqs)

-- | Join a translation and a tensor into a 'Point' on an 'Affine' 'Manifold'.
joinAffine
    :: (Manifold x, Manifold y)
    => d # y
    -> c #> d # f y x
    -> c #> d # Affine f y x
{-# INLINE joinAffine #-}
joinAffine (Point cps) (Point cpqs) = Point $ cps S.++ cpqs


--- Instances ---

-- Tensors --

instance (Manifold x, Manifold y) => Manifold (Tensor y x) where
    type Dimension (Tensor y x) = Dimension x * Dimension y

instance (Manifold x, Manifold y) => Map c d Tensor y x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = Point <$> S.matrixMap (toMatrix pq) (coordinates <$> qs)
    --{-# INLINE (>#>) #-}
    --(>#>) pq qs =
    --    replicatedFromMatrix . S.transpose . S.matrixMatrixMultiply (toMatrix pq) . S.transpose $ replicatedToMatrix qs

instance (Manifold x, Manifold y) => Bilinear Tensor y x where
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


instance (Manifold x, Manifold y, Manifold (f y x)) => Manifold (Affine f y x) where
    type Dimension (Affine f y x) = Dimension y + Dimension (f y x)

instance Map c d f y x => Map c d (Affine f) y x where
    {-# INLINE (>.>) #-}
    (>.>) ppq q =
        let (p,pq) = splitAffine ppq
         in p <+> pq >.> q
    {-# INLINE (>$>) #-}
    (>$>) ppq qs =
        let (p,pq) = splitAffine ppq
         in (p <+>) <$> (pq >$> qs)

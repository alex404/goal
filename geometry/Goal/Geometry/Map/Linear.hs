{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances,UndecidableSuperClasses #-}
-- | This module provides tools for working with tensors, affine transformations, and general
-- multilinear objects.

module Goal.Geometry.Map.Linear
    ( -- * Bilinear Forms
    Bilinear ((>$<),(>.<),transpose)
    , (<.<)
    , (<$<)
    -- * Tensors
    , Tensor
    -- ** Matrix Construction
    , toMatrix
    , fromMatrix
    , toRows
    , toColumns
    , fromRows
    , fromColumns
    -- ** Computation
    --, (<#>)
    , inverse
    , determinant
    -- * Affine Functions
    , Affine
    , Translation ((>+>),anchor)
    , (>.+>)
    , (>$+>)
    , type (<*)
    ) where

--- Imports ---

-- Package --

import Goal.Core

import Goal.Geometry.Manifold
import Goal.Geometry.Vector
import Goal.Geometry.Map

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G


-- Bilinear Forms --


-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class (Bilinear f x y, Manifold x, Manifold y, Manifold (f x y)) => Bilinear f y x where
    -- | Tensor outer product.
    (>.<) :: c # y -> c # x -> c # f y x
    -- | Average of tensor outer products.
    (>$<) :: [c # y] -> [c # x] -> c # f y x
    -- | Tensor transpose.
    transpose :: c # f y x -> c # f x y

-- | Transposed application.
(<.<) :: (Map c f x y, Bilinear f y x) => c #* y -> c # f y x -> c # x
{-# INLINE (<.<) #-}
(<.<) dy f = transpose f >.> dy

-- | Mapped transposed application.
(<$<) :: (Map c f x y, Bilinear f y x) => [c #* y] -> c # f y x -> [c # x]
{-# INLINE (<$<) #-}
(<$<) dy f = transpose f >$> dy


-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Tensor y x

-- | The inverse of a tensor.
inverse
    :: (Manifold x, Manifold y, Dimension x ~ Dimension y)
    => c # Tensor y x -> c #* Tensor x y
{-# INLINE inverse #-}
inverse p = fromMatrix . S.pseudoInverse $ toMatrix p

-- | The determinant of a tensor.
determinant
    :: (Manifold x, Manifold y, Dimension x ~ Dimension y)
    => c # Tensor y x
    -> Double
{-# INLINE determinant #-}
determinant = S.determinant . toMatrix

-- | Converts a point on a 'Tensor manifold into a 'Matrix'.
toMatrix :: (Manifold x, Manifold y) => c # Tensor y x -> S.Matrix (Dimension y) (Dimension x) Double
{-# INLINE toMatrix #-}
toMatrix (Point xs) = G.Matrix xs

-- | Converts a point on a 'Tensor' manifold into a a vector of rows.
toRows :: (Manifold x, Manifold y) => c # Tensor y x -> S.Vector (Dimension y) (c # x)
{-# INLINE toRows #-}
toRows tns = S.map Point . S.toRows $ toMatrix tns

-- | Converts a point on a 'Tensor' manifold into a a vector of rows.
toColumns :: (Manifold x, Manifold y) => c # Tensor y x -> S.Vector (Dimension x) (c # y)
{-# INLINE toColumns #-}
toColumns tns = S.map Point . S.toColumns $ toMatrix tns

-- | Converts a vector of rows into a 'Tensor'.
fromRows :: (Manifold x, Manifold y) => S.Vector (Dimension y) (c # x) -> c # Tensor y x
{-# INLINE fromRows #-}
fromRows rws = fromMatrix . S.fromRows $ S.map coordinates rws

-- | Converts a vector of rows into a 'Tensor'.
fromColumns :: (Manifold x, Manifold y) => S.Vector (Dimension x) (c # y) -> c # Tensor y x
{-# INLINE fromColumns #-}
fromColumns rws = fromMatrix . S.fromColumns $ S.map coordinates rws

-- | Converts a 'Matrix' into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension y) (Dimension x) Double -> c # Tensor y x
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs

---- | Tensor x Tensor multiplication.
--(<#>) :: (Manifold x, Manifold y, Manifold z)
--      => Function d e # Tensor z y
--      -> Function c d # Tensor y x
--      -> Function c e # Tensor z x
--{-# INLINE (<#>) #-}
--(<#>) m1 m2 =
--    fromMatrix $ S.matrixMatrixMultiply (toMatrix m1) (toMatrix m2)


--- Affine Functions ---


-- | An 'Affine' 'Manifold' represents linear transformations followed by a translation.
newtype Affine f y z x = Affine (z,f y x)

deriving instance (Manifold z, Manifold (f y x)) => Manifold (Affine f y z x)
deriving instance (Manifold z, Manifold (f y x)) => Product (Affine f y z x)

-- | Infix synonym for 'Affine'.
type (y <* x) = Affine Tensor y y x
infixr 6 <*

class (Manifold y, Manifold z) => Translation z y where
    (>+>) :: c # z -> c # y -> c # z
    anchor :: c # z -> c # y

-- | Transposed application.

(>.+>) :: (Map c f y x, Translation z x) => c # f y x -> c #* z -> c # y
(>.+>) f w = f >.> anchor w

(>$+>) :: (Map c f y x, Translation z x) => c # f y x -> [c #* z] -> [c # y]
(>$+>) f w = f >$> (anchor <$> w)


--- Instances ---

-- Tensors --

instance (Manifold x, Manifold y) => Manifold (Tensor y x) where
    type Dimension (Tensor y x) = Dimension x * Dimension y

instance (Manifold x, Manifold y) => Map c Tensor y x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = Point <$> S.matrixMap (toMatrix pq) (coordinates <$> qs)

instance (Manifold x, Manifold y) => Bilinear Tensor y x where
    {-# INLINE (>.<) #-}
    (>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs
    {-# INLINE (>$<) #-}
    (>$<) ps qs = fromMatrix . S.averageOuterProduct $ zip (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose (Point xs) = fromMatrix . S.transpose $ G.Matrix xs


-- Affine Maps --

instance Manifold z => Translation z z where
    (>+>) z1 z2 = z1 + z2
    anchor = id

instance (Manifold z, Manifold y) => Translation (y,z) y where
    (>+>) yz y' =
        let (y,z) = split yz
         in join (y + y') z
    anchor = fst . split

instance (Translation z y, Map c f y x) => Map c (Affine f y) z x where
    {-# INLINE (>.>) #-}
    (>.>) fyzx x =
        let (yz,yx) = split fyzx
         in   yz >+> (yx >.> x)
    (>$>) fyzx xs =
        let (yz,yx) = split fyzx
         in (yz >+>) <$> yx >$> xs

instance (KnownNat n, Translation w z)
  => Translation (Replicated n w) (Replicated n z) where
      {-# INLINE (>+>) #-}
      (>+>) w z =
          let ws = splitReplicated w
              zs = splitReplicated z
           in joinReplicated $ S.zipWith (>+>) ws zs
      {-# INLINE anchor #-}
      anchor = mapReplicatedPoint anchor


--instance (Map c f z x) => Map c (Affine f z) z x where
--    {-# INLINE (>.>) #-}
--    (>.>) ppq q =
--        let (p,pq) = split ppq
--         in p + pq >.> q
--    {-# INLINE (>$>) #-}
--    (>$>) ppq qs =
--        let (p,pq) = split ppq
--         in (p +) <$> (pq >$> qs)

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances,UndecidableSuperClasses,TypeApplications #-}
-- | This module provides tools for working with linear and affine
-- transformations.

module Goal.Geometry.Map.Linear
    ( -- * Bilinear Forms
    Bilinear ((>$<),(>.<),transpose, toTensor)
    , (<.<)
    , (<$<)
    , Square (inverse, determinant, inverseLogDeterminant)
    -- * Tensors
    , Tensor
    , Symmetric
    , Diagonal
    , Scale
    -- ** Matrix Construction
    --, toMatrix
    --, fromMatrix
    , toRows
    , toColumns
    , fromRows
    , fromColumns
    -- ** Computation
    --, (<#>)
    -- * Affine Functions
    , Affine (Affine)
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
class (Bilinear f y x, Manifold x, Manifold y, Manifold (f y x)) => Bilinear f x y where
    -- | Tensor outer product.
    (>.<) :: c # x -> c # y -> c # f x y
    -- | Average of tensor outer products.
    (>$<) :: [c # x] -> [c # y] -> c # f x y
    -- | Tensor transpose.
    transpose :: c # f x y -> c # f y x
    -- | Convert to a full tensor
    toTensor :: c # f x y -> c # Tensor x y

-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class (Dimension x ~ Dimension y, Manifold x, Manifold y, Manifold (f x y)) => Square f x y where
    -- | The determinant of a tensor.
    inverse :: c # f x y -> c #* f y x
    -- | The inverse of a tensor.
    determinant :: c # f x y -> Double
    inverseLogDeterminant :: c # f x y -> (c #* f y x, Double, Double)

-- | Transposed application.
(<.<) :: (Map c f y x, Bilinear f x y) => c #* x -> c # f x y -> c # y
{-# INLINE (<.<) #-}
(<.<) dx f = transpose f >.> dx

-- | Mapped transposed application.
(<$<) :: (Map c f y x, Bilinear f x y) => [c #* x] -> c # f x y -> [c # y]
{-# INLINE (<$<) #-}
(<$<) dx f = transpose f >$> dx


-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Tensor y x

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Symmetric x y

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Diagonal x y

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Scale x y


-- | Converts a point on a 'Tensor manifold into a Matrix.
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

-- | Converts a Matrix into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension y) (Dimension x) Double -> c # Tensor y x
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs


--- Affine Functions ---


-- | An 'Affine' 'Manifold' represents linear transformations followed by a
-- translation. The 'First' component is the translation, and the 'Second'
-- component is the linear transformation.
newtype Affine f y z x = Affine (z,f y x)

deriving instance (Manifold z, Manifold (f y x)) => Manifold (Affine f y z x)
deriving instance (Manifold z, Manifold (f y x)) => Product (Affine f y z x)

-- | Infix synonym for simple 'Affine' transformations.
type (y <* x) = Affine Tensor y y x
infixr 6 <*

-- | The 'Translation' class is used to define translations where we only want
-- to translate a subset of the parameters of the given object.
class (Manifold y, Manifold z) => Translation z y where
    -- | Translates the the first argument by the second argument.
    (>+>) :: c # z -> c # y -> c # z
    -- | Returns the subset of the parameters of the given 'Point' that are
    -- translated in this instance.
    anchor :: c # z -> c # y

-- | Operator that applies a 'Map' to a subset of an input's parameters.
(>.+>) :: (Map c f y x, Translation z x) => c # f y x -> c #* z -> c # y
(>.+>) f w = f >.> anchor w

-- | Operator that maps a 'Map' over a subset of the parameters of a list of inputs.
(>$+>) :: (Map c f y x, Translation z x) => c # f y x -> [c #* z] -> [c # y]
(>$+>) f w = f >$> (anchor <$> w)


---- Instances ----


--- Tensors ---


--- Manifold

instance (Manifold x, Manifold y) => Manifold (Tensor y x) where
    type Dimension (Tensor y x) = Dimension x * Dimension y

instance (Manifold x, KnownNat (Triangular (Dimension x))) => Manifold (Symmetric x x) where
    type Dimension (Symmetric x x) = Triangular (Dimension x)

instance Manifold x => Manifold (Diagonal x x) where
    type Dimension (Diagonal x x) = Dimension x

instance Manifold x => Manifold (Scale x x) where
    type Dimension (Scale x x) = 1


--- Map

instance (Manifold x, Manifold y) => Map c Tensor y x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (toMatrix pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = Point <$> S.matrixMap (toMatrix pq) (coordinates <$> qs)

instance (Manifold x, Manifold (Symmetric x x)) => Map c Symmetric x x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ S.matrixVectorMultiply (S.fromLowerTriangular $ coordinates pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = Point <$> S.matrixMap (S.fromLowerTriangular $ coordinates pq) (coordinates <$> qs)

instance (Manifold x, Manifold (Diagonal x x)) => Map c Diagonal x x where
    {-# INLINE (>.>) #-}
    (>.>) (Point xs) (Point ys) = Point $ xs * ys
    {-# INLINE (>$>) #-}
    (>$>) (Point xs) yss = Point <$> S.diagonalMatrixMap xs (coordinates <$> yss)

instance (Manifold x, Manifold (Scale x x)) => Map c Scale x x where
    {-# INLINE (>.>) #-}
    (>.>) (Point x) (Point y) = Point $ S.scale (S.head x) y
    {-# INLINE (>$>) #-}
    (>$>) (Point x) yss = Point . S.scale (S.head x) . coordinates <$> yss


--- Bilinear

instance (Manifold x, Manifold y) => Bilinear Tensor x y where
    {-# INLINE (>.<) #-}
    (>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs
    {-# INLINE (>$<) #-}
    (>$<) ps qs = fromMatrix . S.averageOuterProduct $ zip (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose (Point xs) = fromMatrix . S.transpose $ G.Matrix xs
    toTensor = id

instance Manifold x => Bilinear Symmetric x x where
    {-# INLINE (>.<) #-}
    (>.<) (Point xs) (Point ys) = Point . S.lowerTriangular $ xs `S.outerProduct` ys
    {-# INLINE (>$<) #-}
    (>$<) ps qs = Point . S.lowerTriangular . S.averageOuterProduct $ zip (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose = id
    toTensor (Point mtx) = fromMatrix $ S.fromLowerTriangular mtx

instance Manifold x => Bilinear Diagonal x x where
    {-# INLINE (>.<) #-}
    (>.<) (Point xs) (Point ys) = Point $ xs * ys
    {-# INLINE (>$<) #-}
    (>$<) ps qs = Point . average $ zipWith (*) (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose = id
    toTensor (Point diag) = fromMatrix $ S.diagonalMatrix diag

instance Manifold x => Bilinear Scale x x where
    {-# INLINE (>.<) #-}
    (>.<) (Point x) (Point y) = singleton $ S.dotProduct x y
    {-# INLINE (>$<) #-}
    (>$<) ps qs = singleton . average $ zipWith S.dotProduct (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose = id
    toTensor (Point scl) = fromMatrix . S.diagonalMatrix $ S.scale (S.head scl) 1


-- Square --


instance (Manifold x, Manifold y, Dimension x ~ Dimension y) => Square Tensor y x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse p = fromMatrix . S.pseudoInverse $ toMatrix p
    {-# INLINE determinant #-}
    determinant = S.determinant . toMatrix
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant tns =
        let (imtx,det,sgn) = S.inverseLogDeterminant $ toMatrix tns
         in (fromMatrix imtx, det, sgn)

instance Manifold x => Square Symmetric x x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse (Point trng) =
        let mtx :: S.Matrix (Dimension x) (Dimension x) Double
            mtx = S.fromLowerTriangular trng
         in Point . S.lowerTriangular $ S.pseudoInverse mtx
    {-# INLINE determinant #-}
    determinant (Point trng) =
        let mtx :: S.Matrix (Dimension x) (Dimension x) Double
            mtx = S.fromLowerTriangular trng
         in S.determinant mtx
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant (Point trng) =
        let mtx :: S.Matrix (Dimension x) (Dimension x) Double
            mtx = S.fromLowerTriangular trng
            (imtx,det,sgn) = S.inverseLogDeterminant mtx
         in (Point $ S.lowerTriangular imtx, det, sgn)

instance Manifold x => Square Diagonal x x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse (Point diag) = Point $ recip diag
    {-# INLINE determinant #-}
    determinant (Point diag) = S.product diag
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant sqr =
        let diag = coordinates sqr
            lndet = S.sum $ log diag
         in (inverse sqr, abs lndet, signum lndet)

instance Manifold x => Square Scale x x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse (Point scl) = Point $ recip scl
    {-# INLINE determinant #-}
    determinant (Point scl) = S.head scl ^ (natValInt $ Proxy @(Dimension x))
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant sqr =
        let scl = coordinates sqr
            lndet = (fromIntegral . natVal $ Proxy @(Dimension x)) * log (S.head scl)
         in (inverse sqr, abs lndet, signum lndet)


--- Affine Maps ---


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

--instance (KnownNat n, Translation w z)
--  => Translation (Replicated n w) (Replicated n z) where
--      {-# INLINE (>+>) #-}
--      (>+>) w z =
--          let ws = splitReplicated w
--              zs = splitReplicated z
--           in joinReplicated $ S.zipWith (>+>) ws zs
--      {-# INLINE anchor #-}
--      anchor = mapReplicatedPoint anchor


--instance (Map c f z x) => Map c (Affine f z) z x where
--    {-# INLINE (>.>) #-}
--    (>.>) ppq q =
--        let (p,pq) = split ppq
--         in p + pq >.> q
--    {-# INLINE (>$>) #-}
--    (>$>) ppq qs =
--        let (p,pq) = split ppq
--         in (p +) <$> (pq >$> qs)

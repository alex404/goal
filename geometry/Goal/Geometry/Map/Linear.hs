{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances,TypeApplications #-}
-- | This module provides tools for working with linear and affine
-- transformations.

module Goal.Geometry.Map.Linear
    ( -- * Bilinear Maps
    Bilinear ((>$<),(>.<),transpose, toTensor, fromTensor)
    , LinearlyComposable
    , (<.<)
    , (<$<)
    , Square (inverse, matrixRoot, determinant, inverseLogDeterminant)
    , blockSymmetricMatrixInversion
    -- * Tensors
    , Tensor
    , Symmetric
    , Diagonal
    , Scale
    -- ** Matrix Construction
    , toRows
    , toColumns
    , fromRows
    , fromColumns
    -- ** Operations
    , dualComposition
    , changeOfBasis
    , inverseSchurComplement
    , woodburyMatrix
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
class (Map c f x y, Map c f y x) => Bilinear c f x y where
    -- | Tensor outer product.
    (>.<) :: c # x -> c # y -> c # f x y
    {-# INLINE (>.<) #-}
    (>.<) x y = fromTensor $ x >.< y
    -- | Average of tensor outer products.
    (>$<) :: [c # x] -> [c # y] -> c # f x y
    {-# INLINE (>$<) #-}
    (>$<) x y = fromTensor $ x >$< y
    -- | Tensor transpose.
    transpose :: c # f x y -> c # f y x
    -- | Convert to a full Tensor
    toTensor :: c # f x y -> c # Tensor x y
    -- | Convert a full Tensor to an f, probably throwing out a bunch of elements
    fromTensor :: c # Tensor x y -> c # f x y

-- | A 'Manifold' is 'Bilinear' if its elements are bilinear forms.
class (Bilinear (Dual c) f x x, Bilinear c f x x) => Square c f x where
    -- | The determinant of a tensor.
    inverse :: c # f x x -> c #* f x x
    {-# INLINE inverse #-}
    inverse = fromTensor . inverse . toTensor
    -- | The root of a tensor.
    matrixRoot :: c # f x x -> c # f x x
    {-# INLINE matrixRoot #-}
    matrixRoot = fromTensor . matrixRoot . toTensor
    -- | The inverse of a tensor.
    determinant :: c # f x x -> Double
    {-# INLINE determinant #-}
    determinant = determinant . toTensor
    inverseLogDeterminant :: c # f x x -> (c #* f x x, Double, Double)
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant tns =
        let (inv,lndt,sgn) = inverseLogDeterminant $ toTensor tns
         in (fromTensor inv,lndt,sgn)

dualComposition
    :: ( LinearlyComposable c (Dual c) f Tensor w x z
       , LinearlyComposable (Dual c) c g h x y z)
    => c # f w x
    -> c #* g x y
    -> c # h y z
    -> c # Tensor w z
{-# INLINE dualComposition #-}
dualComposition f g h = unsafeMatrixMultiply f $ unsafeMatrixMultiply g h

--(<#>)
--    :: (LinearlyComposable f g x y z, c ~ Dual c)
--    => c # f x y
--    -> c # g y z
--    -> c # Tensor x z
--{-# INLINE (<#>) #-}
--(<#>) f g = unsafeMatrixMultiply f g

-- | Change of basis formula
changeOfBasis
    :: ( LinearlyComposable c (Dual c) f Tensor y x y, LinearlyComposable (Dual c) c g f x x y )
    => c # f x y -> c #* g x x -> c # Tensor y y
{-# INLINE changeOfBasis #-}
changeOfBasis f g = dualComposition (transpose f) g f

inverseSchurComplement
    :: (LinearlyComposable (Dual c) c f Tensor x x y, Square c f x)
    => c # Tensor y y
    -> c # Tensor x y
    -> c # f x x
    -> c #* Tensor y y
{-# INLINE inverseSchurComplement #-}
inverseSchurComplement br tr tl = inverse $ br - changeOfBasis tr (inverse tl)

woodburyMatrix
    :: ( Primal c, Square c f x, Manifold y
       , LinearlyComposable (Dual c) c f Tensor x x x
       , LinearlyComposable c (Dual c) Tensor f x x x )
    => c # f x x
    -> c # Tensor x y
    -> c #* Tensor y y
    -> c #* Tensor x x
{-# INLINE woodburyMatrix #-}
woodburyMatrix tl tr schr =
    let invtl = inverse tl
        crct = changeOfBasis invtl $ changeOfBasis (transpose tr) schr
     in toTensor invtl + crct

blockSymmetricMatrixInversion
    :: ( Primal c, Square c f x
       , LinearlyComposable (Dual c) c f Tensor x x y
       , LinearlyComposable (Dual c) c f Tensor x x x
       , LinearlyComposable c (Dual c) Tensor f x x x )
    => c # f x x
    -> c # Tensor x y
    -> c # Symmetric y y
    -> (c #* Symmetric x x, c #* Tensor x y, c #* Symmetric y y)
{-# INLINE blockSymmetricMatrixInversion #-}
blockSymmetricMatrixInversion tl tr br =
    let tnsy = toTensor br
        shry = inverseSchurComplement tnsy tr tl
        shrx = woodburyMatrix tl tr shry
        tr' = -dualComposition (inverse tl) tr shry
     in (fromTensor shrx, tr', fromTensor shry)

-- | Transposed application.
(<.<) :: Bilinear c f x y => c #* x -> c # f x y -> c # y
{-# INLINE (<.<) #-}
(<.<) dx f = transpose f >.> dx

-- | Mapped transposed application.
(<$<) :: Bilinear c f x y => [c #* x] -> c # f x y -> [c # y]
{-# INLINE (<$<) #-}
(<$<) dx f = transpose f >$> dx

-- Tensor Products --

-- | 'Manifold' of 'Tensor's given by the tensor product of the underlying pair of 'Manifold's.
data Tensor x y

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

-- | Converts a Matrix into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension y) (Dimension x) Double -> c # Tensor y x
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs


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


--- Internal ---


class (Bilinear c f x y, Bilinear d g y z) => LinearlyComposable c d f g x y z where
    unsafeMatrixMultiply :: c # f x y -> d # g y z -> c # Tensor x z
    {-# INLINE unsafeMatrixMultiply #-}
    unsafeMatrixMultiply f g = fromMatrix
        $ S.matrixMatrixMultiply (toMatrix $ toTensor f) (toMatrix $ toTensor g)

instance (Manifold x, Manifold y, Manifold z) => LinearlyComposable c d Tensor Tensor x y z where
instance (Manifold x, Manifold y) => LinearlyComposable c d Tensor Symmetric x y y where
instance (Manifold x, Manifold y) => LinearlyComposable c d Symmetric Tensor x x y where
instance Manifold x => LinearlyComposable c d Symmetric Symmetric x x x where

instance (Manifold x, Manifold y) => LinearlyComposable c d Diagonal Tensor x x y where
    {-# INLINE unsafeMatrixMultiply #-}
    unsafeMatrixMultiply diag tns =
        fromMatrix . S.diagonalMatrixMatrixMultiply (coordinates diag) $ toMatrix tns

instance (Manifold x, Manifold y) => LinearlyComposable c d Tensor Diagonal x y y where
    {-# INLINE unsafeMatrixMultiply #-}
    unsafeMatrixMultiply tns diag =
        fromMatrix . S.transpose . S.diagonalMatrixMatrixMultiply (coordinates diag)
        . toMatrix $ transpose tns

instance (Manifold x, Manifold y) => LinearlyComposable c d Scale Tensor x x y where
    {-# INLINE unsafeMatrixMultiply #-}
    unsafeMatrixMultiply f g = breakPoint $ S.head (coordinates f) .> g

instance (Manifold x, Manifold y) => LinearlyComposable c d Tensor Scale x y y where
    {-# INLINE unsafeMatrixMultiply #-}
    unsafeMatrixMultiply f g = breakPoint $ S.head (coordinates g) .> f


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

instance (Manifold x, Manifold y) => Bilinear c Tensor x y where
    {-# INLINE (>.<) #-}
    (>.<) (Point pxs) (Point qxs) = fromMatrix $ pxs `S.outerProduct` qxs
    {-# INLINE (>$<) #-}
    (>$<) ps qs = fromMatrix . S.averageOuterProduct $ zip (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose (Point xs) = fromMatrix . S.transpose $ G.Matrix xs
    {-# INLINE toTensor #-}
    toTensor = id
    {-# INLINE fromTensor #-}
    fromTensor = id

instance Manifold x => Bilinear c Symmetric x x where
    {-# INLINE transpose #-}
    transpose = id
    {-# INLINE toTensor #-}
    toTensor (Point trng) = fromMatrix $ S.fromLowerTriangular trng
    {-# INLINE fromTensor #-}
    fromTensor = Point . S.lowerTriangular . toMatrix

instance Manifold x => Bilinear c Diagonal x x where
    {-# INLINE (>.<) #-}
    (>.<) (Point xs) (Point ys) = Point $ xs * ys
    {-# INLINE (>$<) #-}
    (>$<) ps qs = Point . average $ zipWith (*) (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose = id
    {-# INLINE toTensor #-}
    toTensor (Point diag) = fromMatrix $ S.diagonalMatrix diag
    {-# INLINE fromTensor #-}
    fromTensor = Point . S.takeDiagonal . toMatrix



instance Manifold x => Bilinear c Scale x x where
    {-# INLINE (>.<) #-}
    (>.<) (Point x) (Point y) = singleton . S.average $ x * y
    {-# INLINE (>$<) #-}
    (>$<) ps qs = singleton . average
        $ zipWith (\x y -> S.average $ x * y) (coordinates <$> ps) (coordinates <$> qs)
    {-# INLINE transpose #-}
    transpose = id
    {-# INLINE toTensor #-}
    toTensor (Point scl) =
         fromMatrix . S.diagonalMatrix $ S.scale (S.head scl) 1
    {-# INLINE fromTensor #-}
    fromTensor = singleton . S.average . S.takeDiagonal . toMatrix


-- Square --


instance Manifold x => Square c Tensor x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse p = fromMatrix . S.pseudoInverse $ toMatrix p
    {-# INLINE matrixRoot #-}
    matrixRoot p = fromMatrix . S.matrixRoot $ toMatrix p
    {-# INLINE determinant #-}
    determinant = S.determinant . toMatrix
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant tns =
        let (imtx,lndet,sgn) = S.inverseLogDeterminant $ toMatrix tns
         in (fromMatrix imtx, lndet, sgn)

instance Manifold x => Square c Symmetric x where
    -- | The inverse of a tensor.

instance Manifold x => Square c Diagonal x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse (Point diag) = Point $ recip diag
    {-# INLINE matrixRoot #-}
    matrixRoot (Point diag) = Point $ sqrt diag
    {-# INLINE determinant #-}
    determinant (Point diag) = S.product diag
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant sqr =
        let diag = coordinates sqr
            prd = S.product diag
            lndet = log $ abs prd
         in (inverse sqr, lndet, signum prd)

instance Manifold x => Square c Scale x where
    -- | The inverse of a tensor.
    {-# INLINE inverse #-}
    inverse (Point scl) = Point $ recip scl
    {-# INLINE matrixRoot #-}
    matrixRoot (Point scl) = Point $ sqrt scl
    {-# INLINE determinant #-}
    determinant (Point scl) = S.head scl ^ (natValInt $ Proxy @(Dimension x))
    {-# INLINE inverseLogDeterminant #-}
    inverseLogDeterminant sqr =
        let scl = S.head $ coordinates sqr
            lndet = (fromIntegral . natVal $ Proxy @(Dimension x)) * log (abs scl)
         in (inverse sqr, lndet, signum scl)


--- Change of Basis ---


--instance (Manifold x, Manifold y) => Form Tensor x y where
--    changeOfBasis f g =
--        let fmtx = toMatrix $ toTensor f
--            gmtx = toMatrix $ toTensor g
--         in fromMatrix . S.matrixMatrixMultiply (S.transpose fmtx) $ S.matrixMatrixMultiply gmtx fmtx
--
--
----instance Manifold x => Form Symmetric x y where
--
--instance (Manifold x, Manifold y) => Form Diagonal x y where
--    {-# INLINE changeOfBasis #-}
--    changeOfBasis f g =
--        let fmtx = toMatrix f
--            gvec = coordinates g
--         in fromMatrix . S.matrixMatrixMultiply (S.transpose fmtx)
--             $ S.diagonalMatrixMatrixMultiply gvec fmtx
--
--instance (Manifold x, Manifold y) => Form Scale x y where
--    {-# INLINE changeOfBasis #-}
--    changeOfBasis f g =
--        let fmtx = toMatrix f
--            gscl = S.head $ coordinates g
--         in fromMatrix . S.matrixMatrixMultiply (S.transpose fmtx)
--             $ S.withMatrix (S.scale gscl) fmtx


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
         in yz >+> (yx >.> x)
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

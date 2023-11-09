{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | This module provides tools for working with linear and affine
transformations.
-}
module Goal.Geometry.Map.Linear (
    -- * Linear Maps
    Linear,
    KnownLinear (..),
    Tensor,
    Symmetric,
    PositiveDefinite,
    Diagonal,
    Scale,
    Identity,

    -- * Construction
    (>.<),
    (>$<),
    identity,
    toTensor,
    fromTensor,
    toRows,
    toColumns,
    fromRows,
    fromColumns,

    -- * Operations
    (<.<),
    (<$<),
    transpose,
    inverse,
    choleskyDecomposition,
    determinant,
    inverseLogDeterminant,

    -- * Composition
    (<#>),
    dualComposition,
    changeOfBasis,
    schurComplement,

    -- * Affine Maps
    Affine (..),
    LinearSubspace (..),
    (>.+>),
    (>$+>),
    type (<*),
) where

--- Imports ---

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Vector

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Linear Maps ---

-- | A 'Manifold' of 'Linear' operators.
data Linear (t :: L.LinearRep) y x

-- | 'Manifold' of 'Tensor's.
type Tensor y x = Linear L.Full y x

-- | 'Manifold' of 'Symmetric' tensors.
type Symmetric x = Linear L.Symmetric x x

-- | 'Manifold' of 'PositiveDefinite' tensors.
type PositiveDefinite x = Linear L.PositiveDefinite x x

-- | 'Manifold' of 'Diagonal' tensors.
type Diagonal x = Linear L.Diagonal x x

-- | 'Manifold' of 'Scale' transformations.
type Scale x = Linear L.Scale x x

-- | Zero dimensional 'Manifold' of the identity transformation.
type Identity x = Linear L.Identity x x

-- | A class for converting from a 'Linear' 'Manifold' to its non-geometric representation.
class
    ( Manifold x
    , Manifold y
    , KnownNat (L.LinearParameters t (Dimension y) (Dimension x))
    , L.LinearConstruct t (Dimension y) (Dimension x)
    ) =>
    KnownLinear (t :: L.LinearRep) y x
    where
    useLinear :: c # Linear t y x -> L.Linear t (Dimension y) (Dimension x)

--- Construction ---

-- | Tensor outer product.
(>.<) :: forall c t y x. (KnownLinear t y x) => c # y -> c # x -> c # Linear t y x
{-# INLINE (>.<) #-}
(>.<) y x =
    let f :: L.Linear t (Dimension y) (Dimension x)
        f = L.outerProduct (coordinates y) (coordinates x)
     in Point $ L.toVector f

-- | Average of tensor outer products.
(>$<) :: forall c t y x. (KnownLinear t y x) => [c # y] -> [c # x] -> c # Linear t y x
{-# INLINE (>$<) #-}
(>$<) ys xs =
    let f :: L.Linear t (Dimension y) (Dimension x)
        f = L.averageOuterProduct (coordinates <$> ys) (coordinates <$> xs)
     in Point $ L.toVector f

identity :: forall c t x. (KnownLinear t x x) => c # Linear t x x
identity =
    let f :: L.Linear t (Dimension x) (Dimension x)
        f = L.identity
     in Point $ L.toVector f

toTensor :: (KnownLinear t y x) => c # Linear t y x -> c # Tensor y x
toTensor = Point . G.toVector . L.toMatrix . useLinear

fromTensor :: forall c t y x. (KnownLinear t y x) => c # Tensor y x -> c # Linear t y x
fromTensor f =
    let f' :: L.Linear t (Dimension y) (Dimension x)
        f' = L.fromMatrix . L.toMatrix $ useLinear f
     in Point $ L.toVector f'

-- | Converts a point on a 'Tensor' manifold into a a vector of rows.
toRows :: (Manifold x, Manifold y) => c # Tensor y x -> S.Vector (Dimension y) (c # x)
{-# INLINE toRows #-}
toRows tns = S.map Point . S.toRows . L.toMatrix $ useLinear tns

-- | Converts a point on a 'Tensor' manifold into a a vector of rows.
toColumns :: (Manifold x, Manifold y) => c # Tensor y x -> S.Vector (Dimension x) (c # y)
{-# INLINE toColumns #-}
toColumns tns = S.map Point . S.toColumns . L.toMatrix $ useLinear tns

-- | Converts a vector of rows into a 'Tensor'.
fromRows :: forall c x y. (Manifold x, Manifold y) => S.Vector (Dimension y) (c # x) -> c # Tensor y x
{-# INLINE fromRows #-}
fromRows rws =
    let f :: L.Linear L.Full (Dimension y) (Dimension x)
        f = L.fromMatrix . S.fromRows $ S.map coordinates rws
     in Point $ L.toVector f

-- | Converts a vector of rows into a 'Tensor'.
fromColumns :: forall c x y. (Manifold x, Manifold y) => S.Vector (Dimension x) (c # y) -> c # Tensor y x
{-# INLINE fromColumns #-}
fromColumns rws =
    let f :: L.Linear L.Full (Dimension y) (Dimension x)
        f = L.fromMatrix . S.fromColumns $ S.map coordinates rws
     in Point $ L.toVector f

--- Operations ---

transpose ::
    (KnownLinear t y x) =>
    c # Linear t y x ->
    c # Linear t x y
{-# INLINE transpose #-}
transpose = Point . L.toVector . L.transpose . useLinear

-- | Transposed application.
(<.<) :: (KnownLinear t y x, KnownLinear t x y) => c #* y -> c # Linear t y x -> c # x
{-# INLINE (<.<) #-}
(<.<) dx f = transpose f >.> dx

-- | Mapped transposed application.
(<$<) :: (KnownLinear t y x, KnownLinear t x y) => [c #* y] -> c # Linear t y x -> [c # x]
{-# INLINE (<$<) #-}
(<$<) dx f = transpose f >$> dx

-- | 'Linear' operator inversions.
inverse :: (KnownLinear t x x) => c # Linear t x x -> c #* Linear t x x
{-# INLINE inverse #-}
inverse = Point . L.toVector . L.inverse . useLinear

-- | Cholesky decomposition of a 'Symmetric' operator (does not check positive definiteness).
choleskyDecomposition ::
    (KnownLinear L.PositiveDefinite x x) => c # PositiveDefinite x -> c # Tensor x x
{-# INLINE choleskyDecomposition #-}
choleskyDecomposition =
    Point . L.toVector . L.transpose . L.choleskyDecomposition . useLinear

determinant :: (KnownLinear t x x) => c # Linear t x x -> Double
{-# INLINE determinant #-}
determinant = L.determinant . useLinear

inverseLogDeterminant :: (KnownLinear t x x) => c # Linear t x x -> (c #* Linear t x x, Double, Double)
{-# INLINE inverseLogDeterminant #-}
inverseLogDeterminant tns =
    let (inv, lndt, sgn) = L.inverseLogDeterminant $ useLinear tns
     in (Point $ L.toVector inv, lndt, sgn)

--- Composition ---

-- | Composition for boring linear operators.
(<#>) ::
    (KnownLinear t y x, KnownLinear s z y, c ~ Dual c) =>
    c # Linear s z y ->
    c # Linear t y x ->
    c # Linear (L.LinearCompose s t) z x
{-# INLINE (<#>) #-}
(<#>) = unsafeLinearMultiply

-- | Type safe composition of three linear operators that respects duality of coordinate systems.
dualComposition ::
    ( KnownLinear t x w
    , KnownLinear s y x
    , KnownLinear r z y
    , KnownLinear (L.LinearCompose s t) y w
    ) =>
    c # Linear r z y ->
    c #* Linear s y x ->
    c # Linear t x w ->
    c # Linear (L.LinearCompose r (L.LinearCompose s t)) z w
{-# INLINE dualComposition #-}
dualComposition h g f = unsafeLinearMultiply h (unsafeLinearMultiply g f)

-- | Linear change of basis.
changeOfBasis ::
    ( KnownLinear t x y
    , KnownLinear (L.LinearCompose s t) x y
    , KnownLinear t y x
    , KnownLinear s x x
    ) =>
    c # Linear t x y ->
    c #* Linear s x x ->
    c # Linear (L.LinearCompose t (L.LinearCompose s t)) y y
{-# INLINE changeOfBasis #-}
changeOfBasis f g = dualComposition (transpose f) g f

{- | For a block matrix [[A,B],[C,D]], computes the Schur complement of A -- the first argument is the inverse of A, and the subsequent arguments are B, C, and D. The type of the returned matrix is based on the linear type of D,
so that only the necessary parts of the Schur complement are actually computed.
-}
schurComplement ::
    ( KnownLinear f x x
    , KnownLinear g y y
    ) =>
    c #* Linear f x x ->
    c # Tensor x y ->
    c # Tensor y x ->
    c # Linear g y y ->
    c # Linear g y y
schurComplement ainv b c d =
    case useLinear d of
        L.DiagonalLinear _ ->
            let diag = S.zipWith (<.>) (toRows c) (toColumns $ unsafeLinearMultiply ainv b)
             in d - Point diag
        L.ScaleLinear _ ->
            let s = S.singleton . S.average $ S.zipWith (<.>) (toRows c) (toColumns $ unsafeLinearMultiply ainv b)
             in d - Point s
        _ -> d - fromTensor (toTensor $ dualComposition c ainv b)

--- Affine Maps ---

{- | An 'Affine' 'Manifold' represents linear transformations followed by a
translation. The 'First' component is the translation, and the 'Second'
component is the linear transformation.
-}
newtype Affine t y0 y x = Affine (y, Linear t y0 x)

deriving instance (Manifold y, Manifold (Linear t y0 x)) => Manifold (Affine t y0 y x)
deriving instance (Manifold y, Manifold (Linear t y0 x)) => Product (Affine t y0 y x)

-- | Infix synonym for simple 'Affine' transformations.
type y <* x = Affine L.Full y y x

infixr 6 <*

{- | The 'LinearSubspace' class is used to define operations between a larger space 'x' and a subspace 'x0'.
This is based on a non-rigorous interpretation of linear subspaces.
-}
class (Manifold x, Manifold x0) => LinearSubspace x x0 where
    -- | Translates the the first argument by the second argument.
    (>+>) :: c # x -> c # x0 -> c # x

    -- | Returns the subset of the parameters of the given 'Point' that are
    -- translated in this instance.
    projection :: c # x -> c # x0

-- | Operator that applies a 'Map' to a subset of an input's parameters.
(>.+>) :: (Map c f y x0, LinearSubspace x x0) => c # f y x0 -> c #* x -> c # y
(>.+>) f w = f >.> projection w

-- | Operator that maps a 'Map' over a subset of the parameters of a list of inputs.
(>$+>) :: (Map c f y x0, LinearSubspace x x0) => c # f y x0 -> [c #* x] -> [c # y]
(>$+>) f w = f >$> (projection <$> w)

--- Internal ---

unsafeLinearMultiply ::
    (KnownLinear t y x, KnownLinear s z y) =>
    d # Linear s z y ->
    c # Linear t y x ->
    e # Linear (L.LinearCompose s t) z x
{-# INLINE unsafeLinearMultiply #-}
unsafeLinearMultiply g f =
    Point . L.toVector $ L.linearCompose (useLinear g) (useLinear f)

--- Instances ---

instance (KnownLinear t y x, Manifold x, Manifold y) => Manifold (Linear t y x) where
    type Dimension (Linear t y x) = L.LinearParameters t (Dimension y) (Dimension x)

instance (Manifold x, Manifold y) => KnownLinear L.Full y x where
    useLinear (Point xs) = L.FullLinear xs

instance (Manifold x) => KnownLinear L.Symmetric x x where
    useLinear (Point xs) = L.SymmetricLinear xs

instance (Manifold x) => KnownLinear L.PositiveDefinite x x where
    useLinear (Point xs) = L.PositiveDefiniteLinear xs

instance (Manifold x) => KnownLinear L.Diagonal x x where
    useLinear (Point xs) = L.DiagonalLinear xs

instance (Manifold x) => KnownLinear L.Scale x x where
    useLinear (Point xs) = L.ScaleLinear $ S.head xs

instance (Manifold x) => KnownLinear L.Identity x x where
    useLinear (Point _) = L.IdentityLinear

instance (Manifold z) => LinearSubspace z z where
    (>+>) z1 z2 = z1 + z2
    projection = id

instance (Manifold z, Manifold y) => LinearSubspace (y, z) y where
    (>+>) yz y' =
        let (y, z) = split yz
         in join (y + y') z
    projection = fst . split

instance (LinearSubspace z y, KnownLinear t y x) => Map c (Affine t y) z x where
    {-# INLINE (>.>) #-}
    (>.>) fyzx x =
        let (yz, yx) = split fyzx
         in yz >+> (yx >.> x)
    (>$>) fyzx xs =
        let (yz, yx) = split fyzx
         in (yz >+>) <$> yx >$> xs

instance (KnownLinear t y x) => Map c (Linear t) y x where
    {-# INLINE (>.>) #-}
    (>.>) pq (Point xs) = Point $ L.linearVectorMultiply (useLinear pq) xs
    {-# INLINE (>$>) #-}
    (>$>) pq qs = map Point . L.linearMap (useLinear pq) $ coordinates <$> qs

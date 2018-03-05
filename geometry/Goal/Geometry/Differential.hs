{-# LANGUAGE UndecidableInstances #-}

-- | This module provides tools for working with differential and Riemannian
-- geometry.
module Goal.Geometry.Differential (
    -- * Tangent Spaces
      TangentSpace
    , TangentBundle
    -- ** Charts
    , Directional
    , Differential
    -- ** Synonyms
    , TangentVector
    , TangentPair
    , TangentTensor
    , CotangentVector
    , CotangentPair
    , CotangentTensor
    -- ** Manipulation
    , joinTangentPair
    , splitTangentPair
    , projectTangentPair
    , detachTangentVector
    -- ** Replicated Tangent Spaces
    , replicatedJoinTangentPair
    , replicatedSplitTangentPair
    , replicatedJoinTangentSpace
    , replicatedSplitTangentSpace
    -- * Riemannian Manifolds
    , Riemannian (metric, flat, sharp)
    , euclideanDistance
    -- * Differentiation
    , differential
    , differential'
    , hessian
    -- ** Gradient Descent
    , gradientStep
    , gradientStep'
    ) where


--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S

import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Map
import Goal.Geometry.Map.Multilinear

-- Qualified --

import qualified Numeric.AD as D


--- Differentiable Manifolds ---

-- | A tangent bundle on a 'Manifold'.
data TangentBundle c m
-- | A tangent space on a 'Manifold'.
data TangentSpace c m

-- | A 'TangentBundle' or 'TangentSpace' in 'Directional' coordinates represents
-- what is usually called a tangent bundle/space in mathematics; a vector space
-- with a basis given by the derivatives of a set of coordinate functions.
data Directional

-- | A 'TangentBundle' or 'TangentSpace' in 'Differential' coordinates
-- represents a cotangent bundle/space.
data Differential

-- | A synonym for an element of a 'TangentSpace' in 'Directional' coordinates,
-- i.e. a tangent vector. Note that a 'TangentVector' does not include the
-- location of the 'TangentSpace'.
type TangentVector c m x = Point Directional (TangentSpace c m) x

-- | A synonym for an element of a 'TangentSpace' in 'Directional' coordinates,
-- i.e. a cotangent vector.
type CotangentVector c m x = Point Differential (TangentSpace c m) x

-- | A synonym for an element of a 'TangentBundle' in 'Directional' coordinates,
-- i.e. a tangent vector bundled with the location of the tangent space.
type TangentPair c m x = Point Directional (TangentBundle c m) x

-- | A synonym for an element of a 'TangentBundle' in 'Differential' coordinates,
-- i.e. a cotangent vector bundled with the location of the cotangent space.
type CotangentPair c m x = Point Differential (TangentBundle c m) x

-- | A synonym for a 'Tensor' which results from the product of two tangent spaces.
type TangentTensor c m x =
    Point (Function Differential Directional) (Product (TangentSpace c m) (TangentSpace c m)) x

-- | A synonym for a 'Tensor' which results from the product of two cotangent
-- spaces. The 'Riemannian' 'metric' is a form of 'Cotangent' 'Tensor'.
type CotangentTensor c m x
  = Point (Function Directional Differential) (Product (TangentSpace c m) (TangentSpace c m)) x

-- | Computes the differential of a function at a point. This functions returns
-- only the resulting 'CotangentVector', without the corresponding 'Point' where
-- the differential was evaluated.
differential
    :: (Manifold m, Numeric x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point c m x
    -> CotangentVector c m x
{-# INLINE differential #-}
differential f p =
    let p' = fromBPoint . D.grad f $ toBPoint p
     in Point $ coordinates p'

-- | Computes the differential of a function at a point. This functions returns
-- the 'CotangentPair', which includes the 'Point' where the differential was
-- evaluated.
differential'
    :: (Manifold m, Numeric x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point c m x
    -> CotangentPair c m x
{-# INLINE differential' #-}
differential' f p =
    let p' = fromBPoint . D.grad f $ toBPoint p
     in Point $ coordinates p S.++ coordinates p'

-- | Computes the Hessian of a function at a point. This functions returns
-- only the resulting 'CotangentTensor', without the corresponding 'Point' where
-- the Hessian was evaluated.
hessian
    :: (Manifold m, Numeric x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point c m x
    -> CotangentTensor c m x -- ^ The Differential
{-# INLINE hessian #-}
hessian f p =
    let Point xss = fromBPoint $ convert . bCoordinates <$> D.hessian f (toBPoint p)
     in Point . S.toVector $ S.fromRows xss

-- | 'gradientStep' takes a step size, the location of a 'TangentVector', the
-- 'TangentVector' itself, and returns a 'Point' with coordinates that have
-- moved in the direction of the 'TangentVector'.
gradientStep
    :: (Manifold m, Storable x, Num x)
    => x -> Point c m x -> TangentVector c m x -> Point c m x
{-# INLINE gradientStep #-}
gradientStep eps (Point xs) pd =
    Point $ xs + coordinates (eps .> pd)

-- | 'gradientStep' takes a step size and a 'TangentPair', and returns the
-- underlying 'Point' with coordinates shifted in the direction of the
-- 'TangentVector'.
gradientStep'
    :: (Manifold m, Num x, Storable x)
    => x -> TangentPair c m x -> Point c m x
{-# INLINE gradientStep' #-}
gradientStep' eps ppd =
    uncurry (gradientStep eps) $ splitTangentPair ppd

-- | Extract the underlying 'Point' from a 'TangentPair' or 'CotangentPair'.
projectTangentPair
    :: (Manifold m, Storable x)
    => Point d (TangentBundle c m) x
    -> Point c m x
{-# INLINE projectTangentPair #-}
projectTangentPair = fst . splitTangentPair

-- | Detach the 'TangentVector' or 'CotangentVector' from the underlying 'Point'
-- of a pair.
detachTangentVector
    :: (Manifold m, Storable x)
    => Point d (TangentBundle c m) x
    -> Point d (TangentSpace c m) x
{-# INLINE detachTangentVector #-}
detachTangentVector = snd . splitTangentPair

-- | Combine a 'Point' and a 'TangentVector' or 'CotangentVector' into a
-- 'TangentPair' or 'CotangentPair'.
joinTangentPair
    :: (Manifold m, Storable x)
    => Point c m x
    -> Point d (TangentSpace c m) x
    -> Point d (TangentBundle c m) x
{-# INLINE joinTangentPair #-}
joinTangentPair (Point xs) (Point xds) =
    Point $ xs S.++ xds

-- | Split a 'TangentPair' or 'CotangentPair' into a 'Point' and a
-- 'TangentVector' or 'CotangentVector'.
splitTangentPair
    :: (Manifold m, Storable x)
    => Point d (TangentBundle c m) x
    -> (Point c m x, Point d (TangentSpace c m) x)
{-# INLINE splitTangentPair #-}
splitTangentPair (Point xxds) =
    let (x,v) = S.splitAt xxds
     in (Point x, Point v)

-- | Split a 'TangentPair' or 'CotangentPair' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentPair's or 'CotangentPair's.
replicatedSplitTangentPair
    :: (KnownNat k, Manifold m, Storable x)
    => Point d (TangentBundle c (Replicated k m)) x
    -> Vector k (Point d (TangentBundle c m) x)
{-# INLINE replicatedSplitTangentPair #-}
replicatedSplitTangentPair rpv =
    let (rp,rv) = splitTangentPair rpv
        rps = coordinates `S.map` splitReplicated rp
        rvs = coordinates `S.map` replicatedSplitTangentSpace rv
     in Point `S.map` S.zipWith (S.++) rps rvs

-- | Join a 'Vector' of 'TangentPair's or 'CotangentPair's into a 'TangentPair'
-- or 'CotangentPair' on a 'Replicated' 'Manifold'.
replicatedJoinTangentPair
    :: (KnownNat k, Manifold m, Storable x)
    => Vector k (Point d (TangentBundle c m) x)
    -> Point d (TangentBundle c (Replicated k m)) x
{-# INLINE replicatedJoinTangentPair #-}
replicatedJoinTangentPair pvs =
    let p = joinReplicated $ S.map (fst . splitTangentPair) pvs
        v = replicatedJoinTangentSpace $ S.map (snd . splitTangentPair) pvs
     in joinTangentPair p v

-- | Split a 'TangentVector' or 'CotangentVector' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentVector's or 'CotangentVector's.
replicatedSplitTangentSpace
    :: (KnownNat k, Manifold m, Storable x)
    => Point d (TangentSpace c (Replicated k m)) x
    -> Vector k (Point d (TangentSpace c m) x)
{-# INLINE replicatedSplitTangentSpace #-}
replicatedSplitTangentSpace (Point xs) = S.map Point $ S.breakEvery xs

-- | Join a 'Vector' of 'TangentVector's or 'CotangentVector's into a 'TangentVector'
-- or 'CotangentVector' on a 'Replicated' 'Manifold'.
replicatedJoinTangentSpace
    :: (KnownNat k, Manifold m, Storable x)
    => Vector k (Point d (TangentSpace c m) x)
    -> Point d (TangentSpace c (Replicated k m)) x
{-# INLINE replicatedJoinTangentSpace #-}
replicatedJoinTangentSpace ps = Point . S.concat $ S.map coordinates ps

-- | Distance between two 'Point's in 'Euclidean' space.
euclideanDistance
    :: (KnownNat k, Numeric x, Floating x)
    => Point Cartesian (Euclidean k) x
    -> Point Cartesian (Euclidean k) x
    -> x
euclideanDistance (Point xs) (Point ys) = sqrt . S.sum . S.map (^(2 :: Int)) $ xs - ys


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's where associated
-- with each 'Point' in the 'Manifold' is a 'TangentSpace' with a smoothly
-- varying 'CotangentTensor' known as the 'metric'. 'flat' and 'sharp' correspond to applying this 'metric' to elements of the 'TangentBundle' and 'CotangentBundle', respectively.
class Manifold m => Riemannian c m where
    metric :: Numeric x => Point c m x -> CotangentTensor c m x
    flat :: (Numeric x, Field x) => TangentPair c m x -> CotangentPair c m x
    flat pd =
        let (p,v) = splitTangentPair pd
            v' = metric p >.> v
         in joinTangentPair p v'
    sharp :: (Numeric x, Field x) => CotangentPair c m x -> TangentPair c m x
    sharp dp =
        let (p,v) = splitTangentPair dp
            v' = inverse (metric p) >.> v
         in joinTangentPair p v'


--- Instances ---


-- Euclidean --

instance (KnownNat k) => Riemannian Cartesian (Euclidean k) where
    metric _ = fromMatrix S.matrixIdentity
    flat = breakChart
    sharp = breakChart

---- Tangent Spaces --

instance Manifold m => Manifold (TangentSpace c m) where
    type Dimension (TangentSpace c m) = Dimension m

instance Manifold m => Manifold (TangentBundle c m) where
    type Dimension (TangentBundle c m) = 2 * Dimension m

-- Tanget Space Coordinates --

instance Primal Directional where
    type Dual Directional = Differential

instance Primal Differential where
    type Dual Differential = Directional

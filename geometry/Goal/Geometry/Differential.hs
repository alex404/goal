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
import qualified Goal.Core.Vector.Generic as G

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
type TangentVector v c m x = Point v Directional (TangentSpace c m) x

-- | A synonym for an element of a 'TangentSpace' in 'Directional' coordinates,
-- i.e. a cotangent vector.
type CotangentVector v c m x = Point v Differential (TangentSpace c m) x

-- | A synonym for an element of a 'TangentBundle' in 'Directional' coordinates,
-- i.e. a tangent vector bundled with the location of the tangent space.
type TangentPair v c m x = Point v Directional (TangentBundle c m) x

-- | A synonym for an element of a 'TangentBundle' in 'Differential' coordinates,
-- i.e. a cotangent vector bundled with the location of the cotangent space.
type CotangentPair v c m x = Point v Differential (TangentBundle c m) x

-- | A synonym for a 'Tensor' which results from the product of two tangent spaces.
type TangentTensor v c m x =
    Point v (Function Differential Directional) (Product (TangentSpace c m) (TangentSpace c m)) x

-- | A synonym for a 'Tensor' which results from the product of two cotangent
-- spaces. The 'Riemannian' 'metric' is a form of 'Cotangent' 'Tensor'.
type CotangentTensor v c m x
  = Point v (Function Directional Differential) (Product (TangentSpace c m) (TangentSpace c m)) x

-- | Computes the differential of a function at a point. This functions returns
-- only the resulting 'CotangentVector', without the corresponding 'Point' where
-- the differential was evaluated.
differential
    :: (Manifold m, GVector v x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point v c m x
    -> CotangentVector v c m x
{-# INLINE differential #-}
differential f p =
    let p' = convertPoint . D.grad f $ convertPoint p
     in Point $ coordinates p'

-- | Computes the differential of a function at a point. This functions returns
-- the 'CotangentPair', which includes the 'Point' where the differential was
-- evaluated.
differential'
    :: (Manifold m, GVector v x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point v c m x
    -> CotangentPair v c m x
{-# INLINE differential' #-}
differential' f p =
    let p' = convertPoint . D.grad f $ convertPoint p
     in Point $ coordinates p G.++ coordinates p'

-- | Computes the Hessian of a function at a point. This functions returns
-- only the resulting 'CotangentTensor', without the corresponding 'Point' where
-- the Hessian was evaluated.
hessian
    :: (Manifold m, GPoint v c m x, RealFloat x)
    => (forall z. RealFloat z => BPoint c m z -> z)
    -> Point v c m x
    -> CotangentTensor v c m x -- ^ The Differential
{-# INLINE hessian #-}
hessian f p =
    let Point xss = convertPoint $ convert . coordinates <$> D.hessian f (convertPoint p)
     in Point . G.toVector $ G.fromRows xss

-- | 'gradientStep' takes a step size, the location of a 'TangentVector', the
-- 'TangentVector' itself, and returns a 'Point' with coordinates that have
-- moved in the direction of the 'TangentVector'.
gradientStep
    :: (Manifold m, GVector v x, Num x)
    => x -> Point v c m x -> TangentVector v c m x -> Point v c m x
{-# INLINE gradientStep #-}
gradientStep eps (Point xs) pd =
    Point $ xs + coordinates (eps .> pd)

-- | 'gradientStep' takes a step size and a 'TangentPair', and returns the
-- underlying 'Point' with coordinates shifted in the direction of the
-- 'TangentVector'.
gradientStep'
    :: (Manifold m, Num x, GVector v x)
    => x -> TangentPair v c m x -> Point v c m x
{-# INLINE gradientStep' #-}
gradientStep' eps ppd =
    uncurry (gradientStep eps) $ splitTangentPair ppd

-- | Extract the underlying 'Point' from a 'TangentPair' or 'CotangentPair'.
projectTangentPair
    :: (Manifold m, GVector v x)
    => Point v d (TangentBundle c m) x
    -> Point v c m x
{-# INLINE projectTangentPair #-}
projectTangentPair = fst . splitTangentPair

-- | Detach the 'TangentVector' or 'CotangentVector' from the underlying 'Point'
-- of a pair.
detachTangentVector
    :: (Manifold m, GVector v x)
    => Point v d (TangentBundle c m) x
    -> Point v d (TangentSpace c m) x
{-# INLINE detachTangentVector #-}
detachTangentVector = snd . splitTangentPair

-- | Combine a 'Point' and a 'TangentVector' or 'CotangentVector' into a
-- 'TangentPair' or 'CotangentPair'.
joinTangentPair
    :: (Manifold m, GVector v x)
    => Point v c m x
    -> Point v d (TangentSpace c m) x
    -> Point v d (TangentBundle c m) x
{-# INLINE joinTangentPair #-}
joinTangentPair (Point xs) (Point xds) =
    Point $ xs G.++ xds

-- | Split a 'TangentPair' or 'CotangentPair' into a 'Point' and a
-- 'TangentVector' or 'CotangentVector'.
splitTangentPair
    :: (Manifold m, GVector v x)
    => Point v d (TangentBundle c m) x
    -> (Point v c m x, Point v d (TangentSpace c m) x)
{-# INLINE splitTangentPair #-}
splitTangentPair (Point xxds) =
    let (x,v) = G.splitAt xxds
     in (Point x, Point v)

-- | Split a 'TangentPair' or 'CotangentPair' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentPair's or 'CotangentPair's.
replicatedSplitTangentPair
    :: (KnownNat k, Manifold m, GPoint v d (TangentBundle c m) x, GPoint v d (TangentSpace c m) x, GPoint v c m x)
    => Point v d (TangentBundle c (Replicated k m)) x
    -> Vector v k (Point v d (TangentBundle c m) x)
{-# INLINE replicatedSplitTangentPair #-}
replicatedSplitTangentPair rpv =
    let (rp,rv) = splitTangentPair rpv
        rps = coordinates `G.map` splitReplicated rp
        rvs = coordinates `G.map` replicatedSplitTangentSpace rv
     in Point `G.map` G.zipWith (G.++) rps rvs

-- | Join a 'Vector' of 'TangentPair's or 'CotangentPair's into a 'TangentPair'
-- or 'CotangentPair' on a 'Replicated' 'Manifold'.
replicatedJoinTangentPair
    :: (KnownNat k, Manifold m, GPoint v d (TangentBundle c m) x, GPoint v d (TangentSpace c m) x, GPoint v c m x)
    => Vector v k (Point v d (TangentBundle c m) x)
    -> Point v d (TangentBundle c (Replicated k m)) x
{-# INLINE replicatedJoinTangentPair #-}
replicatedJoinTangentPair pvs =
    let p = joinReplicated $ G.map (fst . splitTangentPair) pvs
        v = replicatedJoinTangentSpace $ G.map (snd . splitTangentPair) pvs
     in joinTangentPair p v

-- | Split a 'TangentVector' or 'CotangentVector' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentVector's or 'CotangentVector's.
replicatedSplitTangentSpace
    :: (KnownNat k, Manifold m, GPoint v d (TangentBundle c m) x, GPoint v d (TangentSpace c m) x, GVector v x)
    => Point v d (TangentSpace c (Replicated k m)) x
    -> Vector v k (Point v d (TangentSpace c m) x)
{-# INLINE replicatedSplitTangentSpace #-}
replicatedSplitTangentSpace (Point xs) = G.map Point $ G.breakEvery xs

-- | Join a 'Vector' of 'TangentVector's or 'CotangentVector's into a 'TangentVector'
-- or 'CotangentVector' on a 'Replicated' 'Manifold'.
replicatedJoinTangentSpace
    :: (KnownNat k, Manifold m, GPoint v d (TangentBundle c m) x, GPoint v d (TangentSpace c m) x, GVector v x)
    => Vector v k (Point v d (TangentSpace c m) x)
    -> Point v d (TangentSpace c (Replicated k m)) x
{-# INLINE replicatedJoinTangentSpace #-}
replicatedJoinTangentSpace ps = Point . G.concat $ G.map coordinates ps

-- | Distance between two 'Point's in 'Euclidean' space.
euclideanDistance
    :: (KnownNat k, GVector v x, Floating x)
    => Point v Cartesian (Euclidean k) x
    -> Point v Cartesian (Euclidean k) x
    -> x
euclideanDistance (Point xs) (Point ys) = sqrt . G.sum . G.map (^(2 :: Int)) $ xs - ys


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's where associated
-- with each 'Point' in the 'Manifold' is a 'TangentSpace' with a smoothly
-- varying 'CotangentTensor' known as the 'metric'. 'flat' and 'sharp' correspond to applying this 'metric' to elements of the 'TangentBundle' and 'CotangentBundle', respectively.
class Manifold m => Riemannian v c m x where
    metric :: Point v c m x -> CotangentTensor v c m x
    flat :: TangentPair v c m x -> CotangentPair v c m x
--    flat pd =
--        let (p,v) = splitTangentPair pd
--            v' = metric p >.> v
--         in joinTangentPair p v'
    sharp :: CotangentPair v c m x -> TangentPair v c m x
--    sharp dp =
--        let (p,v) = splitTangentPair dp
--            v' = inverse (metric p) >.> v
--         in joinTangentPair p v'


--- Instances ---


-- Euclidean --

instance (Numeric x, KnownNat k) => Riemannian SVector Cartesian (Euclidean k) x where
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

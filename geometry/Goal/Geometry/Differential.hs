{-# LANGUAGE UndecidableInstances #-}

-- | This module provides tools for modelling the differential and Riemannian
-- geometry of a 'Manifold'.
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
     , primalIsomorphism
     , dualIsomorphism
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
     , Propagate (propagate)
     -- ** Gradient Descent
     , gradientStep
     , gradientStep'
     ) where


--- Imports ---


-- Goal --

import Goal.Core

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
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
type TangentVector c m = Point Directional (TangentSpace c m)

-- | A synonym for an element of a 'TangentSpace' in 'Directional' coordinates,
-- i.e. a cotangent vector.
type CotangentVector c m = Point Differential (TangentSpace c m)

-- | A synonym for an element of a 'TangentBundle' in 'Directional' coordinates,
-- i.e. a tangent vector bundled with the location of the tangent space.
type TangentPair c m = Point Directional (TangentBundle c m)

-- | A synonym for an element of a 'TangentBundle' in 'Differential' coordinates,
-- i.e. a cotangent vector bundled with the location of the cotangent space.
type CotangentPair c m = Point Differential (TangentBundle c m)

-- | A synonym for a 'Tensor' which results from the product of two tangent spaces.
type TangentTensor c m =
    Point (Function Differential Directional) (Tensor (TangentSpace c m) (TangentSpace c m))

-- | A synonym for a 'Tensor' which results from the product of two cotangent
-- spaces. The 'Riemannian' 'metric' is a form of 'Cotangent' 'Tensor'.
type CotangentTensor c m
  = Point (Function Directional Differential) (Tensor (TangentSpace c m) (TangentSpace c m))

-- | Computes the differential of a function of the coordinates at a point. This
-- functions returns only the resulting 'CotangentVector', without the
-- corresponding 'Point' where the differential was evaluated.
differential
    :: Manifold m
    => (forall x. RealFloat x => B.Vector (Dimension m) x -> x)
    -> Point c m
    -> CotangentVector c m
{-# INLINE differential #-}
differential f = Point . G.convert . D.grad f . boxCoordinates

-- | Computes the differential of a function at a point. This functions returns
-- the 'CotangentPair', which includes the 'Point' where the differential was
-- evaluated.
differential'
    :: Manifold m
    => (forall x. RealFloat x => B.Vector (Dimension m) x -> x)
    -> Point c m
    -> CotangentPair c m
{-# INLINE differential' #-}
differential' f p@(Point xs) =
    let dxs = D.grad f $ boxCoordinates p
     in Point $ xs G.++ G.convert dxs

-- | Computes the Hessian of a function at a point. This functions returns
-- only the resulting 'CotangentTensor', without the corresponding 'Point' where
-- the Hessian was evaluated.
hessian
    :: Manifold m
    => (forall x. RealFloat x => B.Vector (Dimension m) x -> x)
    -> Point c m
    -> CotangentTensor c m -- ^ The Differential
{-# INLINE hessian #-}
hessian f p =
    fromMatrix . S.fromRows . G.convert $ G.convert <$> D.hessian f (boxCoordinates p)

-- | A class of functions which can 'propagate' errors. That is, given an error
-- derivative on the output, the input which caused the output, and a
-- function to derive, computes the derivative of the error between the function
-- and target output.
class Map c d f m n => Propagate c d f m n where
    propagate :: [Dual d # m] -- ^ The error derivatives in 'Dual' coordinates
              -> [c # n] -- ^ A vector of inputs
              -> Function c d # f m n -- ^ The function to differentiate
              -> (Function (Dual c) (Dual d) # f m n, [d # m]) -- ^ The derivative, and function output

-- | 'gradientStep' takes a step size, the location of a 'TangentVector', the
-- 'TangentVector' itself, and returns a 'Point' with coordinates that have
-- moved in the direction of the 'TangentVector'.
gradientStep
    :: Manifold m
    => Double -> Point c m -> TangentVector c m -> Point c m
{-# INLINE gradientStep #-}
gradientStep eps (Point xs) pd =
    Point $ xs + coordinates (eps .> pd)

-- | 'gradientStep' takes a step size and a 'TangentPair', and returns the
-- underlying 'Point' with coordinates shifted in the direction of the
-- 'TangentVector'.
gradientStep'
    :: Manifold m
    => Double -> TangentPair c m -> Point c m
{-# INLINE gradientStep' #-}
gradientStep' eps ppd =
    uncurry (gradientStep eps) $ splitTangentPair ppd

-- | Extract the underlying 'Point' from a 'TangentPair' or 'CotangentPair'.
projectTangentPair
    :: Manifold m
    => Point d (TangentBundle c m)
    -> Point c m
{-# INLINE projectTangentPair #-}
projectTangentPair = fst . splitTangentPair

-- | Detach the 'TangentVector' or 'CotangentVector' from the underlying 'Point'
-- of a pair.
detachTangentVector
    :: Manifold m
    => Point d (TangentBundle c m)
    -> Point d (TangentSpace c m)
{-# INLINE detachTangentVector #-}
detachTangentVector = snd . splitTangentPair

-- | Combine a 'Point' and a 'TangentVector' or 'CotangentVector' into a
-- 'TangentPair' or 'CotangentPair'.
joinTangentPair
    :: Manifold m
    => Point c m
    -> Point d (TangentSpace c m)
    -> Point d (TangentBundle c m)
{-# INLINE joinTangentPair #-}
joinTangentPair (Point xs) (Point xds) =
    Point $ xs G.++ xds

-- | Split a 'TangentPair' or 'CotangentPair' into a 'Point' and a
-- 'TangentVector' or 'CotangentVector'.
splitTangentPair
    :: Manifold m
    => Point d (TangentBundle c m)
    -> (Point c m, Point d (TangentSpace c m))
{-# INLINE splitTangentPair #-}
splitTangentPair (Point xxds) =
    let (x,v) = G.splitAt xxds
     in (Point x, Point v)

-- | Split a 'TangentPair' or 'CotangentPair' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentPair's or 'CotangentPair's.
-- NB: Optimize with backpermute
replicatedSplitTangentPair
    :: (KnownNat k, Manifold m)
    => Point d (TangentBundle c (Replicated k m))
    -> S.Vector k (Point d (TangentBundle c m))
{-# INLINE replicatedSplitTangentPair #-}
replicatedSplitTangentPair rpv =
    let (rp,rv) = splitTangentPair rpv
        rps = coordinates `G.map` splitReplicated rp
        rvs = coordinates `G.map` replicatedSplitTangentSpace rv
     in Point `G.map` G.zipWith (G.++) rps rvs

-- | Join a 'Vector' of 'TangentPair's or 'CotangentPair's into a 'TangentPair'
-- or 'CotangentPair' on a 'Replicated' 'Manifold'.
-- NB: Optimize with backpermute
replicatedJoinTangentPair
    :: (KnownNat k, Manifold m)
    => S.Vector k (Point d (TangentBundle c m))
    -> Point d (TangentBundle c (Replicated k m))
{-# INLINE replicatedJoinTangentPair #-}
replicatedJoinTangentPair pvs =
    let p = joinReplicated $ G.map (fst . splitTangentPair) pvs
        v = replicatedJoinTangentSpace $ G.map (snd . splitTangentPair) pvs
     in joinTangentPair p v

-- | Split a 'TangentVector' or 'CotangentVector' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentVector's or 'CotangentVector's.
replicatedSplitTangentSpace
    :: (KnownNat k, Manifold m)
    => Point d (TangentSpace c (Replicated k m))
    -> S.Vector k (Point d (TangentSpace c m))
{-# INLINE replicatedSplitTangentSpace #-}
replicatedSplitTangentSpace (Point xs) = G.map Point $ G.breakEvery xs

-- | Join a 'Vector' of 'TangentVector's or 'CotangentVector's into a 'TangentVector'
-- or 'CotangentVector' on a 'Replicated' 'Manifold'.
replicatedJoinTangentSpace
    :: (KnownNat k, Manifold m)
    => S.Vector k (Point d (TangentSpace c m))
    -> Point d (TangentSpace c (Replicated k m))
{-# INLINE replicatedJoinTangentSpace #-}
replicatedJoinTangentSpace ps = Point . G.concat $ G.map coordinates ps

-- | Distance between two 'Point's based on the 'Euclidean' metric.
euclideanDistance
    :: Manifold m
    => c # m
    -> c # m
    -> Double
euclideanDistance (Point xs) (Point ys) = S.l2Norm xs ys

-- | The 'Dual' space of a 'Manifold' is often isomorphic to its cotangent space, and we often wish to treat the former as the latter.
primalIsomorphism :: Point c m -> CotangentVector (Dual c) m
{-# INLINE primalIsomorphism #-}
primalIsomorphism (Point xs) = Point xs

-- | The inverse of 'primalIsomorphism'.
dualIsomorphism :: CotangentVector c m -> Point (Dual c) m
{-# INLINE dualIsomorphism #-}
dualIsomorphism (Point xs) =  Point xs


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's where associated
-- with each 'Point' in the 'Manifold' is a 'TangentSpace' with a smoothly
-- varying 'CotangentTensor' known as the 'metric'. 'flat' and 'sharp' correspond to applying this 'metric' to elements of the 'TangentBundle' and 'CotangentBundle', respectively.
class Manifold m => Riemannian c m where
    metric :: Point c m -> CotangentTensor c m
    flat :: TangentPair c m -> CotangentPair c m
    {-# INLINE flat #-}
    flat pd =
        let (p,v) = splitTangentPair pd
            v' = metric p >.> v
         in joinTangentPair p v'
    sharp :: CotangentPair c m -> TangentPair c m
    {-# INLINE sharp #-}
    sharp dp =
        let (p,v) = splitTangentPair dp
            v' = inverse (metric p) >.> v
         in joinTangentPair p v'


--- Instances ---


-- Euclidean --

instance KnownNat k => Riemannian Cartesian (Euclidean k) where
    metric _ = fromMatrix S.matrixIdentity
    flat = breakPoint
    sharp = breakPoint

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

-- Replicated Riemannian Manifolds --

instance {-# OVERLAPPABLE #-} (Riemannian c m, KnownNat k) => Riemannian c (Replicated k m) where
    metric = error "Do not call metric on a replicated manifold"
    {-# INLINE flat #-}
    flat = replicatedJoinTangentPair . S.map flat . replicatedSplitTangentPair
    {-# INLINE sharp #-}
    sharp = replicatedJoinTangentPair . S.map sharp . replicatedSplitTangentPair

-- Backprop --

instance Map c d Tensor m n => Propagate c d Tensor m n where
    {-# INLINE propagate #-}
    propagate dps0 qs0 pq =
        let foldfun (dps,qs) (k,dmtx) = (k+1,(dps >.< qs) <+> dmtx)
         in (uncurry (/>) . foldr foldfun (0,zero) $ zip dps0 qs0, pq >$> qs0)

instance (Map c d (Affine f) m n, Propagate c d f m n) => Propagate c d (Affine f) m n where
    {-# INLINE propagate #-}
    propagate dps qs pq =
        let (p,pq') = splitAffine pq
            (dpq',ps') = propagate dps qs pq'
         in (joinAffine (averagePoint dps) dpq', (p <+>) <$> ps')



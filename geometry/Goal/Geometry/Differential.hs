{-# LANGUAGE
    RankNTypes,
    TypeOperators,
    DataKinds,
    MultiParamTypeClasses,
    TypeFamilies,
    FlexibleInstances,
    UndecidableInstances
    #-}

-- | Tools for modelling the differential and Riemannian geometry of a
-- 'Manifold'.
module Goal.Geometry.Differential
    ( -- * Tangent Spaces
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
    , pairTangentFunction
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
    , hessian
    , Propagate (propagate)
    -- ** Gradient Descent
    , vanillaGradient
    , vanillaGradient'
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
data TangentBundle c x
-- | A tangent space on a 'Manifold'.
data TangentSpace c x

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
type TangentVector c x = Point Directional (TangentSpace c x)

-- | A synonym for an element of a 'TangentSpace' in 'Directional' coordinates,
-- i.e. a cotangent vector.
type CotangentVector c x = Point Differential (TangentSpace c x)

-- | A synonym for an element of a 'TangentBundle' in 'Directional' coordinates,
-- i.e. a tangent vector bundled with the location of the tangent space.
type TangentPair c x = Point Directional (TangentBundle c x)

-- | A synonym for an element of a 'TangentBundle' in 'Differential' coordinates,
-- i.e. a cotangent vector bundled with the location of the cotangent space.
type CotangentPair c x = Point Differential (TangentBundle c x)

-- | A synonym for a 'Tensor' which results from the product of two tangent spaces.
type TangentTensor c x =
    Point (Function Differential Directional) (Tensor (TangentSpace c x) (TangentSpace c x))

-- | A synonym for a 'Tensor' which results from the product of two cotangent
-- spaces. The 'Riemannian' 'metric' is a form of 'Cotangent' 'Tensor'.
type CotangentTensor c x
  = Point (Function Directional Differential) (Tensor (TangentSpace c x) (TangentSpace c x))

-- | Computes the differential of a function of the coordinates at a point. This
-- functions returns only the resulting 'CotangentVector', without the
-- corresponding 'Point' where the differential was evaluated.
differential
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> Point c x
    -> CotangentVector c x
{-# INLINE differential #-}
differential f = Point . G.convert . D.grad f . boxCoordinates

-- | Computes the Hessian of a function at a point. This functions returns
-- only the resulting 'CotangentTensor', without the corresponding 'Point' where
-- the Hessian was evaluated.
hessian
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> Point c x
    -> CotangentTensor c x -- ^ The Differential
{-# INLINE hessian #-}
hessian f p =
    fromMatrix . S.fromRows . G.convert $ G.convert <$> D.hessian f (boxCoordinates p)

-- | A class of functions which can 'propagate' errors. That is, given an error
-- derivative on the output, the input which caused the output, and a
-- function to derive, computes the derivative of the error between the function
-- and target output.
class Map c d f y x => Propagate c d f y x where
    propagate :: [Dual d # y] -- ^ The error derivatives in 'Dual' coordinates
              -> [c # x] -- ^ A vector of inputs
              -> Function c d # f y x -- ^ The function to differentiate
              -> (Function (Dual c) (Dual d) # f y x, [d # y]) -- ^ The derivative, and function output

-- | 'gradientStep' takes a step size, the location of a 'TangentVector', the
-- 'TangentVector' itself, and returns a 'Point' with coordinates that have
-- moved in the direction of the 'TangentVector'.
gradientStep
    :: Manifold x
    => Double -> Point c x -> TangentVector c x -> Point c x
{-# INLINE gradientStep #-}
gradientStep eps (Point xs) pd =
    Point $ xs + coordinates (eps .> pd)

-- | 'gradientStep' takes a step size and a 'TangentPair', and returns the
-- underlying 'Point' with coordinates shifted in the direction of the
-- 'TangentVector'.
gradientStep'
    :: Manifold x
    => Double -> TangentPair c x -> Point c x
{-# INLINE gradientStep' #-}
gradientStep' eps ppd =
    uncurry (gradientStep eps) $ splitTangentPair ppd

-- | Extract the underlying 'Point' from a 'TangentPair' or 'CotangentPair'.
projectTangentPair
    :: Manifold x
    => Point d (TangentBundle c x)
    -> Point c x
{-# INLINE projectTangentPair #-}
projectTangentPair = fst . splitTangentPair

-- | Detach the 'TangentVector' or 'CotangentVector' from the underlying 'Point'
-- of a pair.
detachTangentVector
    :: Manifold x
    => Point d (TangentBundle c x)
    -> Point d (TangentSpace c x)
{-# INLINE detachTangentVector #-}
detachTangentVector = snd . splitTangentPair

-- | Combine a 'Point' and a 'TangentVector' or 'CotangentVector' into a
-- 'TangentPair' or 'CotangentPair'.
joinTangentPair
    :: Manifold x
    => Point c x
    -> Point d (TangentSpace c x)
    -> Point d (TangentBundle c x)
{-# INLINE joinTangentPair #-}
joinTangentPair (Point xs) (Point xds) =
    Point $ xs G.++ xds

-- | Turns a function from a Point to a Tangent into a function from a Point to
-- a TangentPair.
pairTangentFunction
    :: Manifold x
    => (c # x -> d # TangentSpace c x)
    -> (c # x -> d # TangentBundle c x)
{-# INLINE pairTangentFunction #-}
pairTangentFunction f p = joinTangentPair p (f p)

-- | Split a 'TangentPair' or 'CotangentPair' into a 'Point' and a
-- 'TangentVector' or 'CotangentVector'.
splitTangentPair
    :: Manifold x
    => Point d (TangentBundle c x)
    -> (Point c x, Point d (TangentSpace c x))
{-# INLINE splitTangentPair #-}
splitTangentPair (Point xxds) =
    let (x,v) = G.splitAt xxds
     in (Point x, Point v)

-- | Split a 'TangentPair' or 'CotangentPair' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentPair's or 'CotangentPair's.
-- NB: Optimize with backpermute
replicatedSplitTangentPair
    :: (KnownNat k, Manifold x)
    => Point d (TangentBundle c (Replicated k x))
    -> S.Vector k (Point d (TangentBundle c x))
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
    :: (KnownNat k, Manifold x)
    => S.Vector k (Point d (TangentBundle c x))
    -> Point d (TangentBundle c (Replicated k x))
{-# INLINE replicatedJoinTangentPair #-}
replicatedJoinTangentPair pvs =
    let p = joinReplicated $ G.map (fst . splitTangentPair) pvs
        v = replicatedJoinTangentSpace $ G.map (snd . splitTangentPair) pvs
     in joinTangentPair p v

-- | Split a 'TangentVector' or 'CotangentVector' on a 'Replicated' 'Manifold' into
-- a 'Vector' of 'TangentVector's or 'CotangentVector's.
replicatedSplitTangentSpace
    :: (KnownNat k, Manifold x)
    => Point d (TangentSpace c (Replicated k x))
    -> S.Vector k (Point d (TangentSpace c x))
{-# INLINE replicatedSplitTangentSpace #-}
replicatedSplitTangentSpace (Point xs) = G.map Point $ G.breakEvery xs

-- | Join a 'Vector' of 'TangentVector's or 'CotangentVector's into a 'TangentVector'
-- or 'CotangentVector' on a 'Replicated' 'Manifold'.
replicatedJoinTangentSpace
    :: (KnownNat k, Manifold x)
    => S.Vector k (Point d (TangentSpace c x))
    -> Point d (TangentSpace c (Replicated k x))
{-# INLINE replicatedJoinTangentSpace #-}
replicatedJoinTangentSpace ps = Point . G.concat $ G.map coordinates ps

-- | Distance between two 'Point's based on the 'Euclidean' metric.
euclideanDistance
    :: Manifold x
    => c # x
    -> c # x
    -> Double
euclideanDistance (Point xs) (Point ys) = S.l2Norm xs ys

-- | Ignore the Riemannian metric.
vanillaGradient
    :: Manifold x
    => CotangentVector c x
    -> TangentVector c x
{-# INLINE vanillaGradient #-}
vanillaGradient = breakPoint

-- | Ignore the Riemannian metric.
vanillaGradient'
    :: Manifold x
    => CotangentPair c x
    -> TangentPair c x
{-# INLINE vanillaGradient' #-}
vanillaGradient' = breakPoint


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's where associated
-- with each 'Point' in the 'Manifold' is a 'TangentSpace' with a smoothly
-- varying 'CotangentTensor' known as the 'metric'. 'flat' and 'sharp' correspond to applying this 'metric' to elements of the 'TangentBundle' and 'CotangentBundle', respectively.
class Manifold x => Riemannian c x where
    metric :: Point c x -> CotangentTensor c x
    flat :: TangentPair c x -> CotangentPair c x
    {-# INLINE flat #-}
    flat pd =
        let (p,v) = splitTangentPair pd
            v' = metric p >.> v
         in joinTangentPair p v'
    sharp :: CotangentPair c x -> TangentPair c x
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

instance Manifold x => Manifold (TangentSpace c x) where
    type Dimension (TangentSpace c x) = Dimension x

instance Manifold x => Manifold (TangentBundle c x) where
    type Dimension (TangentBundle c x) = 2 * Dimension x

-- Tanget Space Coordinates --

instance Primal Directional where
    type Dual Directional = Differential

instance Primal Differential where
    type Dual Differential = Directional

-- Replicated Riemannian Manifolds --

instance {-# OVERLAPPABLE #-} (Riemannian c x, KnownNat k) => Riemannian c (Replicated k x) where
    metric = error "Do not call metric on a replicated manifold"
    {-# INLINE flat #-}
    flat = replicatedJoinTangentPair . S.map flat . replicatedSplitTangentPair
    {-# INLINE sharp #-}
    sharp = replicatedJoinTangentPair . S.map sharp . replicatedSplitTangentPair

-- Backprop --

instance Map c d Tensor y x => Propagate c d Tensor y x where
    {-# INLINE propagate #-}
    propagate dps0 qs0 pq =
        let foldfun (dps,qs) (k,dmtx) = (k+1,(dps >.< qs) <+> dmtx)
         in (uncurry (/>) . foldr foldfun (0,zero) $ zip dps0 qs0, pq >$> qs0)

instance (Map c d (Affine f) y x, Propagate c d f y x) => Propagate c d (Affine f) y x where
    {-# INLINE propagate #-}
    propagate dps qs pq =
        let (p,pq') = splitAffine pq
            (dpq',ps') = propagate dps qs pq'
         in (joinAffine (averagePoint dps) dpq', (p <+>) <$> ps')



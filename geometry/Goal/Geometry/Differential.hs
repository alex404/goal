{-# LANGUAGE TypeApplications,UndecidableInstances,UndecidableSuperClasses #-}

-- | Tools for modelling the differential and Riemannian geometry of a
-- 'Manifold'.
module Goal.Geometry.Differential
    ( -- * Riemannian Manifolds
      Riemannian (metric, flat, sharp)
    , euclideanDistance
    -- * Backpropagation
    , Propagate (propagate)
    -- * Legendre Manifolds
    , Legendre (PotentialCoordinates,potential)
    , DuallyFlat (dualPotential)
    , canonicalDivergence
    -- * Automatic Differentiation
    , differential
    , hessian
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


-- | Computes the differential of a function of the coordinates at a point using
-- automatic differentiation.
differential
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> c # x
    -> c #* x
{-# INLINE differential #-}
differential f = Point . G.convert . D.grad f . boxCoordinates

-- | Computes the Hessian of a function at a point with automatic differentiation.
hessian
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> c # x
    -> c #*> Tensor x x -- ^ The Hessian
{-# INLINE hessian #-}
hessian f p =
    fromMatrix . S.fromRows . G.convert $ G.convert <$> D.hessian f (boxCoordinates p)

-- | A class of 'Map's which can 'propagate' errors. That is, given an error
-- derivative on the output, the input which caused the output, and a
-- 'Map' to derive, return the derivative of the error with respect to the
-- parameters of the 'Map', as well as the output of the 'Map'.
class Map c d f y x => Propagate c d f y x where
    propagate :: [d #* y] -- ^ The error differential
              -> [c # x] -- ^ A vector of inputs
              -> Function c d # f y x -- ^ The function to differentiate
              -> (Function c d #* f y x, [d # y]) -- ^ The derivative, and function output

-- | Distance between two 'Point's based on the 'Euclidean' metric (l2 distance).
euclideanDistance
    :: Manifold x => c # x -> c # x -> Double
{-# INLINE euclideanDistance #-}
euclideanDistance xs ys = S.l2Norm (coordinates $ xs - ys)


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's associated with a
-- smoothly varying 'Tensor' known as the Riemannian 'metric'. 'flat' and
-- 'sharp' correspond to applying this 'metric' to elements of the 'Primal' and
-- 'Dual' spaces, respectively.
class (Primal c, Manifold x) => Riemannian c x where
    metric :: c # x -> c #*> Tensor x x
    flat :: c # x -> c # x -> c #* x
    {-# INLINE flat #-}
    flat p v = metric p >.> v
    sharp :: c # x -> c #* x -> c # x
    {-# INLINE sharp #-}
    sharp p v = inverse (metric p) >.> v


--- Dually Flat Manifolds ---


-- | Although convex analysis is usually developed seperately from differential
-- geometry, it arises naturally out of the theory of dually flat 'Manifold's (<https://books.google.com/books?hl=en&lr=&id=vc2FWSo7wLUC&oi=fnd&pg=PR7&dq=methods+of+information+geometry&ots=4HsxHD_5KY&sig=gURe0tA3IEO-z-Cht_2TNsjjOG8#v=onepage&q=methods%20of%20information%20geometry&f=false Amari and Nagaoka, 2000>).
--
-- A 'Manifold' is 'Legendre' if it is associated with a particular convex
-- function known as a 'potential'.
class ( Primal (PotentialCoordinates x), Manifold x ) => Legendre x where
    type PotentialCoordinates x :: Type
    potential :: PotentialCoordinates x # x -> Double

-- | A 'Manifold' is 'DuallyFlat' when we can describe the 'dualPotential', which
-- is the convex conjugate of 'potential'.
class Legendre x => DuallyFlat x where
    dualPotential :: PotentialCoordinates x #* x -> Double

-- | Computes the 'canonicalDivergence' between two points. Note that relative
-- to the typical definition of the KL-Divergence/relative entropy, the
-- arguments of this function are flipped.
canonicalDivergence
    :: DuallyFlat x => PotentialCoordinates x # x -> PotentialCoordinates x #* x -> Double
{-# INLINE canonicalDivergence #-}
canonicalDivergence pp dq = potential pp + dualPotential dq - (pp <.> dq)

--- Instances ---


-- Euclidean --

instance KnownNat k => Riemannian Cartesian (Euclidean k) where
    {-# INLINE metric #-}
    metric _ = fromMatrix S.matrixIdentity
    {-# INLINE flat #-}
    flat _ = breakPoint
    {-# INLINE sharp #-}
    sharp _ = breakPoint

-- Replicated Riemannian Manifolds --

--instance {-# OVERLAPPABLE #-} (Riemannian c x, KnownNat k) => Riemannian c (Replicated k x) where
--    metric = error "Do not call metric on a replicated manifold"
--    {-# INLINE flat #-}
--    flat = S.map flat
--    {-# INLINE sharp #-}
--    sharp = S.map sharp

-- Backprop --

instance (Bilinear Tensor y x, Primal c) => Propagate c d Tensor y x where
    {-# INLINE propagate #-}
    propagate dps qs pq = (dps >$< qs, pq >$> qs)

--instance (Bilinear Tensor y x, Primal c) => Propagate c d Tensor y x where
--    {-# INLINE propagate #-}
--    propagate dps qs pq =
--        let foldfun (dp,q) (k,dpq) = (k+1,(dp >.< q) + dpq)
--         in (uncurry (/>) . foldr foldfun (0,0) $ zip dps qs, pq >$> qs)

instance (Map c d (Affine f) y x, Propagate c d f y x) => Propagate c d (Affine f) y x where
    {-# INLINE propagate #-}
    propagate dps qs pq =
        let (p,pq') = splitAffine pq
            (dpq',ps') = propagate dps qs pq'
         in (joinAffine (average dps) dpq', (p +) <$> ps')


-- Direct Sums --

instance (Legendre x, Legendre y, PotentialCoordinates x ~ PotentialCoordinates y)
  => Legendre (x,y) where
      type PotentialCoordinates (x,y) = PotentialCoordinates x
      {-# INLINE potential #-}
      potential pmn =
          let (pm,pn) = splitPair pmn
           in potential pm + potential pn

--instance Primal c => Legendre c (Sum '[]) where
--    {-# INLINE potential #-}
--    potential _ = 0
--    --potentialDifferential _ = zero

--instance (Legendre c x, Legendre c (Sum xs)) => Legendre c (Sum (x : xs)) where
--    {-# INLINE potential #-}
--    potential pms =
--        let (pm,pms') = splitSum pms
--         in potential pm + potential pms'

instance (Legendre x, KnownNat k) => Legendre (Replicated k x) where
    type PotentialCoordinates (Replicated k x) = PotentialCoordinates x
    {-# INLINE potential #-}
    potential ps =
        S.sum $ mapReplicated potential ps

instance (DuallyFlat x, KnownNat k) => DuallyFlat (Replicated k x) where
    {-# INLINE dualPotential #-}
    dualPotential ps =
        S.sum $ mapReplicated dualPotential ps

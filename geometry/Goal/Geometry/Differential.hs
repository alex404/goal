{-# LANGUAGE UndecidableInstances #-}

-- | Tools for modelling the differential and Riemannian geometry of a
-- 'Manifold'.
module Goal.Geometry.Differential
    ( -- * Riemannian Manifolds
      Riemannian (metric, flat, sharp)
    , euclideanDistance
    -- * Differentiation
    , differential
    , hessian
    , Propagate (propagate)
    -- * Legendre Manifolds
    , Legendre (PotentialCoordinates,potential)
    , DuallyFlat (dualPotential)
    , canonicalDivergence
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


-- | Computes the differential of a function of the coordinates at a point. This
-- functions returns only the resulting 'CotangentVector', without the
-- corresponding 'Point' where the differential was evaluated.
differential
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> c # x
    -> c #* x
{-# INLINE differential #-}
differential f = Point . G.convert . D.grad f . boxCoordinates

-- | Computes the Hessian of a function at a point. This functions returns
-- only the resulting 'CotangentTensor', without the corresponding 'Point' where
-- the Hessian was evaluated.
hessian
    :: Manifold x
    => (forall a. RealFloat a => B.Vector (Dimension x) a -> a)
    -> c # x
    -> c #> Dual c # Tensor x x -- ^ The Differential
{-# INLINE hessian #-}
hessian f p =
    fromMatrix . S.fromRows . G.convert $ G.convert <$> D.hessian f (boxCoordinates p)

-- | A class of functions which can 'propagate' errors. That is, given an error
-- derivative on the output, the input which caused the output, and a
-- function to derive, computes the derivative of the error between the function
-- and target output.
class Map c d f y x => Propagate c d f y x where
    propagate :: [d #* y] -- ^ The error derivatives in 'Dual' coordinates
              -> [c # x] -- ^ A vector of inputs
              -> Function c d # f y x -- ^ The function to differentiate
              -> (Function c d #* f y x, [d # y]) -- ^ The derivative, and function output

-- | Distance between two 'Point's based on the 'Euclidean' metric (l2 norm).
euclideanDistance
    :: Manifold x => c # x -> c # x -> Double
euclideanDistance (Point xs) (Point ys) = S.l2Norm xs ys


--- Riemannian Manifolds ---


-- | 'Riemannian' 'Manifold's are differentiable 'Manifold's where associated
-- with each 'Point' in the 'Manifold' is a 'TangentSpace' with a smoothly
-- varying 'CotangentTensor' known as the 'metric'. 'flat' and 'sharp' correspond to applying this 'metric' to elements of the 'TangentBundle' and 'CotangentBundle', respectively.
class Manifold x => Riemannian c x where
    metric :: c # x -> c #> Dual c # Tensor x x
    flat :: c # x -> c # x -> c #* x
    {-# INLINE flat #-}
    flat p v = metric p >.> v
    sharp :: c # x -> c #* x -> c # x
    {-# INLINE sharp #-}
    sharp p v = inverse (metric p) >.> v


--- Dually Flat Manifolds ---


-- | Although convex analysis is usually developed seperately from differential
-- geometry, it arrises naturally out of the theory of dually flat 'Manifold's.
--
-- A 'Manifold' is 'Legendre' for a particular coordinated system if it is
-- associated with a particular convex function on points of the manifold known
-- as a 'potential'.
class ( Primal (PotentialCoordinates x), Manifold x ) => Legendre x where
    type PotentialCoordinates x :: *
    potential :: PotentialCoordinates x # x -> Double

class Legendre x => DuallyFlat x where
    dualPotential :: PotentialCoordinates x #* x -> Double

-- | Computes the canonical 'divergence' between two points.
canonicalDivergence
    :: DuallyFlat x => PotentialCoordinates x # x -> PotentialCoordinates x #* x -> Double
{-# INLINE canonicalDivergence #-}
canonicalDivergence pp dq = potential pp + dualPotential dq - (pp <.> dq)


--- Instances ---


-- Euclidean --

instance KnownNat k => Riemannian Cartesian (Euclidean k) where
    metric _ = fromMatrix S.matrixIdentity
    flat _ = breakPoint
    sharp _ = breakPoint

-- Replicated Riemannian Manifolds --

--instance {-# OVERLAPPABLE #-} (Riemannian c x, KnownNat k) => Riemannian c (Replicated k x) where
--    metric = error "Do not call metric on a replicated manifold"
--    {-# INLINE flat #-}
--    flat = S.map flat
--    {-# INLINE sharp #-}
--    sharp = S.map sharp

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


-- Direct Sums --

instance (Legendre x, Legendre y, PotentialCoordinates x ~ PotentialCoordinates y)
  => Legendre (x,y) where
      type PotentialCoordinates (x,y) = PotentialCoordinates x
      {-# INLINE potential #-}
      potential pmn =
          let (pm,pn) = splitPair pmn
           in potential pm + potential pn

--    potentialDifferential pmn =
--        let (pm,pn) = splitPair pmn
--         in joinPair (potentialDifferential pm) (potentialDifferential pn)


--instance Primal c => Legendre c (Sum '[]) where
--    {-# INLINE potential #-}
--    potential _ = 0
--    --potentialDifferential _ = zero

--instance (Legendre c x, Legendre c (Sum xs)) => Legendre c (Sum (x : xs)) where
--    {-# INLINE potential #-}
--    potential pms =
--        let (pm,pms') = splitSum pms
--         in potential pm + potential pms'
--    potentialDifferential pms =
--        let (pm,pms') = splitSum pms
--         in joinSum (potentialDifferential pm) (potentialDifferential pms')

instance (Legendre x, KnownNat k) => Legendre (Replicated k x) where
    type PotentialCoordinates (Replicated k x) = PotentialCoordinates x
    {-# INLINE potential #-}
    potential ps =
        S.sum . S.map potential $ splitReplicated ps
--    potentialDifferential ps =
--        breakPoint $ mapReplicatedPoint potentialDifferential ps

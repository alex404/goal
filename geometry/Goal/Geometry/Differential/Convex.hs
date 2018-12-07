{-# LANGUAGE
    DataKinds,
    MultiParamTypeClasses,
    FlexibleInstances,
    FlexibleContexts,
    TypeOperators
    #-}

-- | Here we provide tools for convex analysis based on differential geometry.
-- The dual structures of convex analysis are equivalent to Riemannian manifolds
-- with certain properties.
module Goal.Geometry.Differential.Convex (
    -- * Legendre Manifolds
      Legendre (potential,potentialDifferential)
    , divergence
      -- ** Util
    , dualTransition
    , primalIsomorphism
    , dualIsomorphism
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Storable as S


--- Dually Flat Manifolds ---


-- | Although convex analysis is usually developed seperately from differential
-- geometry, it arrises naturally out of the theory of dually flat 'Manifold's.
--
-- A 'Manifold' is 'Legendre' for a particular coordinated system if it is
-- associated with a particular convex function on points of the manifold known
-- as a 'potential'.
class (Primal c, Manifold m) => Legendre c m where
    potential :: Point c m -> Double
    potentialDifferential :: Point c m -> CotangentVector c m

-- | Transitions a point to its 'Dual' coordinate system.
dualTransition :: Legendre c m => Point c m -> Point (Dual c) m
{-# INLINE dualTransition #-}
dualTransition p =  breakPoint $ potentialDifferential p

-- | Computes the canonical 'divergence' between two points.
divergence
    :: (Legendre c m, Legendre (Dual c) m)
    => Point c m -> Point (Dual c) m -> Double
{-# INLINE divergence #-}
divergence pp dq = potential pp + potential dq - (pp <.> dq)

-- | The 'metric' for a 'Legendre' 'Manifold'. This function can be used to
-- instatiate 'Riemannian' for a 'Legendre' 'Manifold' in a particular
-- coordinate system.
--legendreMetric :: Legendre c m => Point c m -> Point (c ~> Dual c) (Product m m)
--legendreMetric p =  breakPoint $ potentialHessian p

-- | The 'Dual' space of a 'Convex' 'Manifold' is isomorphic to its cotangent
-- space, and we often wish to treat the former as the latter.
primalIsomorphism :: Point c m -> CotangentVector (Dual c) m
{-# INLINE primalIsomorphism #-}
primalIsomorphism (Point xs) = Point xs

-- | The inverse of 'primalIsomorphism'.
dualIsomorphism :: CotangentVector c m -> Point (Dual c) m
{-# INLINE dualIsomorphism #-}
dualIsomorphism (Point xs) =  Point xs


-- Generic --


-- Direct Sums --

instance (Legendre c m, Legendre c n) => Legendre c (m,n) where
    {-# INLINE potential #-}
    potential pmn =
        let (pm,pn) = splitPair pmn
         in potential pm + potential pn
    potentialDifferential pmn =
        let (pm,pn) = splitPair pmn
         in primalIsomorphism $ joinPair (dualIsomorphism (potentialDifferential pm)) (dualIsomorphism (potentialDifferential pn))


instance Primal c => Legendre c (Sum '[]) where
    {-# INLINE potential #-}
    potential _ = 0
    potentialDifferential _ = zero

instance (Legendre c m, Legendre c (Sum ms)) => Legendre c (Sum (m : ms)) where
    {-# INLINE potential #-}
    potential pms =
        let (pm,pms') = splitSum pms
         in potential pm + potential pms'
    potentialDifferential pms =
        let (pm,pms') = splitSum pms
         in primalIsomorphism $ joinSum (dualIsomorphism (potentialDifferential pm)) (dualIsomorphism (potentialDifferential pms'))

instance {-# OVERLAPPABLE #-} (Legendre c m, KnownNat k) => Legendre c (Replicated k m) where
    {-# INLINE potential #-}
    potential ps =
        S.sum . S.map potential $ splitReplicated ps
    potentialDifferential ps =
        breakPoint $ mapReplicatedPoint potentialDifferential ps

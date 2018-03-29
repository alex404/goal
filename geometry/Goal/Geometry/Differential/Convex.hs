{-# LANGUAGE UndecidableSuperClasses,UndecidableInstances #-}

-- | Here we provide tools for convex analysis based on differential geometry.
-- The dual structures of convex analysis are equivalent to Riemannian manifolds
-- with certain properties.
module Goal.Geometry.Differential.Convex (
    -- * Legendre Manifolds
      Legendre (potential)
    , divergence
      -- ** Util
    , dualTransition
    , legendreMetric
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Map
import Goal.Geometry.Map.Multilinear
import Goal.Geometry.Differential


--- Dually Flat Manifolds ---


-- | Although convex analysis is usually developed seperately from differential
-- geometry, it arrises naturally out of the theory of dually flat 'Manifold's.
--
-- A 'Manifold' is 'Legendre' for a particular coordinated system if it is
-- associated with a particular convex function on points of the manifold known
-- as a 'potential'.
class (Primal c, Manifold m) => Legendre c m where
    potential :: RealFloat x => Point c m x -> x

-- | Transitions a point to its 'Dual' coordinate system.
dualTransition :: (Legendre c m, RealFloat x) => Point c m x -> Point (Dual c) m x
{-# INLINE dualTransition #-}
dualTransition p =  Point . coordinates $ differential potential p

-- | Computes the canonical 'divergence' between two points.
divergence
    :: (Legendre c m, Legendre (Dual c) m, RealFloat x)
    => Point c m x -> Point (Dual c) m x -> x
{-# INLINE divergence #-}
divergence pp dq = potential pp + potential dq - (pp <.> dq)

-- | The 'metric' for a 'Legendre' 'Manifold'. This function can be used to
-- instatiate 'Riemannian' for a 'Legendre' 'Manifold' in a particular
-- coordinate system.
legendreMetric :: (Legendre c m, RealFloat x) => Point c m x -> Point (c ~> Dual c) (Product m m) x
legendreMetric p =  Point . coordinates $ hessian potential p


-- Generic --


-- Direct Sums --

instance (Legendre c m, Legendre c n) => Legendre c (Sum m n) where
    {-# INLINE potential #-}
    potential pmn =
        let (pm,pn) = splitSum pmn
         in potential pm + potential pn

instance (Legendre c m, KnownNat k) => Legendre c (Replicated k m) where
    {-# INLINE potential #-}
    potential ps =
        sum $ potential <$> splitReplicated ps

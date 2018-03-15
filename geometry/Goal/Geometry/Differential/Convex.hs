{-# LANGUAGE UndecidableSuperClasses,UndecidableInstances #-}

-- | Here we provide tools for convex analysis based on differential geometry.
-- The dual structures of convex analysis are equivalent to Riemannian manifolds
-- with certain properties.
module Goal.Geometry.Differential.Convex (
    -- * Legendre Manifolds
      Legendre (potential0)
    , potential
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
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G


--- Dually Flat Manifolds ---


-- | Although convex analysis is usually developed seperately from differential
-- geometry, it arrises naturally out of the theory of dually flat 'Manifold's.
--
-- A 'Manifold' is 'Legendre' for a particular coordinated system if it is
-- associated with a particular convex function on points of the manifold known
-- as a 'potential'.
class (Primal c, Manifold m) => Legendre c m where
    potential0 :: RealFloat x => Point c m -> B.Vector (Dimension m) x -> x

potential :: Legendre c m => Point c m -> Double
potential = unboxFunction potential0

-- | Transitions a point to its 'Dual' coordinate system.
dualTransition :: Legendre c m => Point c m -> Point (Dual c) m
{-# INLINE dualTransition #-}
dualTransition p =  Point . coordinates $ differential (potential0 p) p

-- | Computes the canonical 'divergence' between two points.
divergence
    :: (Legendre c m, Legendre (Dual c) m)
    => Point c m -> Point (Dual c) m -> Double
{-# INLINE divergence #-}
divergence pp dq = potential pp + potential dq - (pp <.> dq)

-- | The 'metric' for a 'Legendre' 'Manifold'. This function can be used to
-- instatiate 'Riemannian' for a 'Legendre' 'Manifold' in a particular
-- coordinate system.
legendreMetric :: Legendre c m => Point c m -> CotangentTensor c m
legendreMetric p =  Point . coordinates $ hessian (potential0 p) p


-- Generic --

-- Direct Sums --

instance (Legendre c m, Legendre c n) => Legendre c (Sum m n) where
    {-# INLINE potential0 #-}
    potential0 p xs =
        let (pm,pn) = splitSum p
            (xsm,xsn) = B.splitAt xs
         in potential0 pm xsm + potential0 pn xsn

instance (Legendre c m, KnownNat k) => Legendre c (Replicated k m) where
    {-# INLINE potential0 #-}
    potential0 p x =
        let ps = G.convert $ splitReplicated p
            xs = G.breakEvery x
         in G.sum $ G.zipWith potential0 ps xs

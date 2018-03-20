{-# LANGUAGE UndecidableSuperClasses,UndecidableInstances #-}

-- | Here we provide tools for convex analysis based on differential geometry.
-- The dual structures of convex analysis are equivalent to Riemannian manifolds
-- with certain properties.
module Goal.Geometry.Differential.Convex (
    -- * Legendre Manifolds
      Legendre (bPotential)
    , potential
    , divergence
    , bDivergence
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
    bPotential :: RealFloat x => BPoint c m x -> x

potential :: Legendre c m => Point c m -> Double
potential = bPotential . boxPoint

-- | Transitions a point to its 'Dual' coordinate system.
dualTransition :: Legendre c m => Point c m -> Point (Dual c) m
{-# INLINE dualTransition #-}
dualTransition p =  Point . coordinates $ differential bPotential p

-- | Computes the canonical 'divergence' between two points.
divergence
    :: (Legendre c m, Legendre (Dual c) m)
    => Point c m -> Point (Dual c) m -> Double
{-# INLINE divergence #-}
divergence pp dq = potential pp + potential dq - (pp <.> dq)

-- | Computes the canonical 'divergence' between two points.
bDivergence
    :: (Legendre c m, Legendre (Dual c) m, RealFloat x)
    => BPoint c m x -> BPoint (Dual c) m x -> x
{-# INLINE bDivergence #-}
bDivergence pp dq = bPotential pp + bPotential dq - G.sum (G.zipWith (*) (bCoordinates pp) (bCoordinates dq))

-- | The 'metric' for a 'Legendre' 'Manifold'. This function can be used to
-- instatiate 'Riemannian' for a 'Legendre' 'Manifold' in a particular
-- coordinate system.
legendreMetric :: Legendre c m => Point c m -> Point (c ~> Dual c) (Product m m)
legendreMetric p =  Point . coordinates $ hessian bPotential p


-- Generic --


splitBSum :: (Manifold m, Manifold n) => BPoint c (Sum m n) x -> (BPoint c m x, BPoint c n x)
{-# INLINE splitBSum #-}
splitBSum (BPoint xs) =
    let (xms,xns) = G.splitAt xs
     in (BPoint xms, BPoint xns)

splitBReplicated
    :: (KnownNat k, Manifold m)
    => BPoint c (Replicated k m) x
    -> B.Vector k (BPoint c m x)
{-# INLINE splitBReplicated #-}
splitBReplicated = G.map BPoint . G.breakEvery . bCoordinates



-- Direct Sums --

instance (Legendre c m, Legendre c n) => Legendre c (Sum m n) where
    {-# INLINE bPotential #-}
    bPotential pmn =
        let (pm,pn) = splitBSum pmn
         in bPotential pm + bPotential pn

instance (Legendre c m, KnownNat k) => Legendre c (Replicated k m) where
    {-# INLINE bPotential #-}
    bPotential ps =
        sum $ bPotential <$> splitBReplicated ps

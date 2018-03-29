-- | The Map module provides tools for developing manifolds of functions.
-- A map is a manifold where the points of the manifold represent
-- parametric functions between manifolds. The defining feature of maps is
-- that they have a particular domain and codomain, which themselves are
-- manifolds.

module Goal.Geometry.Map (
     -- * Charts
       Function
     , type (~>)
     -- * Maps
     , Map (Domain, Codomain)
     , Apply ((>.>),(>>.>),(>$>),(>>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Linear

import qualified Goal.Core.Vector.Boxed as B

-- Charts on Maps --

-- | 'Function' Charts help track Charts on the 'Domain' and 'Codomain'. The
-- first Chart corresponds to the chart of the 'Domain', and the second to the
-- chart of the 'Codomain'.
data Function c d

-- | Infix version of 'Function'.
type (c ~> d) = Function c d
infixl 6 ~>

-- | A 'Manifold' is a 'Map' if it has a 'Domain' and a 'Codomain'.
class (Manifold f, Manifold (Domain f), Manifold (Codomain f)) => Map f where
    -- | 'Domain' of the map.
    type Domain f :: *
    -- | 'Codomain' of the map.
    type Codomain f :: *

-- | A 'Manifold' satisfies 'Apply' if it is associated with a function which maps from the 'Domain' to the 'Codomain' of the 'Map'.
class Map f => Apply c d f where
    -- | 'Map' application restricted to doubles.
    (>.>) :: RealFloat x => Point (Function c d) f x -> Point c (Domain f) x -> Point d (Codomain f) x
    (>.>) f x = B.head $ f >$> B.singleton x
    -- | Non AD version
    (>>.>) :: Function c d # f -> c # Domain f -> d # Codomain f
    -- (>>.>) f x = B.head $ f >$> B.singleton x
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply list-mapping (>.>).
    (>$>) :: (RealFloat x, KnownNat k) => Point (Function c d) f x -> B.Vector k (Point c (Domain f) x) -> B.Vector k (Point d (Codomain f) x)
    (>$>) f = B.map (f >.>)
    -- | Non AD version
    (>>$>) :: KnownNat k => Function c d # f -> B.Vector k (c # Domain f) -> B.Vector k (d # Codomain f)
    -- (>>$>) f = B.map (f >.>)


infix 8 >.>
infix 8 >>.>
infix 8 >$>
infix 8 >>$>

instance (Primal c, Primal d) => Primal (Function c d) where
    type Dual (Function c d) = Function (Dual c) (Dual d)

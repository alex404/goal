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
     , Apply ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold

import qualified Goal.Core.Vector.Storable as S

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
    (>.>) :: Point (Function c d) f -> Point c (Domain f) -> Point d (Codomain f)
    (>.>) f x = S.head $ f >$> S.singleton x
    -- | Non AD version
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply list-mapping (>.>).
    (>$>) :: KnownNat k => Point (Function c d) f -> S.Vector k (Point c (Domain f)) -> S.Vector k (Point d (Codomain f))
    (>$>) f = S.map (f >.>)
    -- | Non AD version

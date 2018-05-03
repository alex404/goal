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
     , Map ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Linear

import qualified Goal.Core.Vector.Storable as S

-- Charts on Maps --

-- | 'Function' Charts help track Charts on the 'Domain' and 'Codomain'. The
-- first Chart corresponds to the chart of the 'Domain', and the second to the
-- chart of the 'Codomain'.
data Function c d

-- | Infix version of 'Function'.
type (c ~> d) = Function c d
infixl 6 ~>

-- | A 'Manifold' satisfies 'Apply' if it is associated with a function which maps from the 'Domain' to the 'Codomain' of the 'Map'.
class (Manifold m, Manifold n, Manifold (f m n)) => Map c d f m n where
    -- | 'Map' application restricted to doubles.
    (>.>) :: Function c d # f m n -> c # n -> d # m
    (>.>) f x = S.head . splitReplicated $ f >$> joinReplicated (S.singleton x)
    -- | Non AD version
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply list-mapping (>.>).
    (>$>) :: KnownNat k
          => Function c d # f m n
          -> c # Replicated k n
          -> d # Replicated k m
    (>$>) f = mapReplicatedPoint (f >.>)
    -- | Non AD version

instance (Primal c, Primal d) => Primal (Function c d) where
    type Dual (Function c d) = Function (Dual c) (Dual d)

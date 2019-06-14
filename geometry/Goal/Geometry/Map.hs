-- | Definitions for working with manifolds of functions, a.k.a. function spaces.

module Goal.Geometry.Map (
     -- * Charts
       Function
     , type (#>)
     , type (#*>)
     -- * Maps
     , Map ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Geometry.Manifold
import Goal.Geometry.Linear

-- Charts on Maps --

-- | 'Function' Charts help track Charts on the domain and codomain of the map it parameterizes.
data Function c d

-- | 'Function' between dual coordinate systems.
type (c #> f) = Function (Dual c) c # f
infixl 3 #>

-- | 'Function' between dual coordinate systems.
type (c #*> f) = Function (Dual c) c #* f
infixl 3 #*>

-- | A 'Manifold' is a 'Map' if it is a binary type-function of two `Manifold's, and can transforms 'Point's on the first 'Manifold' into 'Point's on the second 'Manifold'.
class (Manifold x, Manifold y, Manifold (f y x)) => Map c d f y x where
    -- | 'Map' application restricted.
    (>.>) :: Function c d # f y x -> c # x -> d # y
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply mapping (>.>).
    (>$>) :: Function c d # f y x
          -> [c # x]
          -> [d # y]

instance (Primal c, Primal d) => Primal (Function c d) where
    type Dual (Function c d) = Function (Dual c) (Dual d)

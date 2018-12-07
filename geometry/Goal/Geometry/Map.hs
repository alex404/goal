{-# LANGUAGE
    TypeOperators,
    MultiParamTypeClasses,
    TypeFamilies,
    ExplicitNamespaces
    #-}
-- | Definitions for working with manifolds of functions, a.k.a. function spaces.

module Goal.Geometry.Map (
     -- * Charts
       Function
     , type (#>)
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

-- | Infix version of 'Function'.
type (c #> d) = Function c d
infixl 6 #>

-- | A 'Manifold' is a 'Map' if it is a binary type-function of two `Manifold's, and can transforms 'Point's on the first 'Manifold' into 'Point's on the second 'Manifold'.
class (Manifold m, Manifold n, Manifold (f m n)) => Map c d f m n where
    -- | 'Map' application restricted.
    (>.>) :: Function c d # f m n -> c # n -> d # m
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply mapping (>.>).
    (>$>) :: Function c d # f m n
          -> [c # n]
          -> [d # m]

instance (Primal c, Primal d) => Primal (Function c d) where
    type Dual (Function c d) = Function (Dual c) (Dual d)

-- | Definitions for working with manifolds of functions, a.k.a. function spaces.

module Goal.Geometry.Map (
     Map ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Geometry.Manifold
import Goal.Geometry.Vector

-- Charts on Maps --

-- | A 'Manifold' is a 'Map' if it is a binary type-function of two `Manifold's, and can transforms 'Point's on the first 'Manifold' into 'Point's on the second 'Manifold'.
class (Manifold x, Manifold y, Manifold (f y x)) => Map c f y x where
    -- | 'Map' application restricted.
    (>.>) :: c # f y x -> c #* x -> c # y
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply mapping (>.>).
    (>$>) :: c # f y x
          -> [c #* x]
          -> [c # y]

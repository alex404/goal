-- | Definitions for working with manifolds of functions, a.k.a. function spaces.

module Goal.Geometry.Map (
     Map ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Geometry.Manifold
import Goal.Geometry.Vector

-- Charts on Maps --

-- | 'Map' is a 'Manifold' of functions that maps 'Point's from the @x@ 'Manifold' onto 'Point's on the @y@ 'Manifold' (i.e. the input and output manifolds are read from right to left).
class (Manifold x, Manifold y, Manifold (f y x)) => Map c f y x where
    -- | 'Map' application restricted.
    (>.>) :: c # f y x -> c #* x -> c # y
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply mapping (>.>).
    (>$>) :: c # f y x
          -> [c #* x]
          -> [c # y]

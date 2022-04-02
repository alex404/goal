-- | Definitions for working with manifolds of functions, a.k.a. function spaces.

module Goal.Geometry.Map (
     Map ((>.>),(>$>))
     ) where


--- Imports ---


-- Goal --

import Goal.Geometry.Manifold
import Goal.Geometry.Vector

-- Charts on Maps --

-- | A 'Manifold' is a 'Map' if it is a binary type-function of two `Manifold's, and can transforms 'Point's on the second 'Manifold' into 'Point's on the first 'Manifold' (i.e. the input and output manifolds are read from right to left).
class (Manifold x, Manifold y, Manifold (f x y)) => Map c f x y where
    -- | 'Map' application restricted.
    (>.>) :: c # f x y -> c #* y -> c # x
    -- | 'Map' vector application. May sometimes have a more efficient implementation
    -- than simply mapping (>.>).
    (>$>) :: c # f x y
          -> [c #* y]
          -> [c # x]

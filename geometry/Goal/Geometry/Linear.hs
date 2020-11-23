{-# LANGUAGE UndecidableSuperClasses #-}
-- | The Linear module provides the tools for treating a locally Euclidean patch
-- of a manifold as a linear space.
module Goal.Geometry.Linear
    ( -- * Vector Spaces
      (.>)
    , (/>)
    , convexCombination
    -- * Dual Spaces
    , Primal (Dual)
    , type (#*)
    , (<.>)
    , dotMap
    ) where

--- Imports ---

-- Package --

import Goal.Core
import Goal.Geometry.Manifold

import qualified Goal.Core.Vector.Storable as S

--- Vector Spaces on Manifolds ---


-- | Scalar multiplication of points on a manifold.
(.>) :: Double -> c # x -> c # x
{-# INLINE (.>) #-}
(.>) a (Point xs) = Point $ S.scale a xs
infix 7 .>

-- | Scalar division of points on a manifold.
(/>) :: Double -> c # x -> c # x
{-# INLINE (/>) #-}
(/>) a (Point xs) = Point $ S.scale (recip a) xs
infix 7 />

-- | Combination of two 'Point's. Takes the first argument of the second
-- argument, and (1-first argument) of the third argument.
convexCombination :: Manifold x => Double -> c # x -> c # x -> c # x
convexCombination x p1 p2 = x .> p1 + (1-x) .> p2


--- Dual Spaces ---


-- | 'Primal' charts have a 'Dual' coordinate system.
class (Dual (Dual c) ~ c, Primal (Dual c)) => Primal c where
    type Dual c :: Type

-- | A 'Point' on a 'Manifold' in the 'Dual' coordinates of c.
type (c #* x) = Point (Dual c) x
infix 3 #*

-- | '<.>' is the inner product between a dual pair of 'Point's.
(<.>) :: c # x -> c #* x -> Double
{-# INLINE (<.>) #-}
(<.>) p q = S.dotProduct (coordinates p) (coordinates q)

infix 7 <.>

-- | 'dotMap' computes the inner product over a list of dual elements.
dotMap :: Manifold x => c # x -> [c #* x] -> [Double]
{-# INLINE dotMap #-}
dotMap p qs = S.dotMap (coordinates p) (coordinates <$> qs)

-- Cartesian Spaces --

instance Primal Cartesian where
    type Dual Cartesian = Cartesian

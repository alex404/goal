-- | The Linear module provides the tools for treating a given manifold as a
-- linear space. Note that this is not always sound, as the operations may take
-- the points outside of the coordinate system under consideration. In the
-- future I may incorporate bounds checking into these functions.
module Goal.Geometry.Linear
    ( -- * Vector Spaces
      (<+>)
    , (.>)
    , (<->)
    , (/>)
    , averagePoint
    , convexCombination
    -- * Dual Spaces
    , Primal (Dual)
    , Linear ((<.>))
    ) where

--- Imports ---

-- Package --

import Goal.Core
import Goal.Geometry.Manifold

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

--- Vector Spaces on Manifolds ---


-- | Vector addition of points on a manifold.
(<+>) :: (Manifold m, GVector v x, Num x) => Point v c m x -> Point v c m x -> Point v c m x
{-# INLINE (<+>) #-}
(<+>) (Point xs) (Point xs') = Point $ xs + xs'
infixr 6 <+>

-- | Vector subtraction of points on a manifold.
(<->) :: (Manifold m, GVector v x, Num x) => Point v c m x -> Point v c m x -> Point v c m x
{-# INLINE (<->) #-}
(<->) (Point xs) (Point xs') = Point $ xs - xs'
infixr 6 <->

-- | Scalar multiplication of points on a manifold.
(.>) :: (GVector v x, Num x) => x -> Point v c m x -> Point v c m x
{-# INLINE (.>) #-}
(.>) a (Point xs) = Point $ G.map (*a) xs
infix 7 .>

-- | Scalar division of points on a manifold.
(/>) :: (GVector v x, Fractional x) => x -> Point v c m x -> Point v c m x
{-# INLINE (/>) #-}
(/>) a (Point xs) = Point $ G.map (/a) xs
infix 7 />

-- | Combination of two 'Point's. Takes the first argument of the second
-- argument, and (1-first argument) of the third argument.
convexCombination :: (Manifold m, Fractional x, GVector v x) => x -> Point v c m x -> Point v c m x -> Point v c m x
convexCombination x p1 p2 = x .> p1 <+> (1-x) .> p2

-- | Average 'Point' given a collection of 'Point's.
averagePoint :: (Manifold m, Foldable f, Fractional x, GVector v x) => f (Point v c m x) -> Point v c m x
averagePoint = uncurry (/>) . foldr (\p (s,p') -> (s+1,p <+> p')) (0,zero)



--- Dual Spaces ---


-- | 'Primal' charts have a 'Dual' coordinate system.
class (Dual (Dual c)) ~ c => Primal c where
    type Dual c :: *

class Primal c => Linear v c m x where
    (<.>) :: Point v c m x -> Point v (Dual c) m x -> x

infix 7 <.>

-- Cartesian Spaces --

instance Primal Cartesian where
    type Dual Cartesian = Cartesian

-- | '<.>' is the inner product between a dual pair of 'Point's.
instance (Primal c, Numeric x) => Linear SVector c m x where
    {-# INLINE (<.>) #-}
    (<.>) p q = S.dotProduct (coordinates p) (coordinates q)

instance (Primal c, Num x) => Linear BVector c m x where
    {-# INLINE (<.>) #-}
    (<.>) p q = B.dotProduct (coordinates p) (coordinates q)




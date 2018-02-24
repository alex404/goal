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
    , Primal
    , Dual
    , (<.>)
    ) where

--- Imports ---

-- Package --

import Goal.Core
import Goal.Geometry.Manifold


--- Vector Spaces on Manifolds ---


-- | Vector addition of points on a manifold.
(<+>) :: Num x => Point c m x -> Point c m x -> Point c m x
{-# INLINE (<+>) #-}
(<+>) (Point xs) (Point xs') = Point (zipWithV (+) xs xs')
infixr 6 <+>

-- | Vector subtraction of points on a manifold.
(<->) :: Num x => Point c m x -> Point c m x -> Point c m x
{-# INLINE (<->) #-}
(<->) (Point xs) (Point xs') = Point (zipWithV (-) xs xs')
infixr 6 <->

-- | Scalar multiplication of points on a manifold.
(.>) :: Num x => x -> Point c m x -> Point c m x
{-# INLINE (.>) #-}
(.>) a p = (*a) <$> p
infix 7 .>

-- | Scalar division of points on a manifold.
(/>) :: Fractional x => x -> Point c m x -> Point c m x
{-# INLINE (/>) #-}
(/>) a p = (/a) <$> p
infix 7 />

-- | Combination of two 'Point's. Takes the first argument of the second
-- argument, and (1-first argument) of the third argument.
convexCombination :: (Manifold m, Fractional x) => x -> Point c m x -> Point c m x -> Point c m x
convexCombination x p1 p2 = x .> p1 <+> (1-x) .> p2

-- | Average 'Point' given a collection of 'Point's.
averagePoint :: (Manifold m, Foldable f, Fractional x) => f (Point c m x) -> Point c m x
averagePoint = uncurry (/>) . foldr (\p (s,p') -> (s+1,p <+> p')) (0,zero)



--- Dual Spaces ---


-- | 'Primal' charts have a 'Dual' coordinate system.
class (Dual (Dual c)) ~ c => Primal c where
    type Dual c :: *

-- | '<.>' is the inner product between a dual pair of 'Point's.
(<.>) :: Num x => Point c m x -> Point (Dual c) m x -> x
{-# INLINE (<.>) #-}
(<.>) p q = dotProduct (coordinates p) (coordinates q)

infix 7 <.>


-- Cartesian Spaces --

instance Primal Cartesian where
    type Dual Cartesian = Cartesian

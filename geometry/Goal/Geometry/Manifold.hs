{-# LANGUAGE UndecidableInstances,GeneralizedNewtypeDeriving,DeriveTraversable,GADTs,DataKinds #-}
-- | This module provides the core mathematical definitions used by the rest of Goal. The central
-- object is a 'Point' on a 'Manifold'. A 'Manifold' is an object with a 'Dimension', and a 'Point'
-- represents an element of the 'Manifold' in a particular coordinate system, represented by a
-- chart.
module Goal.Geometry.Manifold
    ( -- * Manifolds
    Manifold (Dimension)
    , dimension
    , Dense
    -- ** Combinators
    , Sum
    , Replicated
    , R
    -- * Points
    , Point (Point,coordinates)
    , type (#)
    , listCoordinates
    , breakChart
    -- ** Reshaping Points
    , splitSum
    , joinSum
    , splitReplicated
    , joinReplicated
    , mapReplicated
    -- * Euclidean Manifolds
    , Continuum
    , Euclidean
    -- ** Charts
    , Cartesian
    , Polar
    -- ** Transition
    , Transition (transition)
    -- ** Constructors
    , zero
    ) where

--- Imports ---


-- Goal --

import Goal.Core

--- Manifolds ---

type Dense x = RealFloat x


-- | A geometric object with a certain 'Dimension'.
class KnownNat (Dimension m) => Manifold m where
    type Dimension m :: Nat

dimension0 :: Manifold m => Proxy (Dimension m) -> Proxy m -> Int
dimension0 prxy _ = natValInt prxy

-- | The 'Dimension' of the given 'Manifold'.
dimension :: Manifold m => Proxy m -> Int
dimension = dimension0 Proxy



--- Points ---

-- | A 'Point' on a 'Manifold'. The phantom type @m@ represents the 'Manifold', and the phantom type
-- @c@ represents the coordinate system, or chart, in which the 'Point' is represented. The variable
-- @x@ indicates the type of the coordinates, and is used to support automatic differentation.
newtype Point c m x =
    Point { coordinates :: Vector (Dimension m) x }
    deriving (Eq,Show,Functor,Foldable,Traversable,NFData)

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # m) = Point c m Double
infix 1 #

-- | Returns the 'Coordinates' of the point in list form.
listCoordinates :: Point c m x -> [x]
listCoordinates = toList . coordinates

-- | Throws away the type-level information about the chart of the given 'Point'.
breakChart :: Point c m x -> Point d m x
breakChart (Point xs) = Point xs


-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @n@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @n@.
data Sum m n

-- | An infix representation of 'Sum' 'Manifold's.
--type (m + n) = Sum m n
--infixr 5 +

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of constituent 'Point's.
splitSum :: (Manifold m, Manifold n) => Point c (Sum m n) x -> (Point c m x, Point c n x)
{-# INLINE splitSum #-}
splitSum (Point xs) =
    let (xms,xns) = splitV xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (Manifold m, Manifold n) => Point c m x -> Point c n x -> Point c (Sum m n) x
{-# INLINE joinSum #-}
joinSum (Point xms) (Point xns) =
    Point $ joinV xms xns

-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k m = Replicated k m

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitReplicated :: (KnownNat k, Manifold m) => Point c (Replicated k m) x -> Vector k (Point c m x)
{-# INLINE splitReplicated #-}
splitReplicated (Point xs) = fmap Point . toRows $ Matrix xs

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinReplicated :: (KnownNat k, Manifold m) => Vector k (Point c m x) -> Point c (Replicated k m) x
{-# INLINE joinReplicated #-}
joinReplicated ps = Point . flattenV $ coordinates <$> ps

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicated :: (KnownNat k, Manifold m) => (Point c m x -> a) -> Point c (Replicated k m) x -> Vector k a
{-# INLINE mapReplicated #-}
mapReplicated f rp = f <$> splitReplicated rp

-- Charts on Euclidean Space --

-- | One dimensional 'Euclidean' space.
data Continuum

-- | @n@-dimensional Euclidean space.
data Euclidean (n :: Nat)

-- | 'Cartesian' coordinates on 'Euclidean' space.
data Cartesian

-- | 'Polar' coordinates on 'Euclidean' space.
data Polar

-- | A 'transition' involves taking a point represented by the chart 'c',
-- and re-representing in terms of the chart 'd'. This will usually require
-- recomputation of the coordinates.
class Transition c d m where
    transition :: Dense x => Point c m x -> Point d m x

-- | Creates a point on the given manifold with coordinates given by the zero vector.
zero :: (Manifold m, Num x) => Point c m x
zero = Point $ replicateV 0


--- Instances ---


-- Transition --

instance Transition c c m where
    {-# INLINE transition #-}
    transition = id

-- Combinators --

instance (Manifold m, Manifold n) => Manifold (Sum m n) where
    type Dimension (Sum m n) = Dimension m + Dimension n

instance (KnownNat k, Manifold m) => Manifold (Replicated k m) where
    type Dimension (Replicated k m) = k * Dimension m

-- Euclidean Space --

instance (KnownNat k) => Manifold (Euclidean k) where
    type Dimension (Euclidean k) = k

instance Transition Polar Cartesian (Euclidean 2) where
    {-# INLINE transition #-}
    transition (Point rphi) =
        let [r,phi] = toList rphi
            x = r * cos phi
            y = r * sin phi
         in Point $ doubleton x y

instance Transition Cartesian Polar (Euclidean 2) where
    {-# INLINE transition #-}
    transition (Point xs) =
        let [x,y] = toList xs
            r = sqrt $ (x*x) + (y*y)
            phi = atan2 y x
         in Point $ doubleton r phi

{-# LANGUAGE StandaloneDeriving,UndecidableInstances,GeneralizedNewtypeDeriving,DeriveTraversable #-}
-- | This module provides the core mathematical definitions used by the rest of Goal. The central
-- object is a 'Point' on a 'Manifold'. A 'Manifold' is an object with a 'Dimension', and a 'Point'
-- represents an element of the 'Manifold' in a particular coordinate system, represented by a
-- chart.
module Goal.Geometry.Manifold
    ( -- * Manifolds
    Manifold (Dimension)
    , dimension
    -- ** Combinators
    , Sum
    , Replicated
    , R
    -- * Points
    , Point (Point,coordinates)
    , GPoint
    , BPoint
    , SPoint
    , type (#)
    , listCoordinates
    , convertPoint
    , toBPoint
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
import qualified Goal.Core.Vector.Generic as G

--- Manifolds ---


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
newtype Point c m v x =
    Point { coordinates :: Vector v (Dimension m) x }
    deriving (Eq,Show,Functor,Foldable,Traversable,NFData)

deriving instance (KnownNat (Dimension m), Storable x) => Storable (SPoint c m x)

type GPoint c m v x = (GVector v x, GVector v (Point c m v x), GVector v (Vector v (Dimension m) x))

type BPoint c m x = Point c m BVector x

type SPoint c m x = Point c m SVector x

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # m) = SPoint c m Double
infix 1 #

-- | Returns the 'Coordinates' of the point in list form.
listCoordinates :: GVector v x => Point c m v x -> [x]
listCoordinates = G.toList . coordinates

convertPoint :: (GVector v x, GVector w x) => Point c m v x -> Point c m w x
convertPoint (Point xs) = Point $ G.convert xs

toBPoint :: GVector v x => Point c m v x -> BPoint c m x
toBPoint (Point xs) = Point $ G.convert xs

-- | Throws away the type-level information about the chart of the given 'Point'.
breakChart :: Point c m v x -> Point d m v x
breakChart (Point xs) = Point xs

-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @n@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @n@.
data Sum m n

-- | An infix representation of 'Sum' 'Manifold's.
--type (m + n) = Sum m n
--infixr 5 +

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of constituent 'Point's.
splitSum :: (GVector v x, Manifold m, Manifold n) => Point c (Sum m n) v x -> (Point c m v x, Point c n v x)
{-# INLINE splitSum #-}
splitSum (Point xs) =
    let (xms,xns) = G.splitAt xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (GVector v x, Manifold m, Manifold n) => Point c m v x -> Point c n v x -> Point c (Sum m n) v x
{-# INLINE joinSum #-}
joinSum (Point xms) (Point xns) =
    Point $ xms G.++ xns

-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k m = Replicated k m

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitReplicated
    :: (GPoint c m v x, KnownNat k, Manifold m)
    => Point c (Replicated k m) v x
    -> Vector v k (Point c m v x)
{-# INLINE splitReplicated #-}
splitReplicated = G.map Point . G.breakEvery . coordinates

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinReplicated
    :: (GPoint c m v x, KnownNat k, Manifold m)
    => Vector v k (Point c m v x)
    -> Point c (Replicated k m) v x
{-# INLINE joinReplicated #-}
joinReplicated ps = Point . G.concat $ coordinates `G.map` ps

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicated
    :: (GPoint c m v x, GVector v a, KnownNat k, Manifold m)
    => (Point c m v x -> a) -> Point c (Replicated k m) v x -> Vector v k a
{-# INLINE mapReplicated #-}
mapReplicated f rp = f `G.map` splitReplicated rp

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
class GVector v x => Transition c d m v x where
    transition :: Point c m v x -> Point d m v x

-- | Creates a point on the given manifold with coordinates given by the zero vector.
zero :: (GVector v x, Manifold m, Num x) => Point c m v x
zero = Point $ G.replicate 0


--- Instances ---


-- Transition --

instance GVector v x => Transition c c m v x where
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

instance (GVector v x, Floating x) => Transition Polar Cartesian (Euclidean 2) v x where
    {-# INLINE transition #-}
    transition (Point rphi) =
        let [r,phi] = G.toList rphi
            x = r * cos phi
            y = r * sin phi
         in Point $ G.doubleton x y

instance (GVector v x, RealFloat x) => Transition Cartesian Polar (Euclidean 2) v x where
    {-# INLINE transition #-}
    transition (Point xs) =
        let [x,y] = G.toList xs
            r = sqrt $ (x*x) + (y*y)
            phi = atan2 y x
         in Point $ G.doubleton r phi

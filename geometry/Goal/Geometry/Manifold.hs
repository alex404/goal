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
    -- Boxed Points
    , BPoint (BPoint, bCoordinates)
    , toBPoint
    , fromBPoint
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import qualified Goal.Core.Vector.Storable as S

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
newtype Point c m x =
    Point { coordinates :: S.Vector (Dimension m) x }

newtype BPoint c m x =
    BPoint { bCoordinates :: BVector (Dimension m) x }
    deriving (Eq,Show,Functor,Foldable,Traversable,NFData)

deriving instance (Manifold m, Storable x) => Storable (Point c m x)

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # m) = Point c m Double
infix 1 #

-- | Returns the 'Coordinates' of the point in list form.
listCoordinates :: Storable x => Point c m x -> [x]
listCoordinates = S.toList . coordinates

-- | Throws away the type-level information about the chart of the given 'Point'.
breakChart :: Point c m x -> Point d m x
breakChart (Point xs) = Point xs

toBPoint :: Storable x => Point c m x -> BPoint c m x
toBPoint (Point xs) = BPoint $ convert xs

fromBPoint :: Storable x => BPoint c m x -> Point c m x
fromBPoint (BPoint xs) = Point $ convert xs

-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @n@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @n@.
data Sum m n

-- | An infix representation of 'Sum' 'Manifold's.
--type (m + n) = Sum m n
--infixr 5 +

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of constituent 'Point's.
splitSum :: (Manifold m, Manifold n, Storable x) => Point c (Sum m n) x -> (Point c m x, Point c n x)
{-# INLINE splitSum #-}
splitSum (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (Manifold m, Manifold n, Storable x) => Point c m x -> Point c n x -> Point c (Sum m n) x
{-# INLINE joinSum #-}
joinSum (Point xms) (Point xns) =
    Point $ xms S.++ xns

-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k m = Replicated k m

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitReplicated :: (KnownNat k, Manifold m, Storable x) => Point c (Replicated k m) x -> S.Vector k (Point c m x)
{-# INLINE splitReplicated #-}
splitReplicated (Point xs) = S.map Point . S.toRows $ S.Matrix xs

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinReplicated :: (KnownNat k, Manifold m, Storable x) => S.Vector k (Point c m x) -> Point c (Replicated k m) x
{-# INLINE joinReplicated #-}
joinReplicated ps = Point . S.concat $ coordinates `S.map` ps

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicated
    :: (KnownNat k, Manifold m, Storable x, Storable a)
    => (Point c m x -> a) -> Point c (Replicated k m) x -> S.Vector k a
{-# INLINE mapReplicated #-}
mapReplicated f rp = f `S.map` splitReplicated rp

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
    transition :: (RealFloat x, Numeric x) => Point c m x -> Point d m x

-- | Creates a point on the given manifold with coordinates given by the zero vector.
zero :: (Manifold m, Num x, Storable x) => Point c m x
zero = Point $ S.replicate 0


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
        let [r,phi] = S.toList rphi
            x = r * cos phi
            y = r * sin phi
         in Point $ S.doubleton x y

instance Transition Cartesian Polar (Euclidean 2) where
    {-# INLINE transition #-}
    transition (Point xs) =
        let [x,y] = S.toList xs
            r = sqrt $ (x*x) + (y*y)
            phi = atan2 y x
         in Point $ S.doubleton r phi

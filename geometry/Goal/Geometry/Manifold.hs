{-# LANGUAGE UndecidableInstances,StandaloneDeriving,GeneralizedNewtypeDeriving #-}
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
    , boxCoordinates
    , breakChart
    , fromBoxed
    -- ** Reshaping Points
    , splitSum
    , joinSum
    , splitReplicated
    , splitBoxedReplicated
    , joinReplicated
    , joinBoxedReplicated
    , mapReplicated
    , mapReplicatedPoint
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
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

--- Manifolds ---


-- | A geometric object with a certain 'Dimension'.
class KnownNat (Dimension m) => Manifold m where
    type Dimension m :: Nat

dimension0 :: Manifold m => Proxy (Dimension m) -> Proxy m -> Int
dimension0 prxy _ = natValInt prxy

-- | The 'Dimension' of the given 'Manifold'.
dimension :: Manifold m => Proxy m -> Int
dimension = dimension0 Proxy

-- | Throws away the type-level information about the chart of the given 'Point'.
breakChart :: Point c m -> Point d m
breakChart (Point xs) = Point xs


--- Points ---


-- | A 'Point' on a 'Manifold'. The phantom type @m@ represents the 'Manifold', and the phantom type
-- @c@ represents the coordinate system, or chart, in which the 'Point' is represented. The variable
-- @x@ indicates the type of the coordinates, and is used to support automatic differentation.
newtype Point c m =
    Point { coordinates :: S.Vector (Dimension m) Double }
    deriving (Eq,Show,NFData)

deriving instance (KnownNat (Dimension m)) => Storable (Point c m)

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # m) = Point c m
infix 1 #

-- | Returns the 'Coordinates' of the point in list form.
listCoordinates :: Point c m -> [Double]
listCoordinates = S.toList . coordinates

boxCoordinates :: Point c m -> B.Vector (Dimension m) Double
boxCoordinates =  G.convert . coordinates

fromBoxed :: B.Vector (Dimension m) Double -> Point c m
{-# INLINE fromBoxed #-}
fromBoxed =  Point . G.convert

-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @n@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @n@.
data Sum m n

-- | An infix representation of 'Sum' 'Manifold's.
--type (m + n) = Sum m n
--infixr 5 +

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of constituent 'Point's.
splitSum :: (Manifold m, Manifold n) => Point c (Sum m n) -> (Point c m, Point c n)
{-# INLINE splitSum #-}
splitSum (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (Manifold m, Manifold n) => Point c m -> Point c n -> Point c (Sum m n)
{-# INLINE joinSum #-}
joinSum (Point xms) (Point xns) =
    Point $ xms S.++ xns

-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k m = Replicated k m

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitBoxedReplicated
    :: (KnownNat k, Manifold m)
    => Point c (Replicated k m)
    -> B.Vector k (Point c m)
{-# INLINE splitBoxedReplicated #-}
splitBoxedReplicated = G.convert . S.map Point . S.breakEvery . coordinates

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitReplicated
    :: (KnownNat k, Manifold m)
    => Point c (Replicated k m)
    -> S.Vector k (Point c m)
{-# INLINE splitReplicated #-}
splitReplicated = S.map Point . S.breakEvery . coordinates

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinReplicated
    :: (KnownNat k, Manifold m)
    => S.Vector k (Point c m)
    -> Point c (Replicated k m)
{-# INLINE joinReplicated #-}
joinReplicated ps = Point $ S.concatMap coordinates ps

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinBoxedReplicated
    :: (KnownNat k, Manifold m)
    => B.Vector k (Point c m)
    -> Point c (Replicated k m)
{-# INLINE joinBoxedReplicated #-}
joinBoxedReplicated ps = Point . S.concatMap coordinates $ G.convert ps

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicated
    :: (Storable a, KnownNat k, Manifold m)
    => (Point c m -> a) -> Point c (Replicated k m) -> S.Vector k a
{-# INLINE mapReplicated #-}
mapReplicated f rp = f `S.map` splitReplicated rp

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicatedPoint
    :: (KnownNat k, Manifold m, Manifold n)
    => (Point c m -> Point d n) -> Point c (Replicated k m) -> Point d (Replicated k n)
{-# INLINE mapReplicatedPoint #-}
mapReplicatedPoint f rp = Point . S.concatMap (coordinates . f) $ splitReplicated rp

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
    transition :: Point c m -> Point d m

-- | Creates a point on the given manifold with coordinates given by the zero vector.
zero :: Manifold m => Point c m
{-# INLINE zero #-}
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

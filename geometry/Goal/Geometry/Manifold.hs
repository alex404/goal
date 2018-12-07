{-# LANGUAGE
    UndecidableInstances,
    StandaloneDeriving,
    GeneralizedNewtypeDeriving,
    ExplicitNamespaces,
    TypeOperators,
    KindSignatures,
    DataKinds,
    TypeFamilies,
    FlexibleContexts,
    NoStarIsType,
    MultiParamTypeClasses,
    FlexibleInstances
    #-}
-- | The core mathematical definitions used by the rest of Goal. The central
-- object is a 'Point' on a 'Manifold'. A 'Manifold' is an object with a
-- 'Dimension', and a 'Point' represents an element of the 'Manifold' in a
-- particular coordinate system, represented by a chart.
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
    , breakPoint
    , listCoordinates
    , boxCoordinates
    , fromBoxed
    -- ** Reshaping Points
    , splitSum
    , joinSum
    , splitPair
    , joinPair
    , fromSingletonSum
    , toSingletonSum
    , splitReplicated
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
    , transition2
    -- ** Constructors
    , zero
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

-- Unqualified --

import Foreign.Storable


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
-- @c@ represents the coordinate system, or chart, in which the 'Point' is represented.
newtype Point c m =
    Point { coordinates :: S.Vector (Dimension m) Double }
    deriving (Eq,Show,NFData)

deriving instance (KnownNat (Dimension m)) => Storable (Point c m)

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # m) = Point c m
infix 3 #

-- | Returns the coordinates of the point in list form.
listCoordinates :: Point c m -> [Double]
listCoordinates = S.toList . coordinates

-- | Returns the coordinates of the point as a boxed vector.
boxCoordinates :: Point c m -> B.Vector (Dimension m) Double
boxCoordinates =  G.convert . coordinates

-- | Constructs a point with coordinates given by a boxed vector.
fromBoxed :: B.Vector (Dimension m) Double -> Point c m
{-# INLINE fromBoxed #-}
fromBoxed =  Point . G.convert

-- | Throws away the type-level information about the chart and manifold of the given 'Point'.
breakPoint :: Dimension m ~ Dimension n => Point c m -> Point d n
breakPoint (Point xs) = Point xs


-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @ms@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @ms@.
data Sum (ms :: [Type])

-- | Conversion to a sum manifold.
toSingletonSum :: Manifold m => c # m -> c # Sum '[m]
toSingletonSum = breakPoint

-- | Conversion from a sum manifold.
fromSingletonSum :: Manifold m => c # Sum '[m] -> c # m
fromSingletonSum = breakPoint

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of head and tail 'Point's.
splitSum :: (Manifold m, Manifold (Sum ms)) => c # Sum (m : ms) -> (c # m, c # Sum ms)
{-# INLINE splitSum #-}
splitSum (Point cs) =
    let (cm,cms) = S.splitAt cs
     in (Point cm, Point cms)

-- | Joins a head and tail sum 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (Manifold m, Manifold (Sum ms)) => c # m -> c # Sum ms -> c # Sum (m : ms)
{-# INLINE joinSum #-}
joinSum (Point cm) (Point cms) =
    Point $ cm S.++ cms

-- | Takes a 'Point' on a pair of 'Manifold's and returns the pair of constituent 'Point's.
splitPair :: (Manifold m, Manifold n) => Point c (m,n) -> (Point c m, Point c n)
{-# INLINE splitPair #-}
splitPair (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a pair 'Manifold'.
joinPair :: (Manifold m, Manifold n) => Point c m -> Point c n -> Point c (m,n)
{-# INLINE joinPair #-}
joinPair (Point xms) (Point xns) =
    Point $ xms S.++ xns


-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k m = Replicated k m

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

-- | A combination of 'splitReplicated' and 'fmap', where the value of the mapped function is also a point.
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
-- and re-representing in terms of the chart 'd'.
class Transition c d m where
    transition :: Point c m -> Point d m

-- | Creates a point on the given manifold with coordinates given by the zero vector.
zero :: Manifold m => Point c m
{-# INLINE zero #-}
zero = Point $ S.replicate 0

-- | Generalizes a function of two points in given coordinate systems to a
-- function on arbitrary coordinate systems.
transition2
    :: (Transition c1 d1 m1, Transition c2 d2 m2)
    => (d1 # m1 -> d2 # m2 -> x)
    -> c1 # m1
    -> c2 # m2
    -> x
{-# INLINE transition2 #-}
transition2 f p q =
   f (transition p) (transition q)


--- Instances ---


-- Transition --


-- Combinators --

instance Manifold (Sum '[]) where
    type Dimension (Sum '[]) = 0

instance (Manifold m, Manifold (Sum ms)) => Manifold (Sum (m : ms)) where
    type Dimension (Sum (m : ms)) = Dimension m + Dimension (Sum ms)

instance (Manifold m, Manifold n) => Manifold (m,n) where
    type Dimension (m,n) = Dimension m + Dimension n

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


--- Transitions ---


instance Transition c d (Sum '[]) where
    {-# INLINE transition #-}
    transition _ = zero

instance (Manifold m, Manifold (Sum ms), Transition c d m, Transition c d (Sum ms))
  => Transition c d (Sum (m : ms)) where
    {-# INLINE transition #-}
    transition pms =
        let (pm,pms') = splitSum pms
         in joinSum (transition pm) (transition pms')

instance (Manifold m, Manifold n, Transition c d m, Transition c d n) => Transition c d (m,n) where
    {-# INLINE transition #-}
    transition pmn =
        let (pm,pn) = splitPair pmn
         in joinPair (transition pm) (transition pn)

instance (KnownNat k, Manifold m, Transition c d m) => Transition c d (Replicated k m) where
    {-# INLINE transition #-}
    transition = mapReplicatedPoint transition

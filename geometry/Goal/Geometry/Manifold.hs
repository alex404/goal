{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    UndecidableInstances,
    StandaloneDeriving,
    GeneralizedNewtypeDeriving
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
class KnownNat (Dimension x) => Manifold x where
    type Dimension x :: Nat

dimension0 :: Manifold x => Proxy (Dimension x) -> Proxy x -> Int
dimension0 prxy _ = natValInt prxy

-- | The 'Dimension' of the given 'Manifold'.
dimension :: Manifold x => Proxy x -> Int
dimension = dimension0 Proxy


--- Points ---


-- | A 'Point' on a 'Manifold'. The phantom type @m@ represents the 'Manifold', and the phantom type
-- @c@ represents the coordinate system, or chart, in which the 'Point' is represented.
newtype Point c x =
    Point { coordinates :: S.Vector (Dimension x) Double }
    deriving (Eq,Ord,Show,NFData)

deriving instance (KnownNat (Dimension x)) => Storable (Point c x)
deriving instance (Manifold x, KnownNat (Dimension x)) => Floating (Point c x)
deriving instance (Manifold x, KnownNat (Dimension x)) => Fractional (Point c x)

-- | An infix version of 'Point', where @x@ is assumed to be of type 'Double'.
type (c # x) = Point c x
infix 3 #

-- | Returns the coordinates of the point in list form.
listCoordinates :: Point c x -> [Double]
listCoordinates = S.toList . coordinates

-- | Returns the coordinates of the point as a boxed vector.
boxCoordinates :: Point c x -> B.Vector (Dimension x) Double
boxCoordinates =  G.convert . coordinates

-- | Constructs a point with coordinates given by a boxed vector.
fromBoxed :: B.Vector (Dimension x) Double -> Point c x
{-# INLINE fromBoxed #-}
fromBoxed =  Point . G.convert

-- | Throws away the type-level information about the chart and manifold of the given 'Point'.
breakPoint :: Dimension x ~ Dimension y => Point c x -> Point d y
breakPoint (Point xs) = Point xs


-- Manifold Combinators --


-- | A 'Sum' type for 'Manifold's, such that the 'Sum' of @m@ and @ms@ has 'Dimension' equal to the
-- sum of 'Dimension' @m@ and 'Dimension' @ms@.
data Sum (ms :: [Type])

-- | Conversion to a sum manifold.
toSingletonSum :: Manifold x => c # x -> c # Sum '[x]
toSingletonSum = breakPoint

-- | Conversion from a sum manifold.
fromSingletonSum :: Manifold x => c # Sum '[x] -> c # x
fromSingletonSum = breakPoint

-- | Takes a 'Point' on a 'Sum' 'Manifold' and returns the pair of head and tail 'Point's.
splitSum :: (Manifold x, Manifold (Sum xs)) => c # Sum (x : xs) -> (c # x, c # Sum xs)
{-# INLINE splitSum #-}
splitSum (Point cs) =
    let (cm,cms) = S.splitAt cs
     in (Point cm, Point cms)

-- | Joins a head and tail sum 'Point's into a 'Point' on a 'Sum' 'Manifold'.
joinSum :: (Manifold x, Manifold (Sum xs)) => c # x -> c # Sum xs -> c # Sum (x : xs)
{-# INLINE joinSum #-}
joinSum (Point cm) (Point cms) =
    Point $ cm S.++ cms

-- | Takes a 'Point' on a pair of 'Manifold's and returns the pair of constituent 'Point's.
splitPair :: (Manifold x, Manifold y) => Point c (x,y) -> (Point c x, Point c y)
{-# INLINE splitPair #-}
splitPair (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

-- | Joins a pair of 'Point's into a 'Point' on a pair 'Manifold'.
joinPair :: (Manifold x, Manifold y) => Point c x -> Point c y -> Point c (x,y)
{-# INLINE joinPair #-}
joinPair (Point xms) (Point xns) =
    Point $ xms S.++ xns


-- | A 'Sum' type for repetitions of the same 'Manifold'.
data Replicated (k :: Nat) m

-- | An abbreviation for 'Replicated'.
type R k x = Replicated k x

-- | Splits a 'Point' on a 'Replicated' 'Manifold' into a 'Vector' of of 'Point's.
splitReplicated
    :: (KnownNat k, Manifold x)
    => Point c (Replicated k x)
    -> S.Vector k (Point c x)
{-# INLINE splitReplicated #-}
splitReplicated = S.map Point . S.breakEvery . coordinates

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinReplicated
    :: (KnownNat k, Manifold x)
    => S.Vector k (Point c x)
    -> Point c (Replicated k x)
{-# INLINE joinReplicated #-}
joinReplicated ps = Point $ S.concatMap coordinates ps

-- | Joins a 'Vector' of of 'Point's into a 'Point' on a 'Replicated' 'Manifold'.
joinBoxedReplicated
    :: (KnownNat k, Manifold x)
    => B.Vector k (Point c x)
    -> Point c (Replicated k x)
{-# INLINE joinBoxedReplicated #-}
joinBoxedReplicated ps = Point . S.concatMap coordinates $ G.convert ps

-- | A combination of 'splitReplicated' and 'fmap'.
mapReplicated
    :: (Storable a, KnownNat k, Manifold x)
    => (Point c x -> a) -> Point c (Replicated k x) -> S.Vector k a
{-# INLINE mapReplicated #-}
mapReplicated f rp = f `S.map` splitReplicated rp

-- | A combination of 'splitReplicated' and 'fmap', where the value of the mapped function is also a point.
mapReplicatedPoint
    :: (KnownNat k, Manifold x, Manifold y)
    => (Point c x -> Point d y) -> Point c (Replicated k x) -> Point d (Replicated k y)
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
class Transition c d x where
    transition :: c # x -> d # x

-- | Generalizes a function of two points in given coordinate systems to a
-- function on arbitrary coordinate systems.
transition2
    :: (Transition cx dx x, Transition cy dy y)
    => (dx # x -> dy # y -> a)
    -> cx # x
    -> cy # y
    -> a
{-# INLINE transition2 #-}
transition2 f p q =
   f (transition p) (transition q)


--- Instances ---


-- Transition --


-- Combinators --

instance Manifold (Sum '[]) where
    type Dimension (Sum '[]) = 0

instance (Manifold x, Manifold (Sum xs)) => Manifold (Sum (x : xs)) where
    type Dimension (Sum (x : xs)) = Dimension x + Dimension (Sum xs)

instance (Manifold x, Manifold y) => Manifold (x,y) where
    type Dimension (x,y) = Dimension x + Dimension y

instance (KnownNat k, Manifold x) => Manifold (Replicated k x) where
    type Dimension (Replicated k x) = k * Dimension x

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
    transition _ = 0

instance (Manifold x, Manifold (Sum xs), Transition c d x, Transition c d (Sum xs))
  => Transition c d (Sum (x : xs)) where
    {-# INLINE transition #-}
    transition cxs =
        let (cx,cxs') = splitSum cxs
         in joinSum (transition cx) (transition cxs')

instance (Manifold x, Manifold y, Transition c d x, Transition c d y) => Transition c d (x,y) where
    {-# INLINE transition #-}
    transition cxy =
        let (cx,cy) = splitPair cxy
         in joinPair (transition cx) (transition cy)

instance (KnownNat k, Manifold x, Transition c d x) => Transition c d (Replicated k x) where
    {-# INLINE transition #-}
    transition = mapReplicatedPoint transition


--- Numeric Classes ---


instance (Manifold x, KnownNat (Dimension x)) => Num (c # x) where
    {-# INLINE (+) #-}
    (+) (Point xs) (Point xs') = Point $ S.add xs xs'
    {-# INLINE (*) #-}
    (*) (Point xs) (Point xs') = Point $ xs * xs'
    {-# INLINE negate #-}
    negate (Point xs) = Point $ S.scale (-1) xs
    {-# INLINE abs #-}
    abs (Point xs) = Point $ abs xs
    {-# INLINE signum #-}
    signum (Point xs) = Point $ signum xs
    {-# INLINE fromInteger #-}
    fromInteger x = Point . S.replicate $ fromInteger x


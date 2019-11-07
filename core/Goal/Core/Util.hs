-- | A collection of generic numerical and list manipulation functions.
module Goal.Core.Util
    ( -- * List Manipulation
      takeEvery
    , breakEvery
    , kFold
    -- * Numeric
    , roundSD
    , toPi
    , circularDistance
    , integrate
    , logistic
    , logit
    , square
    , triangularNumber
    -- ** List Numerics
    , average
    , weightedAverage
    , circularAverage
    , weightedCircularAverage
    , range
    , discretizeFunction
    , logSumExp
    , logIntegralExp
    -- * Tracing
    , traceGiven
    -- * TypeNats
    , finiteInt
    , natValInt
    , Triangular
    -- ** Type Rationals
    , Rat
    , type (/)
    , ratVal
    ) where


--- Imports ---


-- Unqualified --

import Numeric
import Data.Ratio
import Data.Proxy
import Debug.Trace
import Data.Finite
import GHC.TypeNats

-- Qualified --

import qualified Numeric.GSL.Integration as I
import qualified Data.List as L

--- General Functions ---


-- | Takes every nth element, starting with the head of the list.
takeEvery :: Int -> [x] -> [x]
{-# INLINE takeEvery #-}
takeEvery m = map snd . filter (\(x,_) -> mod x m == 0) . zip [0..]

-- | Break the list up into lists of length n.
breakEvery :: Int -> [x] -> [[x]]
{-# INLINE breakEvery #-}
breakEvery _ [] = []
breakEvery n xs = take n xs : breakEvery n (drop n xs)

-- | Runs traceShow on the given element.
traceGiven :: Show a => a -> a
traceGiven a = traceShow a a


--- Numeric ---

-- | Numerically integrates a 1-d function over an interval.
integrate
    :: Double -- ^ Error Tolerance
    -> (Double -> Double) -- ^ Function
    -> Double -- ^ Interval beginning
    -> Double -- ^ Interval end
    -> (Double,Double) -- ^ Integral
integrate errbnd = I.integrateQAGS errbnd 10000

-- | Rounds the number to the specified significant digit.
roundSD :: RealFloat x => Int -> x -> x
{-# INLINE roundSD #-}
roundSD n x =
    let n' :: Int
        n' = round $ 10^n * x
     in fromIntegral n'/10^n

-- | Value of a point on a circle, minus rotations.
toPi :: RealFloat x => x -> x
{-# INLINE toPi #-}
toPi x =
    let xpi = x / (2*pi)
        f = xpi - fromIntegral (floor xpi :: Int)
     in 2 * pi * f

-- | Distance between two points on a circle, removing rotations.
circularDistance :: RealFloat x => x -> x -> x
{-# INLINE circularDistance #-}
circularDistance x y =
    let x' = toPi x
        y' = toPi y
     in min (toPi $ x' - y') (toPi $ y' - x')

-- | A standard sigmoid function.
logistic :: Floating x => x -> x
{-# INLINE logistic #-}
logistic x = 1 / (1 + exp (negate x))

-- | The inverse of the logistic.
logit :: Floating x => x -> x
{-# INLINE logit #-}
logit x = log $ x / (1 - x)

-- | The square of a number (for avoiding endless default values).
square :: Floating x => x -> x
{-# INLINE square #-}
square x = x^(2::Int)

-- | Triangular number.
triangularNumber :: Int -> Int
{-# INLINE triangularNumber #-}
triangularNumber n = flip div 2 $ n * (n+1)


-- Lists --

-- | Average value of a 'Traversable' of 'Fractional's.
average :: (Foldable f, Fractional x) => f x -> x
{-# INLINE average #-}
average = uncurry (/) . foldr (\e (s,c) -> (e+s,c+1)) (0,0)

-- | Weighted Average given a 'Traversable' of (weight,value) pairs.
weightedAverage :: (Foldable f, Fractional x) => f (x,x) -> x
{-# INLINE weightedAverage #-}
weightedAverage = uncurry (/) . foldr (\(w,x) (sm,nrm) -> (sm + w*x,nrm + w)) (0,0)

-- | Circular average value of a 'Traversable' of radians.
circularAverage :: (Traversable f, RealFloat x) => f x -> x
{-# INLINE circularAverage #-}
circularAverage rds =
    let snmu = average $ sin <$> rds
        csmu = average $ cos <$> rds
     in atan2 snmu csmu

-- | Returns (validation,training) pairs
kFold :: Int -> [x] -> [([x],[x])]
{-# INLINE kFold #-}
kFold k xs =
    let nvls = ceiling . (/(fromIntegral k :: Double)) . fromIntegral $ length xs
     in L.unfoldr unfoldFun ([], breakEvery nvls xs)
    where unfoldFun (_,[]) = Nothing
          unfoldFun (hds,tl:tls) = Just ((tl,concat $ hds ++ tls),(tl:hds,tls))

-- | Weighted Circular average value of a 'Traversable' of radians.
weightedCircularAverage :: (Traversable f, RealFloat x) => f (x,x) -> x
{-# INLINE weightedCircularAverage #-}
weightedCircularAverage wxs =
    let snmu = weightedAverage $ sinPair <$> wxs
        csmu = weightedAverage $ cosPair <$> wxs
     in atan2 snmu csmu
    where sinPair (w,rd) = (w,sin rd)
          cosPair (w,rd) = (w,cos rd)

-- | Returns n numbers which uniformly partitions the interval [mn,mx].
range
    :: RealFloat x => x -> x -> Int -> [x]
{-# INLINE range #-}
range _ _ 0 = []
range mn mx 1 = [(mn + mx) / 2]
range mn mx n =
    [ x * mx + (1 - x) * mn | x <- (/ (fromIntegral n - 1)) . fromIntegral <$> [0 .. n-1] ]

-- | Takes range information in the form of a minimum, maximum, and sample count,
-- a function to sample, and returns a list of pairs (x,f(x)) over the specified
-- range.
discretizeFunction :: Double -> Double -> Int -> (Double -> Double) -> [(Double,Double)]
{-# INLINE discretizeFunction #-}
discretizeFunction mn mx n f =
    let rng = range mn mx n
    in zip rng $ f <$> rng

-- | Given a set of values, computes the "soft maximum" by way of taking the
-- exponential of every value, summing the results, and then taking the
-- logarithm. Incorporates some tricks to improve numerical stability.
logSumExp :: (Ord x, Floating x, Traversable f) => f x -> x
{-# INLINE logSumExp #-}
logSumExp xs =
    let mx = maximum xs
     in (+ mx) . log1p . subtract 1 . sum $ exp . subtract mx <$> xs

-- | Given a function, computes the "soft maximum" of the function by computing
-- the integral of the exponential of the function, and taking the logarithm of
-- the result. The maximum is first approximated on a given set of samples to
-- improve numerical stability. Pro tip: If you want to compute the normalizer
-- of a an exponential family probability density, provide the unnormalized
-- log-density to this function.
logIntegralExp
    :: Traversable f
    => Double -- ^ Error Tolerance
    -> (Double -> Double) -- ^ Function
    -> Double -- ^ Interval beginning
    -> Double -- ^ Interval end
    -> f Double -- ^ Samples (for approximating the max)
    -> Double -- ^ Log-Integral-Exp
{-# INLINE logIntegralExp #-}
logIntegralExp err f mnbnd mxbnd xsmps =
    let mx = maximum $ f <$> xsmps
        expf x = exp $ f x - mx
     in (+ mx) . log1p . subtract 1 . fst $ integrate err expf mnbnd mxbnd


--- TypeLits ---


-- | Type-level triangular number.
type Triangular n = Div (n * (n + 1)) 2

-- | Type level rational numbers. This implementation does not currently permit negative numbers.
data Rat (n :: Nat) (d :: Nat)

-- | Infix 'Rat'.
type (/) n d = Rat n d

-- | Recover a rational value from a 'Proxy'.
ratVal :: (KnownNat n, KnownNat d) => Proxy (n / d) -> Rational
{-# INLINE ratVal #-}
ratVal = ratVal0 Proxy Proxy


-- | 'natVal and 'fromIntegral'.
natValInt :: KnownNat n => Proxy n -> Int
{-# INLINE natValInt #-}
natValInt = fromIntegral . natVal

-- | 'getFinite' and 'fromIntegral'.
finiteInt :: KnownNat n => Finite n -> Int
{-# INLINE finiteInt #-}
finiteInt = fromIntegral . getFinite

ratVal0 :: (KnownNat n, KnownNat d) => Proxy n -> Proxy d -> Proxy (n / d) -> Rational
ratVal0 prxyn prxyd _ = fromIntegral (natVal prxyn) % fromIntegral (natVal prxyd)

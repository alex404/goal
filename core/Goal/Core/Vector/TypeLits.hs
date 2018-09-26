-- | A few tools and exports for working with type level values.
module Goal.Core.Vector.TypeLits
    ( -- * TypeLits
      natValInt
    , finiteInt
    -- * Type Rationals
    , Rat
    , type (/)
    , ratVal
    -- * Generation by Proxy
    , withNat
    ) where


--- Imports ---


import GHC.TypeLits

import Data.Proxy
import Data.Ratio

import Data.Finite.Internal
import GHC.TypeLits.Singletons

-- | Type level rational numbers. This implementation does not currently permit negative numbers.
data Rat (n :: Nat) (d :: Nat)

-- | Infix 'Rat'.
type (/) n d = Rat n d

-- | Recover a rational value from a 'Proxy'.
ratVal :: (KnownNat n, KnownNat d) => Proxy (n / d) -> Rational
ratVal = ratVal0 Proxy Proxy

-- | Recover an 'Int' (as opposed to 'Integer') from a 'Proxy'.
natValInt :: forall n proxy . KnownNat n => proxy n -> Int
{-# INLINE natValInt #-}
natValInt = fromInteger . natVal

-- | Recover an 'Int' (as opposed to 'Integer') from a 'Finite' value.
finiteInt :: Finite n -> Int
finiteInt (Finite n) = fromInteger n

ratVal0 :: (KnownNat n, KnownNat d) => Proxy n -> Proxy d -> Proxy (n / d) -> Rational
ratVal0 prxyn prxyd _ = natVal prxyn % natVal prxyd

--generateP0 :: forall i x . KnownNat i => Proxy i -> (forall j . KnownNat j =>  Proxy j -> x) -> [x]
--{-# INLINE generateP0 #-}
--generateP0 prxi f = f prxi : generateP0 (Proxy :: Proxy (i+1)) f
--
--generateP :: (forall j . KnownNat j => Proxy j -> x) -> [x]
--{-# INLINE generateP #-}
--generateP = generateP0 (Proxy :: Proxy 0)
--
--generatePM0
--    :: forall i x m . (KnownNat i, Monad m)
--    => Int
--    -> Proxy i
--    -> (forall j . KnownNat j => Proxy j -> m x)
--    -> m [x]
--{-# INLINE generatePM0 #-}
--generatePM0 k prxi f
--    | k == natValInt (Proxy :: Proxy i) = return []
--    | otherwise = do
--        x <- f prxi
--        (x :) <$> generatePM0 k (Proxy :: Proxy (i+1)) f
--
--generatePM
--    :: forall x m . Monad m
--    => Int
--    -> (forall j . KnownNat j => Proxy j -> m x)
--    -> m [x]
--{-# INLINE generatePM #-}
--generatePM k = generatePM0 k (Proxy :: Proxy 0)
--

--- With Integers ---

withNat
    :: Int
    -> (forall j . KnownNat j => Proxy j -> x)
    -> x
withNat k = withNat0 k PeanoZero

withNat0
    :: forall k x . KnownNat k
    => Int
    -> NatPeano k
    -> (forall j . KnownNat j => Proxy j -> x)
    -> x
withNat0 0 _ f = f (Proxy :: Proxy k)
withNat0 k np f = withNat0 (k-1) (PeanoSucc np) f



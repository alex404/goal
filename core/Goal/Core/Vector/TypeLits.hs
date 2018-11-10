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
    , withNat1
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

--- With Integers ---

-- | Note that this will go into an infinite loop if not fed a natural number >= 0.
withNat
    :: Int
    -> (forall j . KnownNat j => Proxy j -> x)
    -> x
{-# INLINE withNat #-}
withNat k = withNat0 k PeanoZero

withNat0
    :: forall k x . KnownNat k
    => Int
    -> NatPeano k
    -> (forall j . KnownNat j => Proxy j -> x)
    -> x
withNat0 0 _ f = f (Proxy :: Proxy k)
withNat0 k np f = withNat0 (k-1) (PeanoSucc np) f

-- | Note that this will go into an infinite loop if not fed a natural number >= 1.
withNat1
    :: Int
    -> (forall j . (KnownNat j, 1 <= j) => Proxy j -> x)
    -> x
{-# INLINE withNat1 #-}
withNat1 k = withNat10 k (PeanoSucc PeanoZero)

withNat10
    :: forall k x . (KnownNat k, 1 <= k)
    => Int
    -> NatPeano k
    -> (forall j . (KnownNat j, 1 <= j) => Proxy j -> x)
    -> x
withNat10 1 _ f = f (Proxy :: Proxy k)
withNat10 k np f = withNat10 (k-1) (PeanoSucc np) f



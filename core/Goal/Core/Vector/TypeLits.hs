{-# LANGUAGE ExplicitNamespaces,RankNTypes,KindSignatures,TypeOperators,DataKinds #-}
-- | A few tools and exports for working with type level values.
module Goal.Core.Vector.TypeLits
    ( -- * TypeLits
      natValInt
    , finiteInt
    -- * Type Rationals
    , Rat
    , type (/)
    , ratVal
    ) where


--- Imports ---


import GHC.TypeLits

import Data.Proxy
import Data.Ratio

import Data.Finite.Internal

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

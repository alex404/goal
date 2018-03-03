-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLVM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector.TypeLits
    ( -- * TypeLits
      natValInt
    -- * Type Rationals
    , Rat
    , type (/)
    , ratVal
    ) where


--- Imports ---


import GHC.TypeLits
import Data.Proxy
import Data.Ratio


-- | Type level rational numbers. This implementation does not currently permit negative numbers.
data Rat (n :: Nat) (d :: Nat)

-- | Infix 'Rat'.
type (/) n d = Rat n d

-- | Recover a rational value from a 'Proxy'.
ratVal :: (KnownNat n, KnownNat d) => Proxy (n / d) -> Rational
ratVal = ratVal0 Proxy Proxy

-- | Recover an 'Int' (as opposed to 'Integer') from a Proxy.
natValInt :: forall n proxy . KnownNat n => proxy n -> Int
{-# INLINE natValInt #-}
natValInt = fromIntegral . natVal

ratVal0 :: (KnownNat n, KnownNat d) => Proxy n -> Proxy d -> Proxy (n / d) -> Rational
ratVal0 prxyn prxyd _ = natVal prxyn % natVal prxyd

-- | Vectors and Matrices with statically typed dimensions. The 'Vector' and 'Matrix' types are
-- newtypes built on 'Data.Vector', so that GHC reduces all incumbent computations to computations
-- on the highly optimized @vector@ library.
--
-- In my provided benchmarks, my implementation of matrix x matrix multiplication performs about 20%
-- faster than the native implementation provided by the @matrix@ library, and performs within a
-- factor of 2-10 of @hmatrix@. This performance can likely be further improved by compiling with
-- the LLVM backend. Moreover, because the provided 'Vector' and 'Matrix' types are 'Traversable',
-- they may support automatic differentiation with the @ad@ library.
module Goal.Core.Vector.Generic.Mutable
    ( -- * Vector
      module Data.Vector.Generic.Mutable.Sized
    ) where


--- Imports ---


import Data.Vector.Generic.Mutable.Sized

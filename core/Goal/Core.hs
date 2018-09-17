{-# LANGUAGE ExplicitNamespaces #-}

-- | This is the base module of @goal-core@, and is the only module that should be imported in order
-- to use the functionality provided in this library.
module Goal.Core
    ( -- * Module Exports
      module Goal.Core.Plot
    , module Goal.Core.Util
    , module Goal.Core.HList
    , module Goal.Core.Vector.TypeLits
    , module Data.Function
    , module Data.Functor
    , module Data.Foldable
    , module Data.Traversable
    , module Data.Ord
    , module Data.Maybe
    , module Data.Either
    , module Data.Default.Class
    , module Data.Finite
    , module Control.Applicative
    , module Control.Monad
    , module Control.Monad.Primitive
    , module Control.Monad.ST
    , module Control.Arrow
    , module Control.Lens.Type
    , module Control.Lens.Getter
    , module Control.Lens.Setter
    , module Control.Lens.TH
    , module Control.Concurrent
    , module Numeric
    , module Numeric.SpecFunctions
    , module Debug.Trace
    , module System.Directory
    , module Control.DeepSeq
    , module GHC.TypeLits
    , module GHC.Generics
    , module Data.Proxy
    , module Data.Csv
    -- * Types and Classes
    , Matrix (Matrix)
    , Numeric
    , Storable
    , Field
    , convert
    ) where


--- Imports ---


-- Re-exports --

import Goal.Core.Plot hiding (empty,over)
import Goal.Core.Util
import Goal.Core.HList
import Goal.Core.Vector.TypeLits

import Data.Functor
import Data.Foldable
import Data.Traversable
import Data.Ord
import Data.Function hiding ((&))
import Data.Maybe
import Data.Either
import Data.Proxy
import Data.Finite

import Control.Applicative hiding (empty)
import Control.Arrow hiding ((<+>))
import Control.Monad
import Control.Monad.ST
import Control.Lens.Type
import Control.Lens.Getter
import Control.Lens.Setter hiding (Identity,argument)
import Control.Lens.TH
import Control.Concurrent
import Control.DeepSeq hiding (force)
import Control.Monad.Primitive hiding (internal)

import GHC.TypeLits
import GHC.Generics (Generic)

import Debug.Trace
import Data.Default.Class
import System.Directory
import Numeric hiding (log1p,expm1)
import Numeric.SpecFunctions
import Data.Csv (FromNamedRecord,ToNamedRecord,DefaultOrdered)


import Foreign.Storable (Storable)
import Numeric.LinearAlgebra (Numeric,Field)
import Goal.Core.Vector.Generic

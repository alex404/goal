{-# LANGUAGE ExplicitNamespaces #-}

-- | This is the base module of @goal-core@, and is the only module that should be imported in order
-- to use the functionality provided in this library.
module Goal.Core
    ( -- * Module Exports
      module Goal.Core.Plot
    , module Goal.Core.Util
    , module Goal.Core.Vector
    , module Goal.Core.Vector.TypeLits
    , module Data.Function
    , module Data.Functor
    , module Data.Foldable
    , module Data.Traversable
    , module Data.Ord
    , module Data.Maybe
    , module Data.Either
    , module Data.Default.Class
    , module Control.Applicative
    , module Control.Monad
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
    , module GHC.TypeLits.Extra
    , module Data.Proxy
    ) where


--- Imports ---


-- Re-exports --

import Goal.Core.Plot hiding (empty,over)
import Goal.Core.Util
import Goal.Core.Vector
import Goal.Core.Vector.TypeLits

import Data.Functor
import Data.Foldable
import Data.Traversable
import Data.Ord
import Data.Function hiding ((&))
import Data.Maybe
import Data.Either
import Data.Proxy

import Control.Applicative hiding (empty)
import Control.Arrow hiding ((<+>))
import Control.Monad
import Control.Monad.ST
import Control.Lens.Type
import Control.Lens.Getter
import Control.Lens.Setter hiding (Identity)
import Control.Lens.TH
import Control.Concurrent
import Control.DeepSeq hiding (force)
import GHC.TypeLits
import GHC.TypeLits.Extra

import Debug.Trace
import Data.Default.Class
import System.Directory
import Numeric hiding (log1p,expm1)
import Numeric.SpecFunctions

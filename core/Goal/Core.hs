{- | This module re-exports a number of (compatible) modules from across base
and other libraries, as well as most of the modules in @goal-core@. It does
not re-export the Vector modules, which should be imported with
qualification.
-}
module Goal.Core (
    -- * Module Exports
    module Goal.Core.Util,
    module Goal.Core.Project,
    module Goal.Core.Circuit,
    module Data.Function,
    module Data.Functor,
    module Data.Foldable,
    module Data.Traversable,
    module Data.Ord,
    module Data.Maybe,
    module Data.Either,
    module Data.Finite,
    module Data.Csv,
    module Data.Proxy,
    module Data.Kind,
    module Data.Functor.Identity,
    module Data.Type.Equality,
    module Control.Applicative,
    module Control.Monad,
    module Control.Monad.Primitive,
    module Control.Monad.ST,
    module Control.Arrow,
    module Control.Concurrent,
    module Control.DeepSeq,
    module Numeric,
    module Numeric.SpecFunctions,
    module Options.Applicative,
    module GHC.TypeNats,
    module GHC.Generics,
    module Debug.Trace,
    module System.Directory,

    -- * (Re)names
    NatNumber,
    ByteString,
    orderedHeader,
) where

--- Imports ---

-- Re-exports --

import Goal.Core.Circuit
import Goal.Core.Project
import Goal.Core.Util

import Data.Csv hiding (Field, Parser, header)
import Data.Csv qualified as CSV
import Data.Either
import Data.Finite
import Data.Foldable
import Data.Function hiding ((&))
import Data.Functor
import Data.Functor.Identity
import Data.Kind (Type)
import Data.Maybe
import Data.Ord
import Data.Proxy
import Data.Traversable
import Data.Type.Equality

import Control.Applicative hiding (empty)
import Control.Arrow hiding ((<+>))
import Control.Concurrent
import Control.DeepSeq hiding (force)
import Control.Monad hiding (join)
import Control.Monad.Primitive hiding (internal)
import Control.Monad.ST

import Options.Applicative

import GHC.Generics (Generic)
import GHC.TypeNats hiding (Mod, Natural)
import GHC.TypeNats as N

import Debug.Trace
import Numeric hiding (expm1, log1p)
import Numeric.SpecFunctions
import System.Directory

import Data.ByteString (ByteString)

type NatNumber = N.Natural

orderedHeader :: [ByteString] -> Header
orderedHeader = CSV.header

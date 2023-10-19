{- | This module re-exports a number of (compatible) modules from across base
and other libraries, as well as most of the modules in @goal-core@. It does
not re-export the Vector modules, which should be imported with
qualification.
-}
module Goal.Core (
    -- * Module Exports
    module Goal.Core.Util,
    module Goal.Core.Circuit,
    module Goal.Core.Project,
    module GHC.TypeNats,
) where

--- Imports ---

import Goal.Core.Circuit
import Goal.Core.Project
import Goal.Core.Util

--- Re-exports ---

import GHC.TypeNats

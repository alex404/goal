{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    RankNTypes,
    TypeOperators,
    FlexibleContexts,
    ScopedTypeVariables
#-}
-- | The main module of goal-probability. Import this module to use all the
-- types, functions, and classes provided by goal-probability.
module Goal.Graphical
    ( -- * Package Exports
      module Goal.Graphical.Models
    , module Goal.Graphical.Models.Dynamic
    , module Goal.Graphical.Models.Harmonium
    , module Goal.Graphical.Models.Harmonium.FactorAnalysis
    , module Goal.Graphical.Learning
    , module Goal.Graphical.Inference
    ) where


--- Imports ---


-- Re-exports --

import Goal.Graphical.Models
import Goal.Graphical.Models.Dynamic
import Goal.Graphical.Models.Harmonium
import Goal.Graphical.Models.Harmonium.FactorAnalysis
import Goal.Graphical.Learning
import Goal.Graphical.Inference

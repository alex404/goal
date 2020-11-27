{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Graphical.Generative
    ( -- * Latent Variable Models
      LatentVariable
    , ExpectationMaximization (expectationStep)
    ) where

--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Latent Variable Class ---


type LatentVariable x o l =
    ( SamplePoint x ~ HList (o : l) )

class Statistical x => ExpectationMaximization c x where
    expectationStep
        :: LatentVariable x o l
        => [o] -> c # x -> Mean # x


{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Graphical.Generative
    ( Observation
    , Observations
    -- * Latent Variable Models
    , ExpectationMaximization (expectationStep)
    -- ** Hierarchical Models
    , ObservablyContinuous
        ( logObservableDensity
        , logObservableDensities
        , observableDensity
        , observableDensities )
    ) where

--- Imports ---


-- Package --

import Goal.Geometry
import Goal.Probability


--- Latent Variable Class ---

type family Observation f

type Observations f = [Observation f]

class ExpectationMaximization f where
    expectationStep :: Observations f -> Natural # f -> Mean # f

class ObservablyContinuous c f where
    logObservableDensity :: c # f -> Observation f -> Double
    logObservableDensity p = head . logObservableDensities p . (:[])
    logObservableDensities :: c # f -> Observations f -> [Double]
    logObservableDensities p = map (logObservableDensity p)

    observableDensity :: c # f -> Observation f -> Double
    observableDensity p = exp . head . logObservableDensities p . (:[])
    observableDensities :: c # f -> Observations f -> [Double]
    observableDensities p = map (exp . logObservableDensity p)



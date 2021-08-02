{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | A few general definitions for Graphical Models.
module Goal.Graphical.Models
    ( Observation
    , Observations
    -- ** Hierarchical Models
    , ObservablyContinuous ( logObservableDensities, observableDensities )
    , logObservableDensity
    , observableDensity
    ) where

--- Imports ---


-- Package --

import Goal.Geometry

--- Latent Variable Class ---

-- | An observation from a latent variable model.
type family Observation f

-- | A list of observations.
type Observations f = [Observation f]

-- | Probability densities over observations in a latent variable model.
class ObservablyContinuous c f where
    logObservableDensities :: c # f -> Observations f -> [Double]
    logObservableDensities p = map (logObservableDensity p)

    observableDensities :: c # f -> Observations f -> [Double]
    observableDensities p = map (exp . logObservableDensity p)


logObservableDensity :: ObservablyContinuous c f => c # f -> Observation f -> Double
logObservableDensity p = head . logObservableDensities p . (:[])


observableDensity :: ObservablyContinuous c f => c # f -> Observation f -> Double
observableDensity p = exp . head . logObservableDensities p . (:[])

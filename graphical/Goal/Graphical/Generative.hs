{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Graphical.Generative
    ( -- * Latent Variable Models
      ExpectationMaximization (expectationStep)
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


class Statistical x => ExpectationMaximization f z x where
    expectationStep :: Sample z -> Natural # f z x -> Mean # f z x

class Statistical x => ObservablyContinuous c f z x where
    logObservableDensity :: c # f z x -> SamplePoint z -> Double
    logObservableDensity p = head . logObservableDensities p . (:[])
    logObservableDensities :: c # f z x -> Sample z -> [Double]
    logObservableDensities p = map (logObservableDensity p)

    observableDensity :: c # f z x -> SamplePoint z -> Double
    observableDensity p = exp . head . logObservableDensities p . (:[])
    observableDensities :: c # f z x -> Sample z -> [Double]
    observableDensities p = map (exp . logObservableDensity p)



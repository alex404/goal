{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Graphical.Generative
    ( -- * Latent Variable Models
      ExpectationMaximization (expectationStep)
    -- ** Hierarchical Models
    , Observation
    , Observations
    , observableSample
    , observableSamplePoint
    , ObservablyContinuous
        ( logObservableDensity
        , logObservableDensities
        , observableDensity
        , observableDensities )
    ) where

--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry
import Goal.Probability


--- Latent Variable Class ---


class Statistical x => ExpectationMaximization c f z x where
    expectationStep :: Sample z -> c # f z x -> Mean # f z x

-- | A 'SamplePoint' construction for 'HList's.

type family HHead as where
    HHead (HList (a ': _)) = '[a]

type family Head as where
    Head (HList (a ': _)) = a

type Observation x = Head (SamplePoint x)

type Observations x = [Observation x]

observableSample
    :: (Generative c x, SamplePoint x ~ HList (a : as))
    => Int -> c # x -> Random r (Observations x)
observableSample nsmps p = map hHead <$> sample nsmps p

observableSamplePoint
    :: (Generative c x, SamplePoint x ~ HList (a : as))
    => c # x -> Random r (Observation x)
observableSamplePoint p = hHead <$> samplePoint p

class Statistical x => ObservablyContinuous c f z x where
    logObservableDensity :: c # f z x -> SamplePoint z -> Double
    logObservableDensity p = head . logObservableDensities p . (:[])
    logObservableDensities :: c # f z x -> Sample z -> [Double]
    logObservableDensities p = map (logObservableDensity p)

    observableDensity :: c # f z x -> SamplePoint z -> Double
    observableDensity p = exp . head . logObservableDensities p . (:[])
    observableDensities :: c # f z x -> Sample z -> [Double]
    observableDensities p = map (exp . logObservableDensity p)



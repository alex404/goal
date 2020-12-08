{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Conditional.Dynamic
    (
    stateSpaceTransition
    , sampleStateSpaceModel
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Conditional
import Goal.Graphical.Generative.Harmonium


--- Generic ---


stateSpaceTransition
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Generative Natural x, Generative Natural z
       , Bilinear g z x, Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> SamplePoint x
    -> Random s (SamplePoint (z,x))
stateSpaceTransition trns' emsn' x = do
    x' <- samplePoint $ trns' >.>* x
    z' <- samplePoint $ emsn' >.>* x'
    return (z',x')

sampleStateSpaceModel
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Generative Natural x, Generative Natural z
       , Bilinear g z x, Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> Int
    -> Natural # x
    -> Random s (Sample (z,x))
sampleStateSpaceModel trns' emsn' n prr' = do
    x0 <- samplePoint prr'
    z0 <- samplePoint $ emsn' >.>* x0
    iterateM n (stateSpaceTransition trns' emsn' . snd) (z0,x0)


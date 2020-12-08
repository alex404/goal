{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Conditional.Dynamic
    (
    stateSpaceTransition
    , sampleStateSpaceModel
    , stateSpaceDensity
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
stateSpaceTransition trns emsn x = do
    x' <- samplePoint $ trns >.>* x
    z' <- samplePoint $ emsn >.>* x'
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
sampleStateSpaceModel trns emsn n prr = do
    x0 <- samplePoint prr
    z0 <- samplePoint $ emsn >.>* x0
    iterateM n (stateSpaceTransition trns emsn . snd) (z0,x0)

stateSpaceDensity
    :: ( ExponentialFamily z, ExponentialFamily x, Map Natural f x x
       , Map Natural g z x, AbsolutelyContinuous Natural x
       , AbsolutelyContinuous Natural z  )
    => Natural # x
    -> Natural # Affine f x x
    -> Natural # Affine g z x
    -> Sample (z,x)
    -> Double
stateSpaceDensity prr trns emsn zxs =
    let (zs,xs) = unzip zxs
        prrdns = density prr $ head xs
        trnsdnss = zipWith density (trns >$>* xs) $ tail xs
        emsndnss = zipWith density (emsn >$>* xs) zs
     in product $ prrdns : trnsdnss ++ emsndnss

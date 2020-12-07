#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S

-- Unqualified --


trns :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
trns = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.1,0.1)
        , fromTuple (0.8,0.1)
        , fromTuple (0.1,0.8) )
    )  (fromTuple (0.33,0.33))

emsn :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
emsn = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.2,0.2)
        , fromTuple (0.6,0.2)
        , fromTuple (0.2,0.6) )
    )  (fromTuple (0.33,0.33))

prr :: Natural # Categorical 2
prr = toNatural (fromTuple (0.33,0.33) :: Mean # Categorical 2)


transition'
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Generative Natural x, Generative Natural z
       , Bilinear g z x, Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> SamplePoint x
    -> Random s (SamplePoint (z,x))
transition' trns' emsn' x = do
    x' <- samplePoint $ trns' >.>* x
    z' <- samplePoint $ emsn' >.>* x'
    return (z',x')

sampleHMM
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Generative Natural x, Generative Natural z
       , Bilinear g z x, Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> Int
    -> Natural # x
    -> Random s (Sample (z,x))
sampleHMM trns' emsn' n prr' = do
    x0 <- samplePoint prr'
    z0 <- samplePoint $ emsn' >.>* x0
    iterateM n (transition' trns' emsn' . snd) (z0,x0)


--- Main ---


main :: IO ()
main = do

    zxs <- realize $ sampleHMM trns emsn 20 prr
    putStrLn "HMM Simulation:"
    print zxs
    let xs = snd <$> zxs

    let flts = conjugatedFiltering trns emsn prr $ fst <$> zxs
    putStrLn "\nFiltering Probabilities:"
    print . average $ zipWith (!!) (S.toList . categoricalWeights <$> flts) xs

    let smths = conjugatedSmoothing trns emsn prr $ fst <$> zxs
    putStrLn "\nSmoothing Probabilities:"
    print . average $ zipWith (!!) (S.toList . categoricalWeights <$> smths) xs

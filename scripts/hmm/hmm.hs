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


-- True Model --

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

-- Learning

trns0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
trns0 = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.replicate (Point $ S.replicate 0.33)
    )  (fromTuple (0.33,0.33))

emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
emsn0 = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.replicate (Point $ S.replicate 0.33)
    )  (fromTuple (0.33,0.33))

prr0 :: Natural # Categorical 2
prr0 = toNatural (fromTuple (0.33,0.33) :: Mean # Categorical 2)

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)


--- Main ---


main :: IO ()
main = do

    zss <- realize . replicateM 1000 $ map fst <$> sampleStateSpaceModel trns emsn 20 prr

    let em (prr',trns',emsn') = stateSpaceExpectationMaximizationAscent
            alg alg prr' trns' emsn' zss

        hmms = take 10 $ iterate em (prr0,trns0,emsn0)

    putStrLn "True Model:"
    print (prr,trns,emsn)

    putStrLn "Learned Models:"
    print $ last hmms

    --putStrLn "HMM Simulation:"
    --print zxs
    --let xs = snd <$> zxs

    --let flts = conjugatedFiltering trns emsn prr $ fst <$> zxs
    --putStrLn "\nFiltering Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> flts) xs

    --let smths = conjugatedSmoothing prr trns emsn $ fst <$> zxs
    --putStrLn "\nSmoothing Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> smths) xs

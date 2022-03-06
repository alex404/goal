#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L

--- Globals ---


-- Simulation --

nsmps :: Int
nsmps = 20000

nepchs :: Int
nepchs = 20

-- True Model --

emsn :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
emsn = fst . split . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.2,0.2)
        , fromTuple (0.6,0.2)
        , fromTuple (0.2,0.6) )
    )  (fromTuple (0.33,0.33))

prr :: Natural # Categorical 2
prr = toNatural (fromTuple (0.15,0.35) :: Mean # Categorical 2)

hrm :: Natural # Harmonium Tensor (Categorical 2) (Categorical 2)
hrm = joinConjugatedHarmonium emsn prr

-- | A single iteration of EM for 'Harmonium' based models.
expectationMaximization' zs aff' prr' =
    let mxs = toMean . conjugatedBayesRule  aff' prr' <$> zs
        mx = average mxs
     in transition mx

expectationMaximization'' zs aff' prr' =
    let mhrm = expectationStep zs $ joinConjugatedHarmonium aff' prr'
        mprr = snd $ split mhrm
     in snd . splitConjugatedHarmonium $ transition mhrm


--- Main ---


main :: IO ()
main = do


    zxs <- realize $ sample nsmps hrm
    let zs,xs :: [Int]
        (zs,xs) = unzip zxs

    let mx :: Mean # Categorical 2
        mx = averageSufficientStatistic xs
    print . S.toList $ categoricalWeights mx

    prr0 :: Natural # Categorical 2 <- realize $ uniformInitialize (-0.001,0.001)

    let hrm0 = joinConjugatedHarmonium emsn prr0

    let em = expectationMaximization' zs emsn
    let em' = expectationMaximization'' zs emsn
    --let em = latentProcessExpectationMaximizationAscent 1e-3 500 defaultAdamPursuit zss

        prrs = take nepchs $ iterate em prr0
        prrs' = take nepchs $ iterate em' prr0

    let shower vls = L.intercalate ", " $ (\vl -> showFFloat (Just 3) vl "") <$> vls
    let eval prr' =
            let hrm' = joinConjugatedHarmonium emsn prr'
                ll = average $ logObservableDensities hrm' zs
            in ll-- , shower . S.toList $ categoricalWeights prr')

    putStrLn $ "\nTrue LL:" ++ show (eval prr)
    putStrLn "Gradient Ascent:"
    mapM_ print $ zip (eval <$> prrs) (eval <$> prrs')

    --putStrLn "\nModels:"
    --putStrLn "\nInitial:"
    --printHMM $ head hmms
    --putStrLn "\nTarget:"
    --printHMM hmm
    --putStrLn "\nLearned:"
    --printHMM $ last hmms


--- Graveyard ---


    --zs <- realize $ map snd <$> sampleStateSpaceModel trns emsn 2 prr

    --let smths = conjugatedSmoothing prr trns emsn zs
    --    msmths = toMean <$> smths
    --let (hdmsmths,tlmsmths) = unzip . zip msmths $ tail msmths
    --    mxx = tlmsmths >$< hdmsmths
    --    mx' = average tlmsmths
    --    mx = average hdmsmths
    --    trns' :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
    --    trns' = fst . splitBottomHarmonium . toNatural $ joinHarmonium mx' mxx mx

    --    mzs = sufficientStatistic <$> zs
    --    mzx = mzs >$< msmths
    --    emsn' :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
    --    emsn' = fst . splitBottomHarmonium . toNatural $ joinHarmonium (average mzs) mzx (average msmths)

    --mapM_ print $ S.toList . categoricalWeights <$> trns >$>* xspc
    --mapM_ print $ S.toList . categoricalWeights <$> emsn >$>* xspc
    --mapM_ print $ S.toList . categoricalWeights <$> trns' >$>* xspc
    --mapM_ print $ S.toList . categoricalWeights <$> emsn' >$>* xspc

    --let x = 0
    --let xwghts = (`density` x) <$> smths

    --print "foo"

    --let ln = 10

    --zs <- realize $ map fst <$> sampleLatentProcess ln hmm

    --let smths = fst $ conjugatedSmoothing trns emsn prr zs
    --putStrLn "\nSmoothing Probabilities:"
    --mapM_ print $ categoricalWeights <$> smths


    --putStrLn "HMM Simulation:"
    --print zxs
    --let xs = snd <$> zxs

    --let flts = conjugatedFiltering trns emsn prr $ fst <$> zxs
    --putStrLn "\nFiltering Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> flts) xs

    --let smths = conjugatedSmoothing prr trns emsn $ fst <$> zxs
    --putStrLn "\nSmoothing Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> smths) xs

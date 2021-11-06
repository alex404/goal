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


--- Globals ---


-- Simulation --

pthln :: Int
pthln = 8

nsmps :: Int
nsmps = 500

nepchs :: Int
nepchs = 500

-- True Model --

trns :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
trns = fst . split . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.2,0.2)
        , fromTuple (0.6,0.2)
        , fromTuple (0.2,0.6) )
    )  (fromTuple (0.33,0.33))

emsn :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
emsn = fst . split . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.2,0.2)
        , fromTuple (0.6,0.2)
        , fromTuple (0.2,0.6) )
    )  (fromTuple (0.33,0.33))

prr :: Natural # Categorical 2
prr = toNatural (fromTuple (0.2,0.35) :: Mean # Categorical 2)

hmm :: Natural # HiddenMarkovModel 2 2
hmm = joinLatentProcess prr emsn trns

-- Testing

bruteForceMarginalization :: [Int] -> (Int, Int) -> Double
bruteForceMarginalization zs (stp,x0) =
    let dnm = logSumExp $ logDensity hmm . zip zs <$> replicateM pthln xspc
        nmrsqs = do
            hds <- replicateM stp xspc
            tls <- replicateM (pthln - stp - 1) xspc
            return $ hds ++ [x0] ++ tls
        nmr = logSumExp $ logDensity hmm . zip zs <$> nmrsqs
     in exp $ nmr-dnm


-- Learning

printHMM
    :: Natural # HiddenMarkovModel 2 2
    -> IO ()
printHMM hmm' = do
    let (prr',emsn',trns') = splitLatentProcess hmm'
    putStrLn "Prior: "
    print . S.toList $ categoricalWeights prr'
    putStrLn "Transitions: "
    mapM_ print $ S.toList . categoricalWeights <$> trns' >$>* xspc
    putStrLn "Emissions: "
    mapM_ print $ S.toList . categoricalWeights <$> emsn' >$>* xspc

xspc :: [Int]
xspc = [0,1,2]


--- Main ---


main :: IO ()
main = do


    zss <- realize . replicateM nsmps $ map fst <$> sampleLatentProcess pthln hmm
    let zs = head zss

    putStrLn "\n-- Observable Densities --"

    putStrLn "\nAnalytic:"
    let smths = conjugatedSmoothing hmm zs
    mapM_ print $ S.toList . categoricalWeights . toMean <$> smths

    putStrLn "\nBrute Force:"
    mapM_ print [ [ bruteForceMarginalization zs (stp,x) | x <- [0,1,2]] | stp <- [0..pthln-1]]

    putStrLn "\n-- Learning --\n"

    trns0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-0.1,0.1)
    emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-0.1,0.1)
    prr0 :: Natural # Categorical 2 <- realize $ uniformInitialize (-0.1,0.1)

    let hmm0 = joinLatentProcess prr0 emsn0 trns0

    let em = latentProcessExpectationMaximization zss

        hmms = take nepchs $ iterate em hmm0

    let lls hmm' =
            let ll = average $ logObservableDensities hmm' zss
                (prr',_,_) = splitLatentProcess hmm'
            in (ll,show $ toMean prr')

    putStrLn $ "\nTrue LL:" ++ show (fst $ lls hmm)
    putStrLn "Gradient Ascent:"
    mapM_ (print . lls) hmms

    putStrLn "\nModels:"
    putStrLn "\nInitial:"
    printHMM $ head hmms
    putStrLn "\nTarget:"
    printHMM hmm
    putStrLn "\nLearned:"
    printHMM $ last hmms


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

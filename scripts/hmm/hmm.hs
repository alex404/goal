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


-- True Model --

trns :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
trns = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.1,0.2)
        , fromTuple (0.7,0.2)
        , fromTuple (0.2,0.7) )
    )  (fromTuple (0.33,0.33))

emsn :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
emsn = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.1,0.3)
        , fromTuple (0.6,0.1)
        , fromTuple (0.1,0.6) )
    )  (fromTuple (0.33,0.33))

prr :: Natural # Categorical 2
prr = toNatural (fromTuple (0.33,0.33) :: Mean # Categorical 2)

-- Learning

--trns0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
--trns0 = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
--    ( S.replicate (Point $ S.replicate 0.3)
--    )  (fromTuple (0.33,0.33))
--
--emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
--emsn0 = fst . splitBottomHarmonium . toNatural $ joinMeanMixture
--    ( S.replicate (Point $ S.replicate 0.3)
--    )  (fromTuple (0.33,0.33))
--
--prr0 :: Natural # Categorical 2
--prr0 = toNatural (fromTuple (0.33,0.33) :: Mean # Categorical 2)

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)

printHMM (prr',trns',emsn') = do
    putStrLn "Prior: "
    print . S.toList $ categoricalWeights prr'
    putStrLn "Transitions: "
    mapM_ print $ S.toList . categoricalWeights <$> trns' >$>* xspc
    putStrLn "Emissions: "
    mapM_ print $ S.toList . categoricalWeights <$> emsn' >$>* xspc

xspc :: [Int]
xspc = [0,1,2]

bruteForceMarginalization ln zs (stp,x0) =
    let dnm = logSumExp $ stateSpaceLogDensity prr trns emsn . zip zs <$> replicateM ln xspc
        nmrsqs = do
            hds <- replicateM stp xspc
            tls <- replicateM (ln - stp - 1) xspc
            return $ hds ++ [x0] ++ tls
        nmr = logSumExp $ stateSpaceLogDensity prr trns emsn . zip zs <$> nmrsqs
     in exp$nmr-dnm



--- Main ---


main :: IO ()
main = do

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

    --zxs <- realize $ sampleStateSpaceModel trns emsn ln prr

    --putStrLn "Simulation: "
    --print zxs

    --let (zs,xs) = unzip zxs
    --print zxs

    --let (_,smths,_) = unzip3 $ conjugatedSmoothing prr trns emsn zs
    --putStrLn "\nSmoothing Probabilities:"
    --mapM_ print $ categoricalWeights <$> smths

    --let smths' = snd $ conjugatedSmoothing' prr trns emsn zs
    --putStrLn "\nSmoothing Probabilities':"
    --mapM_ print $ categoricalWeights <$> smths'



    --let ixs = zip [0..] xs
    --putStrLn "\nBrute Force:"
    --mapM_ print [ [ bruteForceMarginalization ln zs (stp,x) | x <- [0,1,2]] | stp <- [0..ln-1]]

    trns0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-1,1)
    emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-1,1)
    prr0 :: Natural # Categorical 2 <- realize $ uniformInitialize (-1,1)

    zss <- realize . replicateM 100 $ map fst <$> sampleStateSpaceModel trns emsn 20 prr

    let em (prr',trns',emsn') = stateSpaceExpectationMaximization' prr' trns' emsn' zss

        hmms = take 100 $ iterate em (prr0,trns0,emsn0)

    putStrLn "True Model:"
    printHMM (prr,trns,emsn)

    let lls (prr',trns',emsn') =
            average $ conjugatedFilteringLogDensity trns' emsn' prr' <$> zss

    mapM_ (print . lls) hmms

    putStrLn "\nModels:"
    putStrLn "\nInitial:"
    printHMM $ head hmms
    putStrLn "\nLearned:"
    printHMM $ last hmms

    --putStrLn "HMM Simulation:"
    --print zxs
    --let xs = snd <$> zxs

    --let flts = conjugatedFiltering trns emsn prr $ fst <$> zxs
    --putStrLn "\nFiltering Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> flts) xs

    --let smths = conjugatedSmoothing prr trns emsn $ fst <$> zxs
    --putStrLn "\nSmoothing Probabilities:"
    --print . average $ zipWith (!!) (S.toList . categoricalWeights <$> smths) xs

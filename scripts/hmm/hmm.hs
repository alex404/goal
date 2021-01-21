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

trns :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
trns = fst . split . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.1,0.2)
        , fromTuple (0.7,0.2)
        , fromTuple (0.2,0.7) )
    )  (fromTuple (0.33,0.33))

emsn :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
emsn = fst . split . toNatural $ joinMeanMixture
    ( S.fromTuple
        ( fromTuple (0.1,0.3)
        , fromTuple (0.6,0.1)
        , fromTuple (0.1,0.6) )
    )  (fromTuple (0.33,0.33))

prr :: Natural # Categorical 2
prr = toNatural (fromTuple (0.33,0.33) :: Mean # Categorical 2)

ltnt :: Natural # LatentProcess Tensor Tensor (Categorical 2) (Categorical 2) (Categorical 2) (Categorical 2)
ltnt = joinLatentProcess prr emsn trns

-- Learning

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)

printHMM
    :: Natural # LatentProcess Tensor Tensor (Categorical 2) (Categorical 2) (Categorical 2) (Categorical 2)
    -> IO ()
printHMM ltnt' = do
    let (prr',emsn',trns') = splitLatentProcess ltnt'
    putStrLn "Prior: "
    print . S.toList $ categoricalWeights prr'
    putStrLn "Transitions: "
    mapM_ print $ S.toList . categoricalWeights <$> trns' >$>* xspc
    putStrLn "Emissions: "
    mapM_ print $ S.toList . categoricalWeights <$> emsn' >$>* xspc

xspc :: [Int]
xspc = [0,1,2]

bruteForceMarginalization :: Int -> [Int] -> (Int, Int) -> Double
bruteForceMarginalization ln zs (stp,x0) =
    let dnm = logSumExp $ logDensity ltnt . zip zs <$> replicateM ln xspc
        nmrsqs = do
            hds <- replicateM stp xspc
            tls <- replicateM (ln - stp - 1) xspc
            return $ hds ++ [x0] ++ tls
        nmr = logSumExp $ logDensity ltnt . zip zs <$> nmrsqs
     in exp $ nmr-dnm



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

    --zs <- realize $ map fst <$> sampleLatentProcess ln ltnt

    --let smths = fst $ conjugatedSmoothing trns emsn prr zs
    --putStrLn "\nSmoothing Probabilities:"
    --mapM_ print $ categoricalWeights <$> smths

    --putStrLn "\nBrute Force:"
    --mapM_ print [ [ bruteForceMarginalization ln zs (stp,x) | x <- [0,1,2]] | stp <- [0..ln-1]]

    trns0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-1,1)
    emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 2) (Categorical 2)
        <- realize $ uniformInitialize (-1,1)
    prr0 :: Natural # Categorical 2 <- realize $ uniformInitialize (-1,1)

    let ltnt0 = joinLatentProcess prr emsn trns

    zss <- realize . replicateM 200 $ map fst <$> sampleLatentProcess 200 ltnt

    let em = latentProcessExpectationMaximization zss

        hmms = take 1 $ iterate em ltnt0

    putStrLn "True Model:"
    printHMM ltnt

    let lls ltnt' = average $ logObservableDensities ltnt' zss

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

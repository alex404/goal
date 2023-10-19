#! stack runghc

{-# LANGUAGE ScopedTypeVariables,TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds,KindSignatures,TypeFamilies,NoStarIsType,UndecidableInstances,TypeSynonymInstances,FlexibleInstances,MultiParamTypeClasses
   #-}

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G



--- Globals ---


-- Reproduced tutorial here: https://www.math.umd.edu/~petersd/666/html/iris_pca.html#1

ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/irisdata.dat"

-- Training --

-- Functions --

parseCSV :: String -> [(S.Vector 4 Double,Int)]
parseCSV csvstr = do
    csv <- dropWhile ((== '%') . head) $ lines csvstr
    let x0 = read $ '[' : csv ++ "]"
        (x,cats) = splitAt 4 x0
    return (fromJust $ S.fromList x, round $ head cats)

parseHMoG :: String -> Natural # IsotropicHMoG 4 2 2
parseHMoG hmogstr =
    let cs :: S.Vector (Dimension (IsotropicHMoG 4 2 2)) Double
        cs = fromJust . S.fromList $ read hmogstr
     in Point cs

--printer :: [S.Vector 4 Double] -> (Int,Natural # IsotropicHMoG 4 2 2) -> IO ()
--printer smps (k,nhmog) = do
--    putStrLn $ "\nIteration: " ++ show k
--    let mhmog = expectationStep smps nhmog
--    print $ logLikelihood smps nhmog
--    print $ nhmog
--    print $ mhmog
--    print $ toSource mhmog
--    print $ toNatural mhmog

--tester :: [S.Vector 4 Double] -> (Int,Natural # IsotropicHMoG 4 2 2) -> IO ()
--tester smps (k,nhmog) = do
--    putStrLn $ "\nIteration: " ++ show k
--    let (lkl,nmog) = splitConjugatedHarmonium nhmog
--    yzs <- realize $ sample 100000 nmog
--    let mmog' :: Mean # Mixture (FullNormal 2) 2
--        mmog' = averageSufficientStatistic yzs
--    print $ euclideanDistance nmog (toNatural $ toMean nmog)
--    print $ euclideanDistance mmog' (toMean nmog)
--    let (nx,nxy,nyhrm) = splitHarmonium nhmog
--        (ny,nyz,nz) = splitHarmonium nyhrm
--        lgh = joinHarmonium nx nxy ny
--
--    xyzs <- realize $ sample 100000 nhmog
--    let mhmog' :: Mean # IsotropicHMoG 4 2 2
--        mhmog' = averageSufficientStatistic xyzs
--    print $ listCoordinates nhmog
--    print $ listCoordinates (toNatural $ toMean nhmog)
--    print $ euclideanDistance mhmog' (toMean nhmog)
--    print $ euclideanDistance mhmog' (toMean $ toNatural mhmog')
--    print $ (\val -> showFFloat (Just 4) val "") . square <$> zipWith (-) (listCoordinates nhmog) (listCoordinates . toNatural $ toMean nhmog)

toMatrix :: (Manifold x, Manifold y) => c # Tensor y x -> S.Matrix (Dimension y) (Dimension x) Double
toMatrix (Point xs) = G.Matrix xs

affineMixtureToMixture
    :: (KnownNat k, Manifold x0, Manifold x, Translation x x0)
    => Natural # AffineMixture x0 x k
    -> Natural # Mixture x k
affineMixtureToMixture lmxmdl =
    let (flsk,nk) = split lmxmdl
        (nls,nlk) = split flsk
        nlsk = fromColumns . S.map (0 >+>) $ toColumns nlk
     in join (join nls nlsk) nk

--mixtureToAffineMixture
--    :: (KnownNat k, Manifold x, Manifold x0, Translation x x0)
--    => Mean # Mixture x k
--    -> Mean # AffineMixture x0 x k
--mixtureToAffineMixture mxmdl =
--    let (flsk,mk) = split mxmdl
--        (mls,mlsk) = split flsk
--        mlk = fromColumns . S.map anchor $ toColumns mlsk
--     in join (join mls mlk) mk


toBigMatrix
    :: ( Bilinear c MVNCovariance (MVNMean 2) (MVNMean 2) )
    => c # IsotropicGaussianHarmonium 4 2
    -> S.Matrix 6 6 Double
toBigMatrix lgh =
    let (x,tr0,z) = splitHarmonium lgh
        tr = 2/> tr0
        (_,tl) = split x
        (_,br) = split z
        top = S.horizontalConcat (toMatrix $ toTensor tl) (toMatrix tr)
        btm = S.horizontalConcat (S.transpose $ toMatrix tr) (toMatrix $ toTensor br)
     in S.verticalConcat top btm

type HMoG2 f n m k = AffineMixture (FullNormal m) (LinearGaussianHarmonium f n m) k

hmog1to2
    :: ( KnownNat n, KnownNat m, KnownNat k
       , Manifold (f (MVNMean n) (MVNMean n)) )
    => c # HierarchicalMixtureOfGaussians f n m k
    -> c # HMoG2 f n m k
hmog1to2 hmog =
    let (lkl,mog) = split hmog
        (aff,cats) = split mog
        (mvn,tns) = split aff
     in join (join (lkl `join` mvn) tns) cats



covarianceTester :: [S.Vector 4 Double] -> (Int,Natural # IsotropicHMoG 4 2 2) -> IO ()
covarianceTester _ (_,nhmog) = do
    let (lghs0,nk) = splitNaturalMixture . affineMixtureToMixture $ hmog1to2 nhmog
        (mx,_,_) = splitHarmonium $ toMean nhmog
        mxmu = fst $ split mx
        lghs = S.toList lghs0
        --mtxs = toBigMatrix . toMean <$> lghs
        wghts = S.toList $ categoricalWeights nk
        --mtx = sum $ zipWith (\w mt -> S.withMatrix (S.scale w) mt) wghts mtxs
        --(mlgh,_,_) = splitHarmonium . hmog1to2 $ toMean nhmog
        --ntxs = toBigMatrix <$> lghs
        --ntx = sum $ zipWith (\w mt -> S.withMatrix (S.scale w) mt) wghts ntxs
    --print

    putStrLn "Natural Mus:"

    print $ do
        lgh <- lghs
        let (nx,_,_) = splitHarmonium lgh
        return . fst $ split nx

    putStrLn "Natural to Mean Mus:"

    print $ do
        lgh <- lghs
        let (nx,_,nz) = splitHarmonium lgh
            nxmu = fst $ split nx
            nzmu = fst $ split nz
            ntx = toBigMatrix lgh
            invntx = -0.5 * S.inverse ntx
            nxz' = S.matrixVectorMultiply invntx $ coordinates nxmu S.++ coordinates nzmu
            mx' :: S.Vector 4 Double
            mx' = S.take nxz'
        return mx'

    putStrLn "Mean Mus:"

    print $ do
        lgh <- toMean <$> lghs
        let (mx',_,_) = splitHarmonium lgh
        return . fst $ split mx'
    --print $ toBigMatrix mlgh

    putStrLn "Mean Mu:"

    print mxmu
    print . sum . zipWith (.>) wghts $ do
        lgh <- toMean <$> lghs
        let (mx',_,_) = splitHarmonium lgh
        return . fst $ split mx'
    --print $ toBigMatrix mlgh

    --let avginv = sum . zipWith (.>) wghts $ do
    --    lgh <- lghs
    --    let (nx,_,nz) = splitHarmonium lgh
    --        nxmu = fst $ split nx
    --        nzmu = fst $ split nz
    --        ntx = toBigMatrix lgh
    --        invntx = -0.5 * S.inverse ntx
    --    return invntx

    --print $ avginv >.>


    --let bigMatrix,bigMatrix2 :: S.Matrix (M+N) (M+N) Double
    --    bigMatrix =
    --    bigMatrix2 = G.Matrix cs4


    --print $ listCoordinates lgh
    --print $ listCoordinates (toNatural $ toMean lgh)
    --print $ euclideanDistance lgh (toNatural $ toMean lgh)
    --print . length $ listCoordinates lgh
    --print . length $ listCoordinates nmog

--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    hmogstr <- readFile ("badhmogs/bad-hmog-0")

    let smpcats = parseCSV csvstr
        n = length smpcats
        smps = fst <$> smpcats

    print $ "Number of samples:" ++ show n

    let nhmog = parseHMoG hmogstr
        knhmogs = take 1 . zip [0..] $ iterate (expectationMaximization smps) nhmog
    mapM_ (covarianceTester smps) knhmogs




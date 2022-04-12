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
import qualified Data.List as L

import qualified Data.Map as M


--- Globals ---


-- Reproduced tutorial here: https://www.math.umd.edu/~petersd/666/html/iris_pca.html#1

ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/irisdata.dat"

-- Training --

nepchs :: Int
nepchs = 200

-- Functions --

parseCSV :: String -> [(S.Vector 4 Double,Int)]
parseCSV csvstr = do
    csv <- dropWhile ((== '%') . head) $ lines csvstr
    let x0 = read $ '[' : csv ++ "]"
        (x,cats) = splitAt 4 x0
    return (fromJust $ S.fromList x, round $ head cats)

-- | This is a selfcontained implementation of PCA-EM from Bishop, against
-- which I can benchmark the EF-based version.

type StandardPCA n k = S.Vector (1 + n + n*k) Double

joinStandardPCA
    :: (KnownNat n, KnownNat k)
    => S.Vector n Double -- ^ Mean bias
    -> Double -- ^ Variance
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> StandardPCA n k
joinStandardPCA mus vr mtx =
    mus S.++ S.singleton vr S.++ G.toVector mtx

splitStandardPCA
    :: (KnownNat n, KnownNat k)
    => StandardPCA n k
    -> (S.Vector n Double, Double, S.Matrix n k Double)
splitStandardPCA cs =
    let (mus,cs') = S.splitAt cs
        (vr,mtx) = S.splitAt cs'
     in (mus,S.head vr,G.Matrix mtx)

standardPCAToMultivariateNormal
    :: (KnownNat n, KnownNat k)
    => StandardPCA n k
    -> Source # FullNormal n
standardPCAToMultivariateNormal pca =
    let (mus,vr,mtx) = splitStandardPCA pca
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = realToFrac vr * S.matrixIdentity
     in join (Point mus) . fromTensor . Point . G.toVector $ mtx1 + mtx2

standardPCAExpectationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> StandardPCA n k
    -> ([S.Vector k Double], [S.Matrix k k Double])
standardPCAExpectationStep xs pca =
    let (_,vr,wmtx) = splitStandardPCA pca
        wmtxtr = S.transpose wmtx
        vrmtx = realToFrac vr * S.matrixIdentity
        wmlts = S.matrixMatrixMultiply wmtxtr wmtx
        mmtxinv = S.inverse $ wmlts + vrmtx
        rsds = [ x - average xs | x <- xs ]
        ezs = S.matrixVectorMultiply (S.matrixMatrixMultiply mmtxinv wmtxtr) <$> rsds
        ezzs = [ realToFrac vr * mmtxinv + S.outerProduct ez ez | ez <- ezs ]
     in (ezs,ezzs)

standardPCAMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> [S.Vector k Double]
    -> [S.Matrix k k Double]
    -> StandardPCA n k
standardPCAMaximizationStep xs ezs ezzs =
    let xht = average xs
        rsds = [ x - xht | x <- xs ]
        invsgm = S.pseudoInverse $ sum ezzs
        wmtx0 = sum $ zipWith S.outerProduct rsds ezs
        wmtx = S.matrixMatrixMultiply wmtx0 invsgm
        n = fromIntegral $ natVal (Proxy @n)
        vr = (recip n *) . average $ do
            (ez,ezz,rsd) <- zip3 ezs ezzs rsds
            let rsd2 = S.dotProduct rsd rsd
                intr = -2 * S.dotProduct rsd (S.matrixVectorMultiply wmtx ez)
                trc = S.trace $ S.matrixMatrixMultiply ezz
                       (S.matrixMatrixMultiply (S.transpose wmtx) wmtx)
            return $ rsd2 + intr + trc
     in joinStandardPCA xht vr wmtx

standardPCAExpectationMaximization
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> StandardPCA n k
    -> StandardPCA n k
standardPCAExpectationMaximization xs pca =
    let (ezs,ezzs) = standardPCAExpectationStep xs pca
     in standardPCAMaximizationStep xs ezs ezzs

standardPCAPosteriors
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> StandardPCA n k
    -> [Source # FullNormal k]
standardPCAPosteriors xs pca =
    let (mu,vr,wmtx) = splitStandardPCA pca
        wmtxtr = S.transpose wmtx
        vrmtx = realToFrac vr * S.matrixIdentity
        wmlts = S.matrixMatrixMultiply wmtxtr wmtx
        mmtx = wmlts + vrmtx
        mmtxinv = S.inverse mmtx
        rsds = [ x - mu | x <- xs ]
        muzs = S.matrixVectorMultiply (S.matrixMatrixMultiply mmtxinv wmtxtr) <$> rsds
        cvrz = fromTensor . Point . G.toVector $ realToFrac (recip vr) * mmtx
     in (`join` cvrz) . Point <$> muzs

-- Tests

standardToNaturalPCA
    :: (KnownNat n, KnownNat k)
    => StandardPCA n k
    -> Natural # PrincipleComponentAnalysis n k
standardToNaturalPCA sfa =
    let (cmu,cvr,cwmtx) = splitStandardPCA sfa
        invsg = recip cvr
        thtmu = Point $ realToFrac invsg * cmu
        thtsg = singleton $ (-0.5) * invsg
        imtx = fromTensor . Point . G.toVector $ realToFrac invsg * cwmtx
     in join (join thtmu thtsg) imtx

vr0 :: Double
vr0 = 10

prr0,prr1,prr2 :: Source # FullNormal 2
prr0 = fromTuple (0,0,vr0,0,vr0)
prr1 = fromTuple (0,0,vr0,-0.01,vr0)
prr2 = fromTuple (0,0,vr0,0.01,vr0)

nprr0,nprr1,nprr2 :: Natural # FullNormal 2
nprr0 = transition prr0
nprr1 = transition prr1
nprr2 = transition prr2

mog0 :: Natural # Mixture (FullNormal 2) 2
mog0 = joinNaturalMixture (S.fromTuple (nprr0,nprr1,nprr2)) 0

mog1 :: Natural # Mixture (FullNormal 2) 0
mog1 = breakManifold nprr0

eps :: Double
eps = 1e-3

nstps :: Int
nstps = 2000

partialExpectationStep
    :: forall f x x0 z z0 . ( ExponentialFamily x, LegendreExponentialFamily z
       , Translation x x0, Translation z z0
       , Bilinear Mean f x0 z0, Bilinear Natural f x0 z0 )
    => Sample x -- ^ Model Samples
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> [Natural # z] -- ^ Harmonium expected sufficient statistics
partialExpectationStep xs hrm =
    let mxs :: [Mean # x]
        mxs = sufficientStatistic <$> xs
        mx0s = anchor <$> mxs
        pstr = fst . split $ transposeHarmonium hrm
     in pstr >$> mx0s

correctedExpectationMaximization
    :: forall f x x0 z z0
     . ( ExponentialFamily x, LegendreExponentialFamily z, Translation x x0
       , Translation z z0, DuallyFlatExponentialFamily (AffineHarmonium f x0 z0 x z)
       , Bilinear Mean f x0 z0, Bilinear Natural f x0 z0, ConjugatedLikelihood f x0 z0 x z
       , Transition Mean Natural z )
    => Sample x -- ^ Model Samples
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
correctedExpectationMaximization xs hrm =
    let (lkl1,prr1) = splitConjugatedHarmonium $ expectationMaximization xs hrm
        bnd = 1e-2
        mprr1 :: Mean # z
        mprr1 = Point . S.map (max bnd) . S.map (min (1-bnd)) . coordinates $ toMean prr1
     in joinConjugatedHarmonium lkl1 $ transition mprr1


correctedExpectationMaximization'
    :: ( ExponentialFamily (MultivariateNormal f 4), Square Natural f (MVNMean 4)
       , Square Source f (MVNMean 4), Transition Source Mean (MultivariateNormal f 4)
       , Transition Mean Source (MultivariateNormal f 4)
       , LinearlyComposable Tensor f (MVNMean 4) (MVNMean 4) (MVNMean 4)
       , LinearlyComposable f Tensor (MVNMean 4) (MVNMean 4) (MVNMean 4)
       , LinearlyComposable f Tensor (MVNMean 4) (MVNMean 4) (MVNMean 2) )
    => [S.Vector 4 Double] -- ^ Model Samples
    -> Natural # HierarchicalMixtureOfGaussians f 4 2 2 -- ^ Harmonium
    -> Natural # HierarchicalMixtureOfGaussians f 4 2 2 -- ^ Harmonium
correctedExpectationMaximization' xs hrm =
    let (lkl1,mog1) = splitConjugatedHarmonium $ expectationMaximization xs hrm
        (lkl2, prr1) = splitConjugatedHarmonium mog1
        bnd = 1e-2
        mprr1 :: Mean # Categorical 2
        mprr1 = Point . S.map (max bnd) . S.map (min (1-bnd)) . coordinates $ toMean prr1
     in joinConjugatedHarmonium lkl1 . joinConjugatedHarmonium lkl2 $ transition mprr1

--correctedExpectationMaximization' xs hrm =
--    let mhrm0 = expectationStep xs hrm
--        bnd = 1e-3
--        (mlkl1,mmog0) = split mhrm0
--        (mlkl2,mprr0) = split mmog0
--        mprr = S.map (max bnd) $ S.map (min (1-bnd)) $ coordinates mprr0
--        mhrm = join mlkl $ Point mprr
--     in transition mhrm


tester :: [S.Vector 2 Double] -> Int -> Natural # Mixture (FullNormal 2) 2 -> String
tester smps k nhrm =
    let mhrm = expectationStep smps nhrm
        nhrm' = toNatural mhrm
        lngth = sqrt . sum . map square . listCoordinates
        ncmpts = concat . map listCoordinates . S.toList . fst $ splitNaturalMixture nhrm
        scmpts = concat . map listCoordinates . S.toList . fst . splitSourceMixture $ transition nhrm
        nwghts = listCoordinates . snd $ splitNaturalMixture nhrm
        mwghts = listCoordinates . snd $ splitMeanMixture mhrm
        showf x = showFFloat (Just 3) x ", "
        nzs = partialExpectationStep smps nhrm
        mzs = toMean <$> nzs
     in concat
         [ "Iteration: ", show k, "\n"
         , "\nlog-likelihood: ", showf $ logLikelihood smps nhrm
         , "\nNaturals:\n", unlines $ zipWith (\n str -> concat [show n, ": ", str]) [0..] $ concat . map showf . listCoordinates <$> nzs
         , "\nMeans:\n", unlines $ zipWith (\n str -> concat [show n, ": ", str]) [0..] $ concat . map showf . listCoordinates <$> mzs ]
         --, "\nNatural params: ", concat $ showf <$> listCoordinates nhrm
         --, "\nMean params: ", concat $ showf <$> listCoordinates mhrm ]
         --, "\nNatural components: ", concat $ showf <$> ncmpts
         --, "\nStandard components: ", concat $ showf <$> scmpts ]
         --, "\nExpectations: ", concat $ showf <$> listCoordinates mhrm, "\n"
         --, showFFloat (Just 3) (lngth mhrm) ", "
         --, showFFloat (Just 3) (lngth nhrm') "" ]
--    map (sqrt . sum . map square . listCoordinates) . S.toList . fst . splitNaturalMixture . expectationMaximization smps
--
--ighEM
--    :: [S.Vector 4 Double]
--    -> Natural # IsotropicGaussianHarmonium 4 2
--    -> Natural # IsotropicGaussianHarmonium 4 2
--ighEM zs igh0 =
--    let igh1 = expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps
--        pca = fst $ split igh1
--     in joinConjugatedHarmonium pca nprr0
--
--hmogEM
--    :: [S.Vector 4 Double]
--    -> Natural # IsotropicHMOG 4 2 2
--    -> Natural # IsotropicHMOG 4 2 2
--hmogEM zs igh0 =
--    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps
--
--dghEM
--    :: [S.Vector 4 Double]
--    -> Natural # DiagonalGaussianHarmonium 4 2
--    -> Natural # DiagonalGaussianHarmonium 4 2
--dghEM zs igh0 =
--    let igh1 = expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps
--        pca = fst $ split igh1
--     in joinConjugatedHarmonium pca nprr0
--
--dmogEM
--    :: [S.Vector 4 Double]
--    -> Natural # DiagonalHMOG 4 2 2
--    -> Natural # DiagonalHMOG 4 2 2
--dmogEM zs igh0 =
--    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps

--hmogConj
--    :: Natural # IsotropicHMOG2 4 2 2 -> (Double,Natural # Categorical 2)
--hmogConj hrm = conjugationParameters (fst $ split hrm)

--expectationMaximizationAscent'
--    :: (KnownNat n, KnownNat m, KnownNat k)
--    => Double
--    -> GradientPursuit
--    -> Sample (IsotropicNormal n)
--    -> Natural # IsotropicHMOG n m k
--    -> [Natural # IsotropicHMOG n m k]
--expectationMaximizationAscent' eps gp zs nhrm =
--    let mhrm' = expectationStep zs nhrm
--     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm

hmogEM xs hmog0 =
    expectationMaximizationAscent eps defaultAdamPursuit xs hmog0 !! nstps



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth

    let smpcats = parseCSV csvstr
        n = length smpcats
    print $ "Number of samples:" ++ show n

    (tsmpcats,vsmpcats) <- realize $ splitAt (round $ 0.8 * fromIntegral n) <$> shuffleList smpcats
    let (tsmps,tcats) = unzip tsmpcats
        (vsmps,vcats) = unzip vsmpcats

    let nvx :: Natural # IsotropicNormal 4
        nvx = toNatural $ averageSufficientStatistic tsmps
        svx = toSource nvx

    let nvx' :: Natural # DiagonalNormal 4
        nvx' = toNatural $ averageSufficientStatistic tsmps

        (mus,vr) = split svx

    (lds0 :: Cartesian # Tensor (Replicated 4 NormalMean) (Replicated 2 NormalMean))
        <- realize $ uniformInitialize (-1,1)

    let spca0 :: StandardPCA 4 2
        spca0 = joinStandardPCA (coordinates mus) (S.head $ coordinates vr) (G.Matrix $ coordinates lds0)
        npca0 :: Natural # PrincipleComponentAnalysis 4 2
        npca0 = standardToNaturalPCA spca0
        nfa0 :: Natural # FactorAnalysis 4 2
        nfa0 = join nvx' . snd $ split npca0
        igh0 = joinConjugatedHarmonium npca0 $ transition prr0
        --hmog0 :: Natural # IsotropicHMOG 4 2 2
        --hmog0 = joinConjugatedHarmonium npca0 mog0
    --    hmog1 :: Natural # IsotropicHMOG 4 2 0
    --    hmog1 = joinConjugatedHarmonium npca0 hmogprr1

    --print $ toMean igh0
    --print $ toMean hmog1
    --print $ toMean hmog0
    --smps' <- realize $ sample 100000 igh0
    --let igh0' :: Mean # IsotropicGaussianHarmonium 4 2
    --    igh0' = averageSufficientStatistic smps'
    --print igh0'
    --smps'' <- realize $ sample 100000 hmog0
    --print igh0'
    --let hmog0' :: Mean # IsotropicHMOG 4 2 2
    --    hmog0' = averageSufficientStatistic smps''
    --print hmog0'
    --print $ euclideanDistance (toMean igh0) igh0'
    --print $ euclideanDistance (toMean hmog0) hmog0'

    let emnpcas = take nepchs $ iterate (linearModelExpectationMaximization tsmps) npca0
        emnfas = take nepchs $ iterate (linearModelExpectationMaximization tsmps) nfa0
    let emspcas = take nepchs $ iterate (standardPCAExpectationMaximization tsmps) spca0
    let snrms :: [Natural # FullNormal 4]
        snrms = toNatural . standardPCAToMultivariateNormal <$> emspcas
    --let emighs = take nepchs $ iterate (ighEM smps) igh0
    --let emhmogs = take nepchs $ iterate (hmogEM smps) hmog0

    let lls = zip3
            (logLikelihood vsmps <$> snrms)
            (logLikelihood vsmps . linearModelObservableDistribution <$> emnpcas)
            (logLikelihood vsmps . linearModelObservableDistribution <$> emnfas)

    putStrLn "Bishop/Natural PCA/FA LL Ascent:"
    --mapM_ print lls

    let npca1 = last emnpcas
        nfa1 = last emnfas
        prjctn1 :: Natural # Affine Tensor (MVNMean 2) (FullNormal 2) (MVNMean 4)
        prjctn1 = fst . split . transposeHarmonium $ joinConjugatedHarmonium npca1 standardNormal
        tprjcts1 = coordinates . fst . split . toSource <$> prjctn1 >$>* tsmps
        vprjcts1 = coordinates . fst . split . toSource <$> prjctn1 >$>* vsmps
        prjctn1' :: Natural # Affine Tensor (MVNMean 2) (FullNormal 2) (MVNMean 4)
        prjctn1' = fst . split . transposeHarmonium $ joinConjugatedHarmonium nfa1 standardNormal
        tprjcts1' = coordinates . fst . split . toSource <$> prjctn1' >$>* tsmps
        vprjcts1' = coordinates . fst . split . toSource <$> prjctn1' >$>* vsmps

    let mogs = take nepchs $ iterate (correctedExpectationMaximization tprjcts1) mog0
        mogs' = take nepchs $ iterate (correctedExpectationMaximization tprjcts1') mog0

    let moglls = logLikelihood vprjcts1 <$> mogs
        moglls' = logLikelihood vprjcts1' <$> mogs'
        mogdts = zipWith (tester tprjcts1) [0..] mogs
        mogdts' = zipWith (tester tprjcts1') [0..] mogs'
        hmoglls = logLikelihood tsmps . joinConjugatedHarmonium npca1 <$> mogs
        hmoglls' = logLikelihood tsmps . joinConjugatedHarmonium nfa1 <$> mogs'
    putStrLn "Mog vs HMog LL Ascent:"
    --mapM_ putStrLn mogdts'
    mapM_ print $ L.zip4 moglls hmoglls moglls' hmoglls'
    --mapM_ print $ zip3 [0..] moglls hmoglls
    --mapM_ putStrLn mogdts

    let mog1 = last mogs
        mog1' = last mogs'
        hmog1 = joinConjugatedHarmonium npca0 mog0
        hmog1' = joinConjugatedHarmonium nfa0 mog0

    let emhmogs = take nepchs . iterate (hmogEM tsmps)
            $ joinConjugatedHarmonium npca0 mog0
    let emhmogs' = take nepchs . iterate (hmogEM tsmps)
            $ joinConjugatedHarmonium nfa0 mog0
        hmoglls1 = logLikelihood tsmps <$> emhmogs
        hmoglls1' = logLikelihood tsmps <$> emhmogs'

    putStrLn "True HMog LL Ascent:"
    mapM_ print $ zip hmoglls1 hmoglls1'

    putStrLn "Information Gains:"
    print (last hmoglls1 - last hmoglls,last hmoglls1' - last hmoglls')


    --let hmog2 = last emhmogs
    --    (npca2,mog2) = splitConjugatedHarmonium hmog2
    --let hmog2' = last emhmogs'
    --    (nfa2,mog2') = splitConjugatedHarmonium hmog2'

    --    prjctn2 :: Natural # Affine Tensor (MVNMean 2) (MultivariateNormal 2) (MVNMean 4)
    --    prjctn2 = fst . split . transposeHarmonium $ joinConjugatedHarmonium npca2 $ standardNormal
    --    tprjcts2 = fst . splitMultivariateNormal . toSource <$> prjctn2 >$>* tsmps

    --    prjctn2' :: Natural # Affine Tensor (MVNMean 2) (MultivariateNormal 2) (MVNMean 4)
    --    prjctn2' = fst . split . transposeHarmonium . joinConjugatedHarmonium nfa2 $ standardNormal
    --    tprjcts2' = fst . splitMultivariateNormal . toSource <$> prjctn2' >$>* tsmps


    --let prjctss = [tprjcts1,tprjcts1',tprjcts2,tprjcts2']
    --    clstrmp1 = fst $ split mog1
    --    clstrmp2 = fst $ split mog2
    --    clstrmp1' = fst $ split mog1'
    --    clstrmp2' = fst $ split mog2'
    --    clstrmps = [clstrmp1,clstrmp1',clstrmp2,clstrmp2']
    --    nms = ["standard-pca","standard-fa","joint-pca","joint-fa"]

    --forM_ (zip3 prjctss clstrmps nms) $ \(prjcts,clstrmp,nm) -> do

    --    let catmp :: M.Map Int [S.Vector 2 Double]
    --        catmp = M.fromListWith (++) $ zip tcats $ (:[]) <$> prjcts
    --        kys = M.keys catmp

    --    sequence_ $ do
    --        ky <- kys
    --        let pnts = S.toList <$> catmp M.! ky
    --            mvn = toSource $ clstrmp >.>* ky
    --            xys = bivariateNormalConfidenceEllipse 100 1 mvn
    --        return $ do
    --            goalExport ldpth ("projection-cat-" ++ show ky) $ pnts
    --            goalExport ldpth ("clustering-cat-" ++ show ky) $ xys
    --    runGnuplotWithVariables ldpth "projection" [("nm",nm)]
--
--    let mvn :: Source # MultivariateNormal 5
--        mvn = mle smps
--        npca = last emnpcas
--        spca :: StandardPCA 5 2
--        spca = last emspcas
--        npcamvn = pcaObservableDistribution npca
--        spcamvn = standardPCAToMultivariateNormal spca
--        crrs = getCorrelations mvn
--        npcacrrs = getCorrelations $ transition npcamvn
--        spcacrrs = getCorrelations spcamvn
--
--    putStrLn "Uniquenesses:"
--    print . S.toList $ pcaUniqueness npca
--
--    goalExport ldpth "correlations" $ zip3 crrs npcacrrs spcacrrs
--    runGnuplot ldpth "correlations-scatter"

--    let npca1 = last emnpcas
--        prr1,prr2 :: Source # MultivariateNormal 2
--        prr1 = fromTuple (0,0,1,0,1)
--        prr2 = fromTuple (1,0,2,1,4)
--        igh1 = joinConjugatedHarmonium npca1 $ transition prr1
--        igh2 = joinConjugatedHarmonium npca1 $ transition prr2
--    print igh1
--    print. toNatural $ toMean igh1
--    print igh2
--    print . toNatural $ toMean igh2
--    print $ euclideanDistance (toNatural $ toMean igh2) igh2
--    print $ (toNatural $ toMean igh2) - igh2
--
--    print $ toMean igh2
--    smps' <- realize $ sample 1000000 igh2
--    let igh2' :: Mean # IsotropicGaussianHarmonium 4 2
--        igh2' = averageSufficientStatistic smps'
--    print igh2'
--    print $ euclideanDistance (toMean igh2) igh2'



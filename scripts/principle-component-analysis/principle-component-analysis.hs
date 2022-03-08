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

type HMOG n m k = AffineHarmonium Tensor (MVNMean n) (MVNMean m) (IsotropicNormal n) (Mixture (MultivariateNormal m) k)
type HMOG2 n m k = AffineMixture (MultivariateNormal m) (IsotropicGaussianHarmonium n m) k

hmog1to2
    :: ( KnownNat n, KnownNat m, KnownNat k )
    => c # HMOG n m k
    -> c # HMOG2 n m k
hmog1to2 hmog =
    let (lkl,mog) = split hmog
        (aff,cats) = split mog
        (mvn,tns) = split aff
     in join (join (lkl `join` mvn) tns) cats

hmog2to1
    :: ( KnownNat n, KnownNat m, KnownNat k )
    => c # HMOG2 n m k
    -> c # HMOG n m k
hmog2to1 hmog =
    let (bigaff,cats) = split hmog
        (iso,tns) = split bigaff
        (lkl,mvn) = split iso
     in join lkl (join (join mvn tns) cats)

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
    -> Source # MultivariateNormal n
standardPCAToMultivariateNormal pca =
    let (mus,vr,mtx) = splitStandardPCA pca
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = realToFrac vr * S.matrixIdentity
     in joinMultivariateNormal mus $ mtx1 + mtx2

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

multivariateNormalLogLikelihood
    :: KnownNat n => Source # MultivariateNormal n -> S.Vector n Double -> Double
multivariateNormalLogLikelihood p xs =
    let (mus,sgma) = splitMultivariateNormal p
        nrm = (* (-0.5)) . log . S.determinant $ 2*pi*sgma
        dff = xs - mus
        expval = S.dotProduct dff $ S.matrixVectorMultiply (S.pseudoInverse sgma) dff
     in nrm - expval / 2

standardPCAPosteriors
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> StandardPCA n k
    -> [Source # MultivariateNormal k]
standardPCAPosteriors xs pca =
    let (mu,vr,wmtx) = splitStandardPCA pca
        wmtxtr = S.transpose wmtx
        vrmtx = realToFrac vr * S.matrixIdentity
        wmlts = S.matrixMatrixMultiply wmtxtr wmtx
        mmtx = wmlts + vrmtx
        mmtxinv = S.inverse mmtx
        rsds = [ x - mu | x <- xs ]
        muzs = S.matrixVectorMultiply (S.matrixMatrixMultiply mmtxinv wmtxtr) <$> rsds
        cvrz = realToFrac (recip vr) * mmtx
     in (`joinMultivariateNormal` cvrz) <$> muzs

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
        imtx = fromMatrix $ realToFrac invsg * cwmtx
     in join (join thtmu thtsg) imtx

vr0 :: Double
vr0 = 2
prr0,prr1,prr2 :: Source # MultivariateNormal 2
prr0 = fromTuple (0,0,vr0,0,vr0)
prr1 = fromTuple (0.01,0,vr0,0,vr0)
prr2 = fromTuple (0,0.01,vr0,0,vr0)

nprr0,nprr1,nprr2 :: Natural # MultivariateNormal 2
nprr0 = transition prr0
nprr1 = transition prr1
nprr2 = transition prr2

mog0 :: Natural # Mixture (MultivariateNormal 2) 2
mog0 = joinNaturalMixture (S.fromTuple (nprr0,nprr1,nprr2)) 0

mog1 :: Natural # Mixture (MultivariateNormal 2) 0
mog1 = breakPoint nprr0

eps :: Double
eps = 3e-3

nstps :: Int
nstps = 1000

ighEM
    :: [S.Vector 4 Double]
    -> Natural # IsotropicGaussianHarmonium 4 2
    -> Natural # IsotropicGaussianHarmonium 4 2
ighEM zs igh0 =
    let igh1 = expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps
        pca = fst $ split igh1
     in joinConjugatedHarmonium pca nprr0

hmogEM
    :: [S.Vector 4 Double]
    -> Natural # HMOG 4 2 2
    -> Natural # HMOG 4 2 2
hmogEM zs igh0 =
    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps

--hmogConj
--    :: Natural # HMOG2 4 2 2 -> (Double,Natural # Categorical 2)
--hmogConj hrm = conjugationParameters (fst $ split hrm)

expectationMaximizationAscent'
    :: (KnownNat n, KnownNat m, KnownNat k)
    => Double
    -> GradientPursuit
    -> Sample (IsotropicNormal n)
    -> Natural # HMOG n m k
    -> [Natural # HMOG n m k]
expectationMaximizationAscent' eps gp zs nhrm =
    let mhrm' = expectationStep zs nhrm
     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm


--- Instances ---


instance (KnownNat n, KnownNat m, KnownNat k)
  => Transition Natural Mean (HMOG n m k) where
      transition = hmog2to1 . transition . hmog1to2

instance (KnownNat n, KnownNat m, KnownNat k)
  => Generative Natural (HMOG n m k) where
      sample n hmog = do
          let (pca,mog) = splitConjugatedHarmonium hmog
          yzs <- sample n mog
          xs <- mapM samplePoint $ pca >$>* (fst <$> yzs)
          return $ zip xs yzs



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let (smps,cats) = unzip $ parseCSV csvstr
    --print smps

    let nvx :: Natural # IsotropicNormal 4
        nvx = toNatural $ averageSufficientStatistic smps
        svx = toSource nvx

        (mus,vr) = split svx

    (lds0 :: Cartesian # Tensor (Replicated 4 NormalMean) (Replicated 2 NormalMean))
        <- realize $ uniformInitialize (-1,1)

    let spca0 :: StandardPCA 4 2
        spca0 = joinStandardPCA (coordinates mus) (S.head $ coordinates vr) (toMatrix lds0)
        npca0 :: Natural # PrincipleComponentAnalysis 4 2
        npca0 = standardToNaturalPCA spca0
        igh0 = joinConjugatedHarmonium npca0 $ transition prr0
        hmog0 :: Natural # HMOG 4 2 2
        hmog0 = joinConjugatedHarmonium npca0 mog0
    --    hmog1 :: Natural # HMOG 4 2 0
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
    --let hmog0' :: Mean # HMOG 4 2 2
    --    hmog0' = averageSufficientStatistic smps''
    --print hmog0'
    --print $ euclideanDistance (toMean igh0) igh0'
    --print $ euclideanDistance (toMean hmog0) hmog0'




    let emnpcas = take nepchs $ iterate (pcaExpectationMaximization smps) npca0
    --let emspcas = take nepchs
    --        $ iterate (standardPCAExpectationMaximization smps) spca0
    --let emighs = take nepchs $ iterate (ighEM smps) igh0
    --let emhmogs = take nepchs $ iterate (hmogEM smps) hmog0

    --let lls = do
    --        sz <- standardPCAToMultivariateNormal <$> emspcas
    --        return . average $ multivariateNormalLogLikelihood sz <$> smps


    let nlls = logLikelihood smps . pcaObservableDistribution <$> emnpcas

    putStrLn "Natural PCA LL Ascent:"
    mapM_ print nlls

    let npca1 = last emnpcas
        prjctn1 :: Natural # Affine Tensor (MVNMean 2) (MultivariateNormal 2) (MVNMean 4)
        prjctn1 = fst . split . transposeHarmonium . joinConjugatedHarmonium npca1 . transition
            $ joinMultivariateNormal 0 S.matrixIdentity
        prjcts1 = fst . splitMultivariateNormal . toSource <$> prjctn1 >$>* smps

    let mogs = take nepchs $ iterate (expectationMaximization prjcts1) mog0

    let moglls = logLikelihood prjcts1 <$> mogs
        hmoglls = logLikelihood smps . joinConjugatedHarmonium npca1 <$> mogs
    putStrLn "Mog vs HMog LL Ascent:"
    mapM_ print $ zip moglls hmoglls

    let mog1 = last mogs
        hmog1 = joinConjugatedHarmonium npca1 mog1


    let emhmogs = take nepchs . iterate (hmogEM smps) $ joinConjugatedHarmonium npca1 mog0
        hmoglls' = logLikelihood smps <$> emhmogs

    putStrLn "True HMog LL Ascent:"
    mapM_ print hmoglls'

    let hmog2 = last emhmogs
        (npca2,mog2) = splitConjugatedHarmonium hmog2

        prjctn2 :: Natural # Affine Tensor (MVNMean 2) (MultivariateNormal 2) (MVNMean 4)
        prjctn2 = fst . split . transposeHarmonium . joinConjugatedHarmonium npca2 . transition
            $ joinMultivariateNormal 0 S.matrixIdentity
        prjcts2 = fst . splitMultivariateNormal . toSource <$> prjctn2 >$>* smps


    --let lls = do
    --        (nz,sz,iz,hz) <- L.zip4 (pcaObservableDistribution <$> emnpcas)
    --                       (standardPCAToMultivariateNormal <$> emspcas)
    --                       emighs
    --                       emhmogs
    --        return ( logLikelihood smps nz
    --               , average $ multivariateNormalLogLikelihood sz <$> smps
    --               , logLikelihood smps iz
    --               , logLikelihood smps hz )

    --putStrLn "LL Ascent:"
    --mapM_ print lls

    let catmp1,catmp2 :: M.Map Int [S.Vector 2 Double]
        catmp1 = M.fromListWith (++) $ zip cats $ (:[]) <$> prjcts1
        catmp2 = M.fromListWith (++) $ zip cats $ (:[]) <$> prjcts2
        kys = M.keys catmp1
        clstrmp1 = fst $ split mog1
        clstrmp2 = fst $ split mog2

    --let prjcts = S.toList . fst . splitMultivariateNormal <$> standardPCAPosteriors smps (last emspcas)
    ----    catprjcts = L.transpose $ cats : L.transpose prjcts
    sequence_ $ do
        ky <- kys
        let pnts1 = S.toList <$> catmp1 M.! ky
            pnts2 = S.toList <$> catmp2 M.! ky
            mvn1 = toSource $ clstrmp1 >.>* ky
            xys1 = bivariateNormalConfidenceEllipse 100 1 mvn1
            mvn2 = toSource $ clstrmp2 >.>* ky
            xys2 = bivariateNormalConfidenceEllipse 100 1 mvn2
        return $ do
            print $ splitMultivariateNormal mvn1
            print $ splitMultivariateNormal mvn2
            goalExport ldpth ("projection1-cat-" ++ show ky) $ pnts1
            goalExport ldpth ("clustering1-cat-" ++ show ky) $ xys1
            goalExport ldpth ("projection2-cat-" ++ show ky) $ pnts2
            goalExport ldpth ("clustering2-cat-" ++ show ky) $ xys2
    runGnuplot ldpth "projection"
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



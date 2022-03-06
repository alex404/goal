#! stack runghc

{-# LANGUAGE ScopedTypeVariables,TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds,KindSignatures,TypeFamilies,NoStarIsType,UndecidableInstances
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



--- Instances ---



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let (smps,cats) = unzip $ parseCSV csvstr
    print smps

    let mvx :: Source # IsotropicNormal 4
        mvx = toSource . toNatural $ averageSufficientStatistic smps

        (mus,vr) = split mvx

    (lds0 :: Cartesian # Tensor (Replicated 4 NormalMean) (Replicated 2 NormalMean))
        <- realize $ uniformInitialize (-1,1)

    let sfa0 :: StandardPCA 4 2
        sfa0 = joinStandardPCA (coordinates mus) (S.head $ coordinates vr) (toMatrix lds0)
--        nfa0 = standardToNaturalPCA sfa0
--
--    let emnfas = take nepchs $ iterate (pcaExpectationMaximization' smps) nfa0
    let emsfas = take nepchs
            $ iterate (standardPCAExpectationMaximization smps) sfa0

    let lls = do
            sz <- standardPCAToMultivariateNormal <$> emsfas
            return . average $ multivariateNormalLogLikelihood sz <$> smps


--    let lls = do
--            (nz,sz) <- zip (pcaObservableDistribution <$> emnfas)
--                           (standardPCAToMultivariateNormal <$> emsfas)
--            return ( logLikelihood smps nz
--                   , average $ multivariateNormalLogLikelihood sz <$> smps)
--
    putStrLn "LL Ascent:"
    mapM_ print lls


    let prjcts = S.toList . fst . splitMultivariateNormal <$> standardPCAPosteriors smps (last emsfas)
    --    catprjcts = L.transpose $ cats : L.transpose prjcts
        catmp :: M.Map Int ([[Double]])
        catmp = M.fromListWith (++) $ zip cats $ (:[]) <$> prjcts
        kys = M.keys catmp

    sequence_ $ do
        ky <- kys
        let pnts = catmp M.! ky
        return . goalExport ldpth ("projection-cat-" ++ show ky) $ pnts
    runGnuplot ldpth "projection"
--
--    let mvn :: Source # MultivariateNormal 5
--        mvn = mle smps
--        nfa = last emnfas
--        sfa :: StandardPCA 5 2
--        sfa = last emsfas
--        nfamvn = pcaObservableDistribution nfa
--        sfamvn = standardPCAToMultivariateNormal sfa
--        crrs = getCorrelations mvn
--        nfacrrs = getCorrelations $ transition nfamvn
--        sfacrrs = getCorrelations sfamvn
--
--    putStrLn "Uniquenesses:"
--    print . S.toList $ pcaUniqueness nfa
--
--    goalExport ldpth "correlations" $ zip3 crrs nfacrrs sfacrrs
--    runGnuplot ldpth "correlations-scatter"

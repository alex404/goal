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



--- Globals ---


-- Data file pulled from:
-- https://userpage.fu-berlin.de/soga/300/30100_data_sets/food-texture.csv
ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/food-texture.dat"

-- Training --

nepchs :: Int
nepchs = 1

-- Functions --

parseCSV :: String -> Sample (MultivariateNormal 5)
parseCSV csvstr = do
    csv <- tail $ lines csvstr
    let lst = read $ '[' : drop 5 csv ++ "]"
    return . fromJust $ S.fromList lst

getCorrelations
    :: Source # MultivariateNormal 5
    -> [Double]
getCorrelations mvn = do
    let crrs = multivariateNormalCorrelations mvn
    (idx,crw) <- zip [1 :: Int ..] . S.toList $ S.toRows crrs
    drop idx $ S.toList crw

-- Source Factor Analysis --

data SourceFactorAnalysis (n :: Nat) (k :: Nat)

joinSourceFactorAnalysis0
    :: (KnownNat n, KnownNat k)
    => Source # Replicated n Normal -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> Source # SourceFactorAnalysis n k
joinSourceFactorAnalysis0 muvrs mtx =
    let (mus,vrs) = S.toPair . S.toColumns . S.fromRows
            . S.map coordinates $ splitReplicated muvrs
     in Point $ mus S.++ vrs S.++ G.toVector mtx

joinSourceFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => S.Vector n Double -- ^ Mean bias
    -> S.Vector n Double -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> Source # SourceFactorAnalysis n k
joinSourceFactorAnalysis mus vrs mtx =
    Point $ mus S.++ vrs S.++ G.toVector mtx

splitSourceFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => Source # SourceFactorAnalysis n k
    -> (S.Vector n Double, S.Vector n Double, S.Matrix n k Double)
splitSourceFactorAnalysis (Point cs) =
    let (mus,cs') = S.splitAt cs
        (vrs,mtx) = S.splitAt cs'
     in (mus,vrs,G.Matrix mtx)

sourceFAToMultivariateNormal
    :: (KnownNat n, KnownNat k)
    => Source # SourceFactorAnalysis n k
    -> Source # MultivariateNormal n
sourceFAToMultivariateNormal fan =
    let (mus,vrs,mtx) = splitSourceFactorAnalysis fan
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = S.diagonalMatrix vrs
     in joinMultivariateNormal mus $ mtx1 + mtx2

sourceFactorAnalysisExpectationMaximization
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Source # SourceFactorAnalysis n k
    -> Source # SourceFactorAnalysis n k
sourceFactorAnalysisExpectationMaximization xs fan =
    let (mu,vrs,wmtx) = splitSourceFactorAnalysis fan
        wmtxtr = S.transpose wmtx
        vrinv = S.pseudoInverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtxtr vrinv) wmtx
        gmtx = S.pseudoInverse $ S.matrixIdentity + mlts
        xht = average xs
        rsds = [ x - mu | x <- xs ]
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtxtr) vrinv
        muhts = S.matrixVectorMultiply mlts' <$> rsds
        invsgm = S.pseudoInverse . (gmtx +) . S.averageOuterProduct $ zip muhts muhts
        wmtx0 = S.averageOuterProduct $ zip rsds muhts
        wmtx' = S.matrixMatrixMultiply wmtx0 invsgm
        vrs0 = S.matrixMatrixMultiply wmtx (S.transpose wmtx0)
        smtx = S.averageOuterProduct $ zip rsds rsds
        vrs' = S.takeDiagonal $ smtx - vrs0
     in joinSourceFactorAnalysis xht vrs' wmtx'

sourceFactorAnalysisExpectationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Source # SourceFactorAnalysis n k
    -> Mean # MultivariateNormal k
sourceFactorAnalysisExpectationStep xs fan =
    let (mu,vrs,wmtx) = splitSourceFactorAnalysis fan
        wmtxtr = S.transpose wmtx
        vrinv = S.pseudoInverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtxtr vrinv) wmtx
        gmtx = S.pseudoInverse $ S.matrixIdentity + mlts
        --xht = average xs
        rsds = [ x - mu | x <- xs ]
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtxtr) vrinv
        muhts = S.matrixVectorMultiply mlts' <$> rsds
        sghts = (gmtx +) <$> zipWith S.outerProduct muhts muhts
     in joinMeanMultivariateNormal (average muhts) (average sghts)

sourceFactorAnalysisMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => Mean # LinearGaussianHarmonium n k
    -> Source # SourceFactorAnalysis n k
sourceFactorAnalysisMaximizationStep hrm =
    let (nz,nzx,nx) = splitHarmonium hrm
        (muz,etaz) = splitMeanMultivariateNormal nz
        (mux,etax) = splitMeanMultivariateNormal nx
        outrs = toMatrix nzx - S.outerProduct muz mux
        wmtx = S.matrixMatrixMultiply outrs $ S.inverse etax
        zcvr = etaz - S.outerProduct muz muz
        vrs = S.takeDiagonal $ zcvr - S.matrixMatrixMultiply wmtx (S.transpose outrs)
     in joinSourceFactorAnalysis muz vrs wmtx

multivariateNormalLogLikelihood
    :: KnownNat n => Source # MultivariateNormal n -> S.Vector n Double -> Double
multivariateNormalLogLikelihood p xs =
    let (mus,sgma) = splitMultivariateNormal p
        nrm = (* (-0.5)) . log . S.determinant $ 2*pi*sgma
        dff = xs - mus
        expval = S.dotProduct dff $ S.matrixVectorMultiply (S.pseudoInverse sgma) dff
     in nrm - expval / 2

-- Tests

sourceToNaturalFA
    :: (KnownNat n, KnownNat k)
    => Source # SourceFactorAnalysis n k
    -> Natural # FactorAnalysis n k
sourceToNaturalFA sfa =
    let (cmu,cvr,cwmtx) = splitSourceFactorAnalysis sfa
        invsg = recip cvr
        thtmu = invsg * cmu
        thtsg = -0.5 * invsg
        imtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) cwmtx
        nrms = joinReplicated $ S.zipWith (curry fromTuple) thtmu thtsg
     in join nrms $ fromMatrix imtx

--principalComponentAnalysis :: KnownNat k => Source # MultivariateNormal k -> [(Double,Double)]
--principalComponentAnalysis mnrm =
--    let cvr = snd $ splitMultivariateNormal mnrm
--        eigs = reverse . L.sort . map realPart . S.toList . fst $ S.eigens cvr
--        expvr = scanl1 (+) eigs
--     in zip eigs $ (/last expvr) <$> expvr



--- Instances ---


instance (KnownNat n, KnownNat k) => Manifold (SourceFactorAnalysis n k) where
    type Dimension (SourceFactorAnalysis n k) = 2*n + k*n

instance (KnownNat n, KnownNat k) => Statistical (SourceFactorAnalysis n k) where
    type SamplePoint (SourceFactorAnalysis n k) = (S.Vector n Double, S.Vector k Double)



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let smps = parseCSV csvstr

    let mvx :: Mean # Replicated 5 Normal
        mvx = averageSufficientStatistic smps

    (lds0 :: Cartesian # Tensor (Replicated 5 NormalMean) (Replicated 2 NormalMean)) <- realize $ uniformInitialize (-1,1)

    let sfa0 :: Source # SourceFactorAnalysis 5 2
        sfa0 = joinSourceFactorAnalysis0 (transition mvx) (toMatrix lds0)
        fa0 = sourceToNaturalFA sfa0

    let emfas = take nepchs $ iterate (factorAnalysisExpectationMaximization smps) fa0
    let emsfas = take nepchs
            $ iterate (sourceFactorAnalysisExpectationMaximization smps) sfa0


    let lls = do
            (nz,sz) <- zip (factorAnalysisObservableDistribution <$> emfas)
                           (sourceFAToMultivariateNormal <$> emsfas)
            return ( logLikelihood smps nz
                   , average $ multivariateNormalLogLikelihood sz <$> smps)

    putStrLn "LL Ascent:"
    --mapM_ print lls


    let mvn :: Source # MultivariateNormal 5
        mvn = mle smps
        fa = last emfas
        sfa = last emsfas
        famvn = factorAnalysisObservableDistribution fa
        sfamvn = sourceFAToMultivariateNormal sfa
        crrs = getCorrelations mvn
        facrrs = getCorrelations $ transition famvn
        sfacrrs = getCorrelations sfamvn

    goalExport ldpth "correlations" $ zip3 crrs facrrs sfacrrs
    runGnuplot ldpth "correlations-scatter"

    --let (a,b,c) = splitHarmonium . toMean . toNatural $ expectationStep smps lgh
    --print . L.sort $ (listCoordinates c)
    --print . L.sort . listCoordinates $ sourceFactorAnalysisExpectationStep smps sfa0


    --let (a,b,c) = splitHarmonium . toMean . toNatural $ expectationStep smps lgh
    --print . L.sort . listCoordinates . factorAnalysisFromLinearGaussianHarmonium $ expectationMaximization smps lgh
    --print . L.sort . listCoordinates . sourceToNaturalFA $ sourceFactorAnalysisExpectationMaximization smps sfa0


    let lgh0 = naturalFactorAnalysisToLGH fa0
        smpmu = S.toList $ average smps

    let fa1 = factorAnalysisExpectationMaximization smps fa0
        sfa1 = sourceFactorAnalysisExpectationMaximization smps sfa0
        sfa1' = sourceFactorAnalysisMaximizationStep $ expectationStep smps lgh0

    --print . L.sort . listCoordinates $ fa0 - sourceToNaturalFA sfa0
    --print . L.sort . listCoordinates $ fa1 - sourceToNaturalFA sfa1
    print . listCoordinates $ sfa1 - sfa1'

naturalFactorAnalysisToLGH
    :: (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> Natural # LinearGaussianHarmonium n k
naturalFactorAnalysisToLGH fa =
    let ltnt = toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1
        (nzs,tns) = split fa
        (mus,vrs) = S.toPair . S.toColumns . S.fromRows
            . S.map coordinates $ splitReplicated nzs
        cvr = S.diagonalMatrix vrs
        mvn = joinNaturalMultivariateNormal mus cvr
        fa' = join mvn tns
     in joinConjugatedHarmonium fa' ltnt



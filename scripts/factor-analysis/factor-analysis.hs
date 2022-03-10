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


-- Reproduced tutorial here: https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/factor-analysis/A-simple-example-of-FA/index.html
-- Data file pulled from:
-- https://userpage.fu-berlin.de/soga/300/30100_data_sets/food-texture.csv
ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/food-texture.dat"

-- Training --

nepchs :: Int
nepchs = 100

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

-- | This is an older, selfcontained implementation of Factor Analysis, against
-- which I can benchmark the EF-based version.

type StandardFactorAnalysis n k = S.Vector (2*n + n*k) Double

joinStandardFactorAnalysis0
    :: (KnownNat n, KnownNat k)
    => Source # Replicated n Normal -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> StandardFactorAnalysis n k
joinStandardFactorAnalysis0 muvrs mtx =
    let (mus,vrs) = S.toPair . S.toColumns . S.fromRows
            . S.map coordinates $ splitReplicated muvrs
     in mus S.++ vrs S.++ G.toVector mtx


joinStandardFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => S.Vector n Double -- ^ Mean bias
    -> S.Vector n Double -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> StandardFactorAnalysis n k
joinStandardFactorAnalysis mus vrs mtx =
    mus S.++ vrs S.++ G.toVector mtx

splitStandardFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => StandardFactorAnalysis n k
    -> (S.Vector n Double, S.Vector n Double, S.Matrix n k Double)
splitStandardFactorAnalysis cs =
    let (mus,cs') = S.splitAt cs
        (vrs,mtx) = S.splitAt cs'
     in (mus,vrs,G.Matrix mtx)

standardFAToMultivariateNormal
    :: (KnownNat n, KnownNat k)
    => StandardFactorAnalysis n k
    -> Source # MultivariateNormal n
standardFAToMultivariateNormal fan =
    let (mus,vrs,mtx) = splitStandardFactorAnalysis fan
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = S.diagonalMatrix vrs
     in joinMultivariateNormal mus $ mtx1 + mtx2

standardFAExpectationMaximization
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> StandardFactorAnalysis n k
    -> StandardFactorAnalysis n k
standardFAExpectationMaximization xs fan =
    let (mu,vrs,wmtx) = splitStandardFactorAnalysis fan
        wmtxtr = S.transpose wmtx
        vrinv = S.pseudoInverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtxtr vrinv) wmtx
        gmtx = S.pseudoInverse $ S.matrixIdentity + mlts
        xht = average xs
        rsds = [ x - xht | x <- xs ]
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtxtr) vrinv
        muhts = S.matrixVectorMultiply mlts' <$> rsds
        invsgm = S.pseudoInverse . (gmtx +) . S.averageOuterProduct $ zip muhts muhts
        wmtx0 = S.averageOuterProduct $ zip rsds muhts
        wmtx' = S.matrixMatrixMultiply wmtx0 invsgm
        vrs0 = S.matrixMatrixMultiply wmtx (S.transpose wmtx0)
        smtx = S.averageOuterProduct $ zip rsds rsds
        vrs' = S.takeDiagonal $ smtx - vrs0
     in joinStandardFactorAnalysis xht vrs' wmtx'

multivariateNormalLogLikelihood
    :: KnownNat n => Source # MultivariateNormal n -> S.Vector n Double -> Double
multivariateNormalLogLikelihood p xs =
    let (mus,sgma) = splitMultivariateNormal p
        nrm = (* (-0.5)) . log . S.determinant $ 2*pi*sgma
        dff = xs - mus
        expval = S.dotProduct dff $ S.matrixVectorMultiply (S.pseudoInverse sgma) dff
     in nrm - expval / 2

-- Tests

standardToNaturalFA
    :: (KnownNat n, KnownNat k)
    => StandardFactorAnalysis n k
    -> Natural # FactorAnalysis n k
standardToNaturalFA sfa =
    let (cmu,cvr,cwmtx) = splitStandardFactorAnalysis sfa
        invsg = recip cvr
        thtmu = invsg * cmu
        thtsg = -0.5 * invsg
        imtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) cwmtx
        nrms = join (Point thtmu) (Point thtsg)
     in join nrms $ fromMatrix imtx


--- Instances ---



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth
    let smps = parseCSV csvstr

    let mvx :: Mean # Replicated 5 Normal
        mvx = averageSufficientStatistic smps

    (lds0 :: Cartesian # Tensor (Replicated 5 NormalMean) (Replicated 2 NormalMean))
        <- realize $ uniformInitialize (-1,1)

    let sfa0 :: StandardFactorAnalysis 5 2
        sfa0 = joinStandardFactorAnalysis0 (transition mvx) (toMatrix lds0)
        nfa0 = standardToNaturalFA sfa0

    let emnfas = take nepchs $ iterate (factorAnalysisExpectationMaximization smps) nfa0
    let emsfas = take nepchs
            $ iterate (standardFAExpectationMaximization smps) sfa0


    let lls = do
            (nz,sz) <- zip (factorAnalysisObservableDistribution <$> emnfas)
                           (standardFAToMultivariateNormal <$> emsfas)
            return ( logLikelihood smps nz
                   , average $ multivariateNormalLogLikelihood sz <$> smps)

    putStrLn "LL Ascent:"
    mapM_ print lls

    let mvn :: Source # MultivariateNormal 5
        mvn = mle smps
        nfa = last emnfas
        sfa :: StandardFactorAnalysis 5 2
        sfa = last emsfas
        nfamvn = factorAnalysisObservableDistribution nfa
        sfamvn = standardFAToMultivariateNormal sfa
        crrs = getCorrelations mvn
        nfacrrs = getCorrelations $ transition nfamvn
        sfacrrs = getCorrelations sfamvn

    putStrLn "Uniquenesses:"
    print . S.toList $ factorAnalysisUniqueness nfa

    goalExport ldpth "correlations" $ zip3 crrs nfacrrs sfacrrs
    runGnuplot ldpth "correlations-scatter"

-- Reproduced tutorial here: https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/factor-analysis/A-simple-example-of-FA/index.html
-- Data file pulled from:
-- https://userpage.fu-berlin.de/soga/300/30100_data_sets/food-texture.csv
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Maybe (fromJust)

--- Globals ---

nepchs :: Int
nepchs = 100

type N = 5
type K = 2

-- Functions --

parseCSV :: String -> Sample (FullNormal N)
parseCSV csvstr = do
    csv <- tail $ lines csvstr
    let lst = read $ '[' : drop 5 csv ++ "]"
    return . fromJust $ S.fromList lst

getCorrelations ::
    Source # FullNormal 5 ->
    [Double]
getCorrelations mvn = do
    let crrs = multivariateNormalCorrelations mvn
    (idx, crw) <- zip [1 :: Int ..] . S.toList $ toRows crrs
    drop idx . S.toList $ coordinates crw

{- | This is an older, selfcontained implementation of Factor Analysis, against
which I can benchmark the EF-based version.
-}
type StandardFactorAnalysis n k = S.Vector (2 * n + n * k) Double

joinStandardFactorAnalysis0 ::
    (KnownNat n, KnownNat k) =>
    -- | Variances
    Source # DiagonalNormal n ->
    -- | Interaction Parameters
    S.Matrix n k Double ->
    StandardFactorAnalysis n k
joinStandardFactorAnalysis0 muvrs mtx =
    let (mus, vrs) = split muvrs
     in coordinates mus S.++ coordinates vrs S.++ G.toVector mtx

joinStandardFactorAnalysis ::
    (KnownNat n, KnownNat k) =>
    -- | Mean bias
    S.Vector n Double ->
    -- | Variances
    S.Vector n Double ->
    -- | Interaction Parameters
    S.Matrix n k Double ->
    StandardFactorAnalysis n k
joinStandardFactorAnalysis mus vrs mtx =
    mus S.++ vrs S.++ G.toVector mtx

splitStandardFactorAnalysis ::
    (KnownNat n, KnownNat k) =>
    StandardFactorAnalysis n k ->
    (S.Vector n Double, S.Vector n Double, S.Matrix n k Double)
splitStandardFactorAnalysis cs =
    let (mus, cs') = S.splitAt cs
        (vrs, mtx) = S.splitAt cs'
     in (mus, vrs, G.Matrix mtx)

standardFAToFullNormal ::
    forall n k.
    (KnownNat n, KnownNat k) =>
    StandardFactorAnalysis n k ->
    Source # FullNormal n
standardFAToFullNormal fan =
    let mus :: S.Vector n Double
        vrs :: S.Vector n Double
        mtx :: S.Matrix n k Double
        (mus, vrs, mtx) = splitStandardFactorAnalysis fan
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = S.diagonalMatrix vrs
     in join (Point mus) . fromTensor . Point . G.toVector $ mtx1 + mtx2

standardFAExpectationMaximization ::
    forall n k.
    (KnownNat n, KnownNat k) =>
    [S.Vector n Double] ->
    StandardFactorAnalysis n k ->
    StandardFactorAnalysis n k
standardFAExpectationMaximization xs fan =
    let (_, vrs, wmtx) = splitStandardFactorAnalysis fan
        wmtxtr = S.transpose wmtx
        vrinv = S.pseudoInverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtxtr vrinv) wmtx
        gmtx = S.pseudoInverse $ S.matrixIdentity + mlts
        xht = average xs
        rsds = [x - xht | x <- xs]
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtxtr) vrinv
        muhts = S.matrixVectorMultiply mlts' <$> rsds
        invsgm = S.pseudoInverse . (gmtx +) . S.averageOuterProduct $ zip muhts muhts
        wmtx0 = S.averageOuterProduct $ zip rsds muhts
        wmtx' = S.matrixMatrixMultiply wmtx0 invsgm
        vrs0 = S.matrixMatrixMultiply wmtx (S.transpose wmtx0)
        smtx = S.averageOuterProduct $ zip rsds rsds
        vrs' = S.takeDiagonal $ smtx - vrs0
     in joinStandardFactorAnalysis xht vrs' wmtx'

-- multivariateNormalLogLikelihood
--    :: KnownNat n => Source # FullNormal n -> S.Vector n Double -> Double
-- multivariateNormalLogLikelihood p xs =
--    let (mus,sgma) = splitMultivariateNormal p
--        nrm = (* (-0.5)) . log . S.determinant $ 2*pi*sgma
--        dff = xs - mus
--        expval = S.dotProduct dff $ S.matrixVectorMultiply (S.pseudoInverse sgma) dff
--     in nrm - expval / 2

-- Tests

standardToNaturalFA ::
    (KnownNat n, KnownNat k) =>
    StandardFactorAnalysis n k ->
    Natural # FactorAnalysis n k
standardToNaturalFA sfa =
    let (cmu, cvr, cwmtx) = splitStandardFactorAnalysis sfa
        invsg = recip cvr
        thtmu = invsg * cmu
        thtsg = -0.5 * invsg
        imtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) cwmtx
        nrms = join (Point thtmu) (Point thtsg)
     in join nrms . Point $ G.toVector imtx

--- Instances ---

--- Main ---

main :: IO ()
main = do
    csvpth <- dataFilePath "food-texture.dat"
    csvstr <- readFile csvpth
    let smps = parseCSV csvstr

    let mvx :: Mean # DiagonalNormal N
        mvx = averageSufficientStatistic smps

    (lds0 :: Cartesian # Tensor (StandardNormal N) (StandardNormal K)) <-
        realize $ uniformInitialize (-1, 1)

    let sfa0 :: StandardFactorAnalysis N K
        sfa0 = joinStandardFactorAnalysis0 (transition mvx) (L.toMatrix $ useLinear lds0)
        nfa0 = standardToNaturalFA sfa0

    let emnfas = take nepchs $ iterate (linearModelExpectationMaximization smps) nfa0
    let emsfas :: [StandardFactorAnalysis 5 2]
        emsfas =
            take nepchs
                $ iterate (standardFAExpectationMaximization smps) sfa0

        ems :: [(Natural # FullNormal 5, Natural # FullNormal 5)]
        ems =
            zip
                (linearModelObservableDistribution <$> emnfas)
                (toNatural . standardFAToFullNormal <$> emsfas)

    let (nlls, slls) = unzip $ do
            (nz, sz) <- ems
            return (logLikelihood smps nz, logLikelihood smps sz)

    -- putStrLn "LL Ascent:"
    -- mapM_ print lls

    let mvn :: Source # FullNormal 5
        mvn = mle smps
        nfa = last emnfas
        sfa :: StandardFactorAnalysis 5 2
        sfa = last emsfas
        nfamvn = linearModelObservableDistribution nfa
        sfamvn = standardFAToFullNormal sfa
        crrs = getCorrelations mvn
        nfacrrs = getCorrelations $ transition nfamvn
        sfacrrs = getCorrelations sfamvn

    putStrLn "Uniquenesses:"
    print . S.toList $ factorAnalysisUniqueness nfa

    let json =
            toJSON
                [ "standard-log-likelihood" .= slls
                , "ef-log-likelihood" .= nlls
                , "data-correlations" .= crrs
                , "standard-factor-analysis-correlations" .= sfacrrs
                , "ef-factor-analysis-correlations" .= nfacrrs
                ]

    --- Process data
    rsltsfl <- resultsFilePath "factor-analysis.json"
    exportJSON rsltsfl json

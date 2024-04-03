{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S

--- Globals ---

--- Factor Analysis

type N = 10
type K = 3

ldngs :: S.Vector K (Source # StandardNormal N)
ldngs =
    S.fromTuple
        ( fromTuple (0.9, -0.1, 0.1, 0.4, -0.8, -0.5, -0.2, 0.3, 0.2, -0.6)
        , fromTuple (0.1, 0.8, -0.1, -0.4, 0.3, 0.2, 0.5, -0.4, 0.2, 0.7)
        , fromTuple (0.1, 0.1, 0.7, -0.2, 0.8, 0.3, -0.3, 0.3, -0.6, 0.4)
        )

diag :: Source # DiagonalNormal N
diag =
    join (Point . S.generate $ \k -> 5 - (fromIntegral k / 3)) $
        fromTuple (0.1, 0.1, 0.1, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3)

sfa :: Source # FactorAnalysis N K
sfa = join diag $ fromColumns ldngs

nfa :: Natural # FactorAnalysis N K
nfa = toNatural sfa

--- Affine CoM Mixture

type CoMBasedMixture n k = AffineMixture (Replicated n Poisson) (Replicated n CoMPoisson) k

cbldngs0 :: Random (Natural # Tensor (Replicated N Poisson) (Categorical K))
cbldngs0 = uniformInitialize (-0.01, 0.01)

cbdiag0 :: Random (Source # Replicated N CoMPoisson)
cbdiag0 = uniformInitialize (0.9, 1.1)

wghts0 :: Random (Natural # Categorical K)
wghts0 = uniformInitialize (-0.01, 0.01)

initializeCBM :: Random (Natural # CoMBasedMixture N K)
initializeCBM = do
    cbldngs <- cbldngs0
    cbdiag <- cbdiag0
    wghts <- wghts0
    let ncbdiag = toNatural cbdiag
    return $ joinHarmonium ncbdiag cbldngs wghts

--- Training

smpn :: Int
smpn = 1000

eps :: Double
eps = 3e-3

nstps, nepchs :: Int
nstps = 500
nepchs = 50

gp :: GradientPursuit
gp = defaultAdamPursuit

--- Helper Functions

discretize :: (KnownNat n) => S.Vector n Double -> S.Vector n Int
discretize = S.map (round . max 0)

nsmpcrl :: Int
nsmpcrl = 10000

cbmSecondOrderStatistics ::
    forall n k.
    (KnownNat n, KnownNat k) =>
    Natural # CoMBasedMixture n k ->
    Random (Source # FullNormal n)
cbmSecondOrderStatistics cbm = do
    xzsmp <- sample nsmpcrl cbm
    let xsmp0 :: [S.Vector n Int]
        xsmp0 = fst <$> xzsmp
        xsmp = G.convert . G.map realToFrac <$> xsmp0
        mu = Point $ average xsmp
        cvr = S.lowerTriangular . S.averageOuterProduct $ zip xsmp xsmp
        mcvr = Point cvr
        mnrm :: Mean # FullNormal n
        mnrm = join mu mcvr
    return $ toSource mnrm

loggingEMStep ::
    (KnownNat n, KnownNat k) =>
    [S.Vector n Int] ->
    (Int, Natural # CoMBasedMixture n k) ->
    IO (Int, Natural # CoMBasedMixture n k)
loggingEMStep xs (k, cbm) = do
    let gbhrms = expectationMaximizationAscent eps gp xs cbm
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show $ logLikelihood xs cbm
            ]
    return (k + 1, gbhrms !! nstps)

gatherStatistics ::
    (KnownNat n) =>
    Source # FullNormal n ->
    (S.Vector n Double, S.Vector n Double, S.Vector (Triangular n) Double, S.Vector n (S.Vector n Double))
gatherStatistics mvn =
    let (mu, cvr) = split mvn
        cvrs = coordinates cvr
        vrs = S.triangularTakeDiagonal cvrs
        ffs = vrs / coordinates mu
        mus = coordinates mu
        crrs = S.map coordinates . toRows $ multivariateNormalCorrelations mvn
     in (mus, ffs, cvrs, crrs)

--- Main ---

main :: IO ()
main = do
    --- Mixtures
    let fadns = linearModelObservableDistribution nfa
        sfadns = toSource fadns
        (famus, faffs, facvrs, facrrs) = gatherStatistics sfadns

    xs <- realize $ sample smpn fadns
    putStrLn "Percent of data that is not strictly positive:"
    let cnt = fromIntegral . length . filter (any (< 0)) $ S.toList <$> xs
    putStrLn $ show ((100 :: Double) * cnt / fromIntegral smpn) ++ "%"

    let ns = discretize <$> xs
        mmvn1 :: Mean # FullNormal N
        mmvn1 = mle $ S.map fromIntegral <$> ns
        mvn1 = toSource mmvn1
        (smpmus, smpffs, smpcvrs, smpcrrs) = gatherStatistics mvn1

    cbm0 <- realize initializeCBM
    kcbms <- iterateM nepchs (loggingEMStep ns) (0, cbm0)

    let cbms = snd <$> kcbms
        cbm1 = last cbms
    cbmmvn <- realize $ cbmSecondOrderStatistics cbm1
    let (cbmmus, cbmffs, cbmcvrs, cbmcrrs) = gatherStatistics cbmmvn

    --- Export
    let json =
            toJSON
                [ "fa-means" .= famus
                , "cbm-means" .= cbmmus
                , "sample-means" .= smpmus
                , "fa-fano-factors" .= faffs
                , "cbm-fano-factors" .= cbmffs
                , "sample-fano-factors" .= smpffs
                , "fa-covariances" .= facvrs
                , "cbm-covariances" .= cbmcvrs
                , "sample-covariances" .= smpcvrs
                , "fa-correlation-matrix" .= facrrs
                , "cbm-correlation-matrix" .= cbmcrrs
                , "sample-correlation-matrix" .= smpcrrs
                ]

    --- Process data
    flnm <- resultsFilePath "mixture-com-based.json"

    exportJSON flnm json

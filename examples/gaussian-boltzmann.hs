{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MonoLocalBinds #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoStarIsType #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Finite (Finite, natToFinite)
import Data.Proxy (Proxy (..))

--- Globals ---

--- Types

type BN = 6
type OS = 2

--- Initialization

unibnd :: Double
unibnd = 1

initializeLatent ::
    Random (Natural # GaussianBoltzmannHarmonium L.PositiveDefinite OS BN)
initializeLatent = do
    bltz :: Natural # Boltzmann BN <- uniformInitialize (-unibnd, unibnd)
    shfts :: Natural # Tensor (StandardNormal OS) (Replicated BN Bernoulli) <-
        uniformInitialize (-unibnd, unibnd)
    let mvn = standardNormal
        lmdl :: Natural # BoltzmannLinearModel L.PositiveDefinite OS BN
        lmdl = join mvn shfts
    return $ joinConjugatedHarmonium lmdl bltz

--- Data generation

circle2d :: Double -> Double -> SamplePoint (FullNormal 2)
circle2d rds t = S.fromTuple (rds * cos t, rds * sin t)

nssd1, nssd2, nssd3 :: Double
nssd1 = 0.2
nssd2 = 0.1
nssd3 = 0.4

noisyCircle :: Double -> Double -> Int -> Random (Sample (FullNormal 2))
noisyCircle nssd rds nsmps = do
    let nsmdl :: Source # IsotropicNormal 2
        nsmdl = fromTuple (0, 0, nssd ^ 2)
    let intrvl = range 0 (2 * pi) nsmps
    mapM (noisyFunction nsmdl (circle2d rds)) intrvl

--- Training

eps :: Double
eps = 0.005

rto1, rto2, rto3 :: Double
rto1 = 0.05
rto2 = 0.35
rto3 = 0.6

nobs, nobs1, nobs2, nobs3 :: Int
nobs = 1000
nobs1 = round $ rto1 * fromIntegral nobs
nobs2 = round $ rto2 * fromIntegral nobs
nobs3 = round $ rto3 * fromIntegral nobs

rds1, rds2, rds3 :: Double
rds1 = 0.5
rds2 = 4
rds3 = 8

nstps, nepchs :: Int
nstps = 500
nepchs = 1000

gp :: GradientPursuit
gp = defaultAdamPursuit

loggingEMStep ::
    (KnownNat n, KnownNat k, 1 <= k) =>
    [S.Vector n Double] ->
    (Int, Natural # GaussianBoltzmannHarmonium L.PositiveDefinite n k) ->
    IO (Int, Natural # GaussianBoltzmannHarmonium L.PositiveDefinite n k)
loggingEMStep xs (k, gbhrm) = do
    let gbhrms = expectationMaximizationAscent eps gp xs gbhrm
    putStrLn $
        concat
            [ "Iteration "
            , show k
            , " Log-Likelihood: "
            , show $ logLikelihood xs gbhrm
            ]
    return (k + 1, gbhrms !! nstps)

--- Statistics

momentMatrixRows ::
    Natural # Boltzmann BN ->
    ([S.Vector BN Double], [S.Vector BN Double])
momentMatrixRows nbltz =
    let blss = pointSampleSpace nbltz
        blss' = S.map (fromIntegral . fromEnum) <$> blss
        prbs = densities nbltz blss
        mcvr = S.weightedAverageOuterProduct $ zip3 prbs blss' blss'
        mmu = S.takeDiagonal mcvr
        mnrm :: Mean # FullNormal BN
        mnrm = join (Point mmu) (Point $ S.lowerTriangular mcvr)
     in ( S.toList $ S.toRows mcvr
        , map coordinates . S.toList . toRows . multivariateNormalCorrelations $ toSource mnrm
        )

--- Plotting

pltres :: Int
pltres = 100

dnspltmn, dnspltmx :: Double
dnspltmn = -10
dnspltmx = 10

pltrng :: [Double]
pltrng = range dnspltmn dnspltmx pltres

dnspltxs :: Sample (FullNormal 2)
dnspltxs = [S.fromTuple (x1, x2) | x1 <- pltrng, x2 <- pltrng]

--- Main ---

main :: IO ()
main = do
    --- Initialization
    gbltz0 <- realize initializeLatent
    xs1 <- realize $ noisyCircle nssd1 rds1 nobs1
    xs2 <- realize $ noisyCircle nssd2 rds2 nobs2
    xs3 <- realize $ noisyCircle nssd3 rds3 nobs3
    xs <- realize . shuffleList $ xs1 ++ xs2 ++ xs3

    -- pstxs0 <- realize $ take 2 <$> shuffleList xs1
    -- pstxs0' <- realize $ take 2 <$> shuffleList xs2
    let pstxs1 = [head xs1, xs1 !! round (fromIntegral nobs1 / 2)]
    let pstxs2 = [xs2 !! round (fromIntegral nobs2 / 4), xs2 !! round (3 * fromIntegral nobs2 / 4)]
    let pstxs = pstxs1 ++ pstxs2

    --- Training
    kgbltzs <- iterateM nepchs (loggingEMStep xs) (0, gbltz0)

    let gbltzs = snd <$> kgbltzs
        gbltz1 = last gbltzs
        lrndns = observableDensities gbltz1 dnspltxs

    --- Confidence ellipses

    let bn :: Finite BN
        bn = natToFinite (Proxy :: Proxy (BN - 1))

    let sngspks :: [S.Vector BN Bool]
        sngspks = [S.generate (== j) | j <- [0 .. bn]]
        lrndlkl = fst $ splitConjugatedHarmonium gbltz1

    let cnfs = bivariateNormalConfidenceEllipse 1000 1 . toSource <$> lrndlkl >$>* sngspks
        prrmtx = momentMatrixRows . snd $ splitConjugatedHarmonium gbltz1
        (lkl1, prr1) = splitConjugatedHarmonium gbltz1
        prrmtxs = unzip [momentMatrixRows $ conjugatedBayesRule lkl1 prr1 pstx | pstx <- pstxs]

    --- Export
    let jsonData =
            toJSON
                [ "observations" .= xs
                , "plot-range-x" .= pltrng
                , "plot-range-y" .= pltrng
                , "learned-density" .= lrndns
                , "component-confidence-ellipses" .= cnfs
                , "prior-moment-matrix" .= prrmtx
                , "posterior-observations" .= pstxs
                , "posterior-moment-matrices" .= prrmtxs
                ]

    rsltfl <- resultsFilePath "gaussian-boltzmann.json"
    exportJSON rsltfl jsonData

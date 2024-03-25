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
        lmdl = join mvn shfts
    return $ joinConjugatedHarmonium lmdl bltz

--- Data generation

circle2d :: Double -> Double -> SamplePoint (FullNormal 2)
circle2d rds t = S.fromTuple (rds * cos t, rds * sin t)

nsmdl :: Source # FullNormal 2
nsmdl = fromTuple (0, 0, 0.2, 0, 0.2)

noisyCircle :: Double -> Int -> Random (Sample (FullNormal 2))
noisyCircle rds nsmps = do
    let intrvl = range 0 (2 * pi) nsmps
    mapM (noisyFunction nsmdl (circle2d rds)) intrvl

--- Training

eps :: Double
eps = 3e-3

nstps, ndatsmps, nepchs :: Int
nstps = 2000
ndatsmps = 1000
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
dnspltmn = -5
dnspltmx = 5

pltrng :: [Double]
pltrng = range dnspltmn dnspltmx pltres

dnspltxs :: Sample (FullNormal 2)
dnspltxs = [S.fromTuple (x1, x2) | x1 <- pltrng, x2 <- pltrng]

--- Main ---

main :: IO ()
main = do
    --- Initialization
    gbltz0 <- realize initializeLatent
    let frc :: Double
        frc = 0.85
        ndatsmps1 = round $ frc * fromIntegral ndatsmps
        ndatsmps2 = ndatsmps - ndatsmps1

    xs' <- realize $ noisyCircle 3 ndatsmps1
    xs'' <- realize $ noisyCircle 0.1 ndatsmps2
    xs <- realize . shuffleList $ xs' ++ xs''

    -- pstxs0 <- realize $ take 2 <$> shuffleList xs'
    -- pstxs0' <- realize $ take 2 <$> shuffleList xs''
    let pstxs0 = [head xs', xs' !! round (fromIntegral ndatsmps1 / 2)]
    let pstxs0' = [xs'' !! round (fromIntegral ndatsmps2 / 4), xs'' !! round (3 * fromIntegral ndatsmps2 / 4)]
    let pstxs = pstxs0 ++ pstxs0'

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

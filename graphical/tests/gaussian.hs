{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# OPTIONS_GHC -Wno-orphans #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Control.Monad (replicateM)
import Data.Maybe (fromJust)
import Data.Proxy (Proxy (..))

import System.Exit (exitFailure)
import Test.QuickCheck

--- Globals ---

type N = 3
type K = 2

genBounds :: (Double, Double)
genBounds = (-10, 10)

numericalTolerance :: Double
numericalTolerance = 1e-1

--- Helpers Functions

approxEqual :: Double -> Double -> Bool
approxEqual x y = abs (x - y) < numericalTolerance

mvnApproxEqual ::
    c # MultivariateNormal f (N + K) ->
    c # MultivariateNormal f (N + K) ->
    Bool
mvnApproxEqual mvn1 mvn2 =
    and $ zipWith approxEqual (listCoordinates mvn1) (listCoordinates mvn2)

meanGaussianToHarmonium ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Mean # MultivariateNormal f (N + K) ->
    Mean # LinearGaussianHarmonium f N K
meanGaussianToHarmonium mvn =
    let (mu, sgma) = split mvn
        (mux, muz) = S.splitAt $ coordinates mu
        sgmarws = toRows $ toTensor sgma
        lftrws = S.map (S.take . coordinates) sgmarws
        rghtrws = S.map (S.drop . coordinates) sgmarws
        sgmax = fromRows . S.map Point $ S.take lftrws
        sgmaxz = fromRows . S.map Point $ S.take rghtrws
        sgmaz = fromRows . S.map Point $ S.drop rghtrws
        nx = join (Point mux) $ fromTensor sgmax
        nz = join (Point muz) $ fromTensor sgmaz
        lkl = join nx sgmaxz
     in join lkl nz

meanHarmoniumToGaussian ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Mean # LinearGaussianHarmonium f N K ->
    Mean # MultivariateNormal f (N + K)
meanHarmoniumToGaussian lgm =
    let (lkl, nz) = split lgm
        (nx, sgmaxz) = split lkl
        (muz, sgmaz) = split nz
        (mux, sgmax) = split nx
        lftrws = toRows (toTensor sgmax) S.++ toRows (toTensor $ transpose sgmaxz)
        rghtrws = toRows (toTensor sgmaxz) S.++ toRows (toTensor sgmaz)
        sgma = fromRows $ S.zipWith (\lft rght -> Point $ coordinates lft S.++ coordinates rght) lftrws rghtrws
        mu = Point $ coordinates mux S.++ coordinates muz
     in join mu $ fromTensor sgma

naturalGaussianToHarmonium ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Natural # MultivariateNormal f (N + K) ->
    Natural # LinearGaussianHarmonium f N K
naturalGaussianToHarmonium mvn =
    let (mu, sgma) = splitNaturalNormal mvn
        (mux, muz) = S.splitAt $ coordinates mu
        sgmarws = toRows $ toTensor sgma
        lftrws = S.map (S.take . coordinates) sgmarws
        rghtrws = S.map (S.drop . coordinates) sgmarws
        sgmax = fromRows . S.map Point $ S.take lftrws
        sgmaxz = fromRows . S.map Point $ S.take rghtrws
        sgmaz = fromRows . S.map Point $ S.drop rghtrws
        nx = joinNaturalNormal (Point mux) $ fromTensor sgmax
        nz = joinNaturalNormal (Point muz) $ fromTensor sgmaz
        lkl = join nx $ sgmaxz * 2
     in join lkl nz

naturalHarmoniumToGaussian ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Natural # LinearGaussianHarmonium f N K ->
    Natural # MultivariateNormal f (N + K)
naturalHarmoniumToGaussian lgm =
    let (lkl, nz) = split lgm
        (nx, sgmaxz0) = split lkl
        sgmaxz = sgmaxz0 / 2
        (muz, sgmaz) = splitNaturalNormal nz
        (mux, sgmax) = splitNaturalNormal nx
        lftrws = toRows (toTensor sgmax) S.++ toRows (toTensor $ transpose sgmaxz)
        rghtrws = toRows (toTensor sgmaxz) S.++ toRows (toTensor sgmaz)
        sgma = fromRows $ S.zipWith (\lft rght -> Point $ coordinates lft S.++ coordinates rght) lftrws rghtrws
        mu = Point $ coordinates mux S.++ coordinates muz
     in joinNaturalNormal mu $ fromTensor sgma

--- QuickCheck ---

--- Instances

instance (KnownNat n) => Arbitrary (S.Vector n Double) where
    arbitrary = do
        let n = natValInt $ Proxy @n
        xs <- vectorOf n (choose genBounds)
        return . fromJust $ S.fromList xs

instance (KnownNat n) => Arbitrary (L.Linear L.PositiveDefinite n n) where
    arbitrary = do
        let n = natValInt $ Proxy @n
        xs :: [S.Vector n Double] <- replicateM n arbitrary
        let pdmtx = sum [x `S.outerProduct` x | x <- xs]
        return . L.PositiveDefiniteLinear $ S.lowerTriangular pdmtx

instance (KnownNat n) => Arbitrary (Natural # FullNormal n) where
    arbitrary = do
        xs <- arbitrary :: Gen (S.Vector n Double)
        pdmtx <- arbitrary :: Gen (L.Linear L.PositiveDefinite n n)
        let sgma = Point $ L.toVector pdmtx
            mu = Point xs
            smvn :: Source # FullNormal n
            smvn = join mu sgma
        return $ toNatural smvn

--- Properties

gaussianConversion :: Natural # FullNormal (N + K) -> Bool
gaussianConversion mvn =
    let lgh :: Natural # FullGaussianHarmonium N K
        lgh = naturalGaussianToHarmonium mvn
        mvn' :: Natural # FullNormal (N + K)
        mvn' = naturalHarmoniumToGaussian lgh
        mmvn = toMean mvn
        mlgh = meanGaussianToHarmonium mmvn
        mmvn' = meanHarmoniumToGaussian mlgh
     in mvnApproxEqual mvn mvn' && mvnApproxEqual mmvn mmvn'

gaussianDensities :: (S.Vector (N + K) Double, Natural # FullNormal (N + K)) -> Bool
gaussianDensities (xz, mvn) =
    let nvmdns = density mvn xz
        lgh :: Natural # FullGaussianHarmonium N K
        lgh = naturalGaussianToHarmonium mvn
        (x, z) = S.splitAt xz
        lghdns = density lgh (x, z)
     in approxEqual nvmdns lghdns

gaussianNaturalToMean :: Natural # FullNormal (N + K) -> Bool
gaussianNaturalToMean mvn =
    let mmvn = meanHarmoniumToGaussian . toMean $ naturalGaussianToHarmonium mvn
     in mvnApproxEqual (toMean mvn) mmvn

harmoniumDualTransition :: Natural # FullNormal (N + K) -> Bool
harmoniumDualTransition mvn =
    let nlgh = naturalGaussianToHarmonium mvn
     in mvnApproxEqual mvn . naturalHarmoniumToGaussian . toNatural $ toMean nlgh

--- Main ---

main :: IO ()
main = do
    -- Run tests
    putStrLn "Running tests..."
    putStrLn "\nHarmonium-MVN Conversion...\n"
    result1 <- quickCheckResult gaussianConversion
    putStrLn "\nDensity Comparison...\n"
    result2 <- quickCheckResult gaussianDensities
    putStrLn "\nNatural-Mean Transition...\n"
    result3 <- quickCheckResult gaussianNaturalToMean
    putStrLn "\nLGH Transition Duality...\n"
    result4 <- quickCheckResult harmoniumDualTransition
    if all isSuccess [result1, result2, result3, result4]
        then putStrLn "All tests successful!"
        else do
            putStrLn "Some tests failed!"
            exitFailure

-- test :: IO ()
-- test = do
--     let mvn :: Natural # FullNormal (N + K)
--         mvn =
--             Point
--                 . fromJust
--                 $ S.fromList
--                     [ 8.952390101693199e-2
--                     , -0.7135014606442407
--                     , -0.3085838217345925
--                     , 0.44558533724440613
--                     , 0.759498822358389
--                     , -0.37271257758700244
--                     , 0.8210596958380142
--                     , -1.0494211977115449
--                     , -0.19580456080653005
--                     , 0.23675623645297522
--                     , -0.5973628125916197
--                     , 0.5411771550235183
--                     , -1.0570757794957901
--                     , 0.6046950870510807
--                     , -0.7059999860958122
--                     , -0.16932819119877926
--                     , 0.4375449311891449
--                     , 0.21042914224202172
--                     , -0.3245843146251926
--                     , -0.4015447047333673
--                     ]
--         nlgh = naturalGaussianToHarmonium mvn
--         mlgh = toMean nlgh
--         (mfxz, mz) = split mlgh
--         (mx, mxzvr) = split mfxz
--         (mzmu, mzvr) = split mz
--         (mxmu, mxvr) = split mx
--         sgmaxz = mxzvr - mxmu >.< mzmu
--         sgmax = mxvr - mxmu >.< mxmu
--         sgmaz = mzvr - mzmu >.< mzmu
--         sgmaxinv = inverse sgmax
--         prsz = inverse $ toTensor sgmaz - changeOfBasis sgmaxz sgmaxinv
--         prsxz = -dualComposition sgmaxinv sgmaxz prsz
--         nzmu = prsz >.> mzmu + mxmu <.< prsxz
--         nxmu0 = sgmaxinv >.> mxmu
--         nxmu = prsxz >.> mzmu + (nxmu0 - prsxz >.> (nxmu0 <.< sgmaxz))
--         nxvr = -0.5 * extractObservableCovariance prsxz (transpose sgmaxz) sgmaxinv
--         nzvr = -0.5 * prsz
--         nxzvr = -prsxz
--         nx = joinNaturalNormal nxmu nxvr
--         nz = joinNaturalNormal nzmu $ fromTensor nzvr
--         nfxz = join nx nxzvr
--         nlgh' :: Natural # FullGaussianHarmonium N K
--         nlgh' = join nfxz nz
--     -- Origin
--     putStrLn "Origin"
--     let (nlkl0, nz0) = split nlgh
--     let (nx0, nxzvr0) = split nlkl0
--     print nz0
--     print nxzvr0
--     print nx0
--     -- Transition
--     putStrLn "Transition"
--     print nz
--     print nxzvr
--     print nx
--
-- extractObservableCovariance ::
--     (Primal c, KnownCovariance f n, KnownNat k) =>
--     c # Tensor (StandardNormal n) (StandardNormal k) ->
--     Dual c # Tensor (StandardNormal k) (StandardNormal n) ->
--     c # CovarianceMatrix f n ->
--     c # CovarianceMatrix f n
-- extractObservableCovariance sgmaxz nzxvr nxvrinv =
--     let ainv = useLinear nxvrinv
--      in case ainv of
--             L.PositiveDefiniteLinear _ ->
--                 nxvrinv - fromTensor (dualComposition sgmaxz nzxvr nxvrinv)
--             -- L.DiagonalLinear _ ->
--             --     extractDiagonalCovariance sgmaxz nzxvr nxvrinv
--             -- L.ScaleLinear _ ->
--             --     extractScaleCovariance sgmaxz nzxvr nxvrinv
--             _ -> error "extractObservableCovariance: unsupported covariance type"

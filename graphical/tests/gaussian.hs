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

import Control.Concurrent (threadDelay)
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

pointApproxEqual ::
    (Manifold x) =>
    c # x ->
    c # x ->
    Bool
pointApproxEqual mvn1 mvn2 =
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

instance (KnownCovariance f n) => Arbitrary (Source # CovarianceMatrix f n) where
    arbitrary = do
        let n = natValInt $ Proxy @n
        xs :: [S.Vector n Double] <- replicateM n arbitrary
        return $ sum [x >.< x | x <- Point <$> xs]

instance
    (KnownCovariance f n, Transition Source c (MultivariateNormal f n)) =>
    Arbitrary (c # MultivariateNormal f n)
    where
    arbitrary = do
        xs <- arbitrary :: Gen (S.Vector n Double)
        sgma <- arbitrary :: Gen (Source # CovarianceMatrix f n)
        let mu = Point xs
            smvn = join mu sgma
        return $ transition smvn

instance Arbitrary (Natural # DiagonalGaussianHarmonium N K) where
    arbitrary = do
        lgh <- arbitrary :: Gen (Natural # FullGaussianHarmonium N K)
        nx <- arbitrary :: Gen (Natural # DiagonalNormal N)
        let (_, nxz, nz) = splitHarmonium lgh
        return $ joinHarmonium nx nxz nz

instance Arbitrary (Natural # FullGaussianHarmonium N K) where
    arbitrary = do
        nmvn <- arbitrary :: Gen (Natural # FullNormal (N + K))
        return $ naturalGaussianToHarmonium nmvn

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
     in pointApproxEqual mvn mvn' && pointApproxEqual mmvn mmvn'

gaussianDensities :: (S.Vector (N + K) Double, Natural # FullNormal (N + K)) -> Bool
gaussianDensities (xz, mvn) =
    let nvmdns = density mvn xz
        lgh :: Natural # FullGaussianHarmonium N K
        lgh = naturalGaussianToHarmonium mvn
        (x, z) = S.splitAt xz
        lghdns = density lgh (x, z)
     in approxEqual nvmdns lghdns

toMeanComparison :: Natural # FullNormal (N + K) -> Bool
toMeanComparison mvn =
    let mmvn = meanHarmoniumToGaussian . toMean $ naturalGaussianToHarmonium mvn
     in pointApproxEqual (toMean mvn) mmvn

harmoniumDualTransition :: Natural # FullGaussianHarmonium N K -> Bool
harmoniumDualTransition nlgh =
    pointApproxEqual nlgh . toNatural $ toMean nlgh

diagonalHarmoniumDualTransition :: Natural # DiagonalGaussianHarmonium N K -> Bool
diagonalHarmoniumDualTransition nlgh = pointApproxEqual nlgh . toNatural $ toMean nlgh

--- Main ---

isSuccessful :: Result -> IO ()
isSuccessful rslt = putStrLn $ "\nIs Successful: " ++ show (isSuccess rslt) ++ "\n\n"

main :: IO ()
main = do
    -- Run tests
    putStrLn "Running tests..."
    -- putStrLn "\nHarmonium-MVN Conversion...\n"
    -- result1 <- quickCheckResult gaussianConversion
    -- isSuccessful result1
    --
    -- putStrLn "\nDensity Comparison...\n"
    -- result2 <- quickCheckResult gaussianDensities
    -- isSuccessful result2
    --
    -- putStrLn "\nTo-Mean Comparison...\n"
    -- result3 <- quickCheckResult toMeanComparison
    -- isSuccessful result3
    --
    -- putStrLn "\nLGH Transition Duality...\n"
    -- result4 <- quickCheckResult harmoniumDualTransition
    -- isSuccessful result4
    --
    putStrLn "\nLGH Diagonal Transition Duality...\n"
    verboseCheck diagonalHarmoniumDualTransition

--
-- if all isSuccess [result1, result2, result3, result4, result5]
--     then putStrLn "All tests successful!"
--     else do
--         putStrLn "Some tests failed!"
--         exitFailure
--
-- prettyPrintCovariance :: Natural # FullNormal N -> IO ()
-- prettyPrintCovariance = print . useLinear . (negate . (* 2)) . snd . splitNaturalNormal
--
test :: IO ()
test = do
    -- let nlgh :: Natural # DiagonalGaussianHarmonium N K
    --     nlgh =
    --         Point
    --             . fromJust
    --             $ S.fromList
    --                 [ -6.437022872967686e-2
    --                 , -0.13118813252859815
    --                 , -5.017864723174083e-2
    --                 , -6.245326315352055e-3
    --                 , -1.574892070360161e-2
    --                 , -4.257950042932042e-3
    --                 , 1.4938681381595426e-2
    --                 , -7.736837092725053e-3
    --                 , -5.017876978995941e-3
    --                 , 4.160826807797723e-3
    --                 , -2.7449575192035505e-3
    --                 , 4.638875166836262e-3
    --                 , -2.0051343463945777e-2
    --                 , 9.794554238701977e-2
    --                 , -1.4372071236777272e-2
    --                 , 8.825852570395452e-3
    --                 , -7.385179978603333e-3
    --                 ]
    --
    let nlgh :: Natural # DiagonalGaussianHarmonium N K
        nlgh =
            Point
                . fromJust
                $ S.fromList
                    [ -9.886514367608074e-2
                    , 0.21929425582500522
                    , -4.095328590508397e-3
                    , -5.967725112497226e-3
                    , -5.199665945968467e-2
                    , -3.976534788311693e-2
                    , 2.61828854697067e-3
                    , -1.0025529959967491e-2
                    , 4.8289316674747794e-4
                    , -1.149337767786303e-3
                    , 9.281132916373775e-3
                    , -7.631558559740613e-4
                    , 0.2557526715425806
                    , -0.18486490069115694
                    , -1.1511923006576708e-2
                    , 1.9238188473302963e-2
                    , -1.5500083483247438e-2
                    ]

    -- print nlgh
    let diff = subtract nlgh . toNatural $ toMean nlgh
    print diff
    print $ euclideanDistance nlgh $ toNatural $ toMean nlgh

-- harmoniumConversion ::
--     (KnownCovariance g N, KnownCovariance f N) =>
--     Natural # LinearGaussianHarmonium g N K ->
--     Natural # LinearGaussianHarmonium f N K
-- harmoniumConversion lgm =
--     let (lkl, nz) = split lgm
--         (nx, nfxz) = split lkl
--         (nxmu, nxvr) = splitNaturalNormal nx
--         nxvr' = fromTensor $ toTensor nxvr
--         nx' = joinNaturalNormal nxmu nxvr'
--         lkl' = join nx' nfxz
--      in join lkl' nz
--
-- meanHarmoniumConversion ::
--     (KnownCovariance g N, KnownCovariance f N) =>
--     Mean # LinearGaussianHarmonium g N K ->
--     Mean # LinearGaussianHarmonium f N K
-- meanHarmoniumConversion lgm =
--     let (lkl, nz) = split lgm
--         (nx, nfxz) = split lkl
--         (nxmu, nxvr) = split nx
--         nxvr' = fromTensor $ toTensor nxvr
--         nx' = join nxmu nxvr'
--         lkl' = join nx' nfxz
--      in join lkl' nz
--
--
-- print mlgh'

-- print nlgh2

-- print nvrx
--
--
-- let (sx, svrxz, sz) = splitHarmonium slgh
--     (smux, svrx) = split sx
--     (smuz, svrz) = split sz
--     invsvrz = inverse svrz
--     nvrx0 = inverse $ svrx - fromTensor (changeOfBasis (transpose svrxz) invsvrz)
--     nvrx = -0.5 .> nvrx0
--     nxz = dualComposition nvrx0 svrxz invsvrz
--     nmux = nvrx0 >.> smux - nxz >.> smuz
--     nmvn = join nmux nvrx
--     nlkl = join nmvn nxz
--     nz = toNatural sz
--  in joinConjugatedHarmonium (breakChart nlkl) nz
--

-- print $ splitHarmonium nlgh
-- let mlgh = toMean nlgh
-- let (mx, mxz, mz) = splitHarmonium mlgh

-- print $ toNatural mlgh'
-- print $ toNatural mlgh'''

-- print $ splitHarmonium nlgh'
-- let (mx', mxz', mz') = splitHarmonium $ toMean nlgh'

-- print mx
-- print mx'
-- print mxz
-- print mxz'
-- print mz
-- print mz'
--
-- let (mx, msgmaxz, mz) = splitHarmonium mlgh
--     (mux, msgmax) = split mx
--     (muz, msgmaz) = split mz
--     sgmaz = msgmaz - muz >.< muz
--     sgmax = msgmax - mux >.< mux
--     sgmaxz = msgmaxz - mux >.< muz
--     sgmaxinv = inverse sgmax
--     prsz = inverse $ toTensor sgmaz - changeOfBasis sgmaxz sgmaxinv
--     prsxz = negate $ dualComposition sgmaxinv sgmaxz prsz
--     nxvr = -0.5 .> extractObservableCovariance prsxz (transpose sgmaxz) sgmaxinv
--     nmux0 = sgmaxinv >.> mux
--     nmux = (nmux0 - prsxz >.> (nmux0 <.< sgmaxz)) + prsxz >.> muz
--     nx = joinNaturalNormal nmux nxvr
--     nlkl = join nx $ negate prsxz
--     nlgh'' :: Natural # DiagonalGaussianHarmonium N K
--     nlgh'' = joinConjugatedHarmonium nlkl $ toNatural mz
--
-- print nlgh'
-- print nlgh''
--
-- sourceGaussianPrior ::
--     (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
--     Source # MultivariateNormal f (N + K) ->
--     Source # FullNormal K
-- sourceGaussianPrior mvn =
--     let (mu, sgma) = split mvn
--         muz = S.drop $ coordinates mu
--         sgmarws = toRows $ toTensor sgma
--         rghtrws = S.map (S.drop . coordinates) sgmarws
--         sgmaz = fromRows . S.map Point $ S.drop rghtrws
--      in join (Point muz) $ fromTensor sgmaz

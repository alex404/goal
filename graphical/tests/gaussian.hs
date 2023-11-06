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
genBounds = (-1, 1)

numericalTolerance :: Double
numericalTolerance = 1e-12

--- Helpers Functions

approxEqual :: Double -> Double -> Bool
approxEqual x y = abs (x - y) < numericalTolerance

mvnApproxEqual ::
    Natural # MultivariateNormal f (N + K) ->
    Natural # MultivariateNormal f (N + K) ->
    Bool
mvnApproxEqual mvn1 mvn2 =
    and $ zipWith approxEqual (listCoordinates mvn1) (listCoordinates mvn2)

gaussianToHarmonium ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Natural # MultivariateNormal f (N + K) ->
    Natural # LinearGaussianHarmonium f N K
gaussianToHarmonium mvn =
    let (mu, sgma) = splitNaturalNormal mvn
        (mux, muz) = S.splitAt $ coordinates mu
        sgmarws = toRows $ toTensor sgma
        lftrws :: S.Vector (N + K) (S.Vector N Double)
        lftrws = S.map (S.take . coordinates) sgmarws
        rghtrws :: S.Vector (N + K) (S.Vector K Double)
        rghtrws = S.map (S.drop . coordinates) sgmarws
        sgmax = fromRows . S.map Point $ S.take lftrws
        sgmaxz = fromRows . S.map Point $ S.take rghtrws
        sgmaz = fromRows . S.map Point $ S.drop rghtrws
        nx = joinNaturalNormal (Point mux) $ fromTensor sgmax
        nz = joinNaturalNormal (Point muz) $ fromTensor sgmaz
        lkl = join nx $ sgmaxz * 2
     in join lkl nz

harmoniumToGaussian ::
    (KnownCovariance f (N + K), KnownCovariance f K, KnownCovariance f N) =>
    Natural # LinearGaussianHarmonium f N K ->
    Natural # MultivariateNormal f (N + K)
harmoniumToGaussian lgm =
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
    let lgm :: Natural # FullGaussianHarmonium N K
        lgm = gaussianToHarmonium mvn
        mvn' :: Natural # FullNormal (N + K)
        mvn' = harmoniumToGaussian lgm
     in mvnApproxEqual mvn mvn'

gaussianDensities :: (S.Vector (N + K) Double, Natural # FullNormal (N + K)) -> Bool
gaussianDensities (xz, mvn) =
    let nvmdns = density mvn xz
        lgm :: Natural # FullGaussianHarmonium N K
        lgm = gaussianToHarmonium mvn
        (x, z) = S.splitAt xz
        lgmdns = density lgm (x, z)
     in approxEqual nvmdns lgmdns

--- Main ---

test :: IO ()
test = do
    let mvn :: Natural # FullNormal (N + K)
        mvn =
            Point
                . fromJust
                $ S.fromList
                    [ 8.952390101693199e-2
                    , -0.7135014606442407
                    , -0.3085838217345925
                    , 0.44558533724440613
                    , 0.759498822358389
                    , -0.37271257758700244
                    , 0.8210596958380142
                    , -1.0494211977115449
                    , -0.19580456080653005
                    , 0.23675623645297522
                    , -0.5973628125916197
                    , 0.5411771550235183
                    , -1.0570757794957901
                    , 0.6046950870510807
                    , -0.7059999860958122
                    , -0.16932819119877926
                    , 0.4375449311891449
                    , 0.21042914224202172
                    , -0.3245843146251926
                    , -0.4015447047333673
                    ]
        lgh = gaussianToHarmonium mvn
        smvn = toSource mvn
        slgh = toSource lgh
        (smu, ssgma) = split smvn
        (slkl, sz) = split slgh
        (sx, ssgmaxz) = split slkl
    putStrLn "MVN"
    print smu
    print ssgma
    putStrLn "LGH"
    print $ split sx
    print $ split sz
    print ssgmaxz

main :: IO ()
main = do
    -- Run tests
    putStrLn "Running tests..."
    putStrLn "\nDensity Comparison...\n"
    result2 <- verboseCheckResult gaussianDensities
    putStrLn "\nHarmonium-MVN Conversion...\n"
    result1 <- verboseCheckResult gaussianConversion
    if isSuccess result1 && isSuccess result2
        then putStrLn "All tests successful!"
        else do
            putStrLn "Some tests failed!"
            exitFailure

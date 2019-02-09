{-# LANGUAGE
    DataKinds,
    ScopedTypeVariables,
    DeriveGeneric,
    TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

--- Program ---


-- Globals --

nsmps :: Int
nsmps = 100000

mu,kap :: Double
mu = 2
kap = 2

tru :: Source # VonMises
tru = Point $ S.doubleton mu kap

-- Plot

mn,mx :: Double
(mn,mx) = (0,2*pi)

xs :: [Double]
xs = range mn mx 200

nb :: Int
nb = 50

expmnt :: Experiment
expmnt = Experiment "probability" "von-mises"

-- CSV

data Histogram = Histogram
    { radiansBin :: Double
    , histogramDensity :: Double }
    deriving (Generic, Show)

instance FromNamedRecord Histogram
instance ToNamedRecord Histogram
instance DefaultOrdered Histogram

data Density = Density
    { radians :: Double
    , trueDensity :: Double }
    deriving (Generic, Show)

instance FromNamedRecord Density
instance ToNamedRecord Density
instance DefaultOrdered Density





-- Main --

main :: IO ()
main = do

    smps <- realize $ sample nsmps tru
    let cosht = average $ cos <$> smps
        sinht = average $ sin <$> smps

    let (cosht',sinht') = S.toPair . coordinates . dualTransition $ toNatural tru

    putStrLn "Expected Value of Cos (Samples):"
    print cosht
    putStrLn "Expected Value of Cos (Bessel Approx.):"
    print cosht'

    putStrLn "Expected Value of Sin (Samples):"
    print sinht
    putStrLn "Expected Value of Sin (Bessel Approx.):"
    print sinht'

    let (xbns,_,[dnss]) = histograms nb (Just (mn,mx)) [smps]

    goalExportNamed True expmnt Nothing $ zipWith Histogram xbns dnss
    goalExportNamed False expmnt Nothing . zipWith Density xs $ density tru <$> xs

    runGnuplot expmnt Nothing defaultGnuplotOptions "von-mises.gpi"

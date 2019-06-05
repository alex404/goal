{-# LANGUAGE
    DataKinds,
    ScopedTypeVariables,
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Program ---


-- Globals --

alphs :: S.Vector 3 Double
alphs = S.fromTuple (3,7,5)

tru :: Natural # Dirichlet 3
tru = Point alphs

mn,mx :: Double
mn = 1e-5
mx = 1 - mn

expmnt :: Experiment
expmnt = Experiment "probability" "multivariate"

isoexpmnt :: Maybe Analysis
isoexpmnt = Just $ Analysis "isolines" "dirichlet"

sgdexpmnt :: Maybe Analysis
sgdexpmnt = Just $ Analysis "descent" "dirichlet"


-- CSV

newtype DirichletSGD = DirichletSGD
    { ascent :: Double }
    deriving (Generic, Show)

instance ToNamedRecord DirichletSGD where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered DirichletSGD where
    headerOrder = goalCSVOrder

-- Training

eps :: Double
eps = 0.01

nsmps :: Int
nsmps = 10

nepchs :: Int
nepchs = 5000

drch0 :: Natural # Dirichlet 3
drch0 = Point $ S.fromTuple (1,1,1)

-- Functions

fitDirichlet
    :: Sample (Dirichlet 3)
    -> [Natural # Dirichlet 3]
fitDirichlet xyzs =
     vanillaGradientSequence (logLikelihoodDifferential xyzs) eps defaultAdamPursuit drch0

density2d
    :: Natural # Dirichlet 3
    -> (Double,Double)
    -> Double
density2d drch (x,y) =
    let z = 1 - x - y
     in if x + y < 0.995
           then density drch $ S.fromTuple (x,y,z)
           else 0


-- Main --


main :: IO ()
main = do


    -- Fit --

    xyzs <- realize $ sample nsmps tru

    let drchs = take nepchs $ fitDirichlet xyzs
        csts = logLikelihood xyzs <$> drchs

    -- Simple Statistics --

    let mnxs = S.map (/S.sum alphs) alphs
        sdxs = let alph0 = S.sum alphs
                   f alphi = alphi * (alph0 - alphi) / (square alph0 * (alph0 + 1))
                in S.map f alphs

    putStrLn "Dirichlet Means:"
    print mnxs
    putStrLn "Dirichlet Mean Logs:"
    print . coordinates $ toMean tru
    putStrLn "Dirichlet Variances:"
    print sdxs

    let dxyzs drch = do
            x <- range mn mx 100
            y <- range mn mx 100
            return (x,y,density2d drch (x,y))

        trups = dxyzs tru
        lrnps = dxyzs $ last drchs

    --let

    goalExport True expmnt isoexpmnt $ S.toList <$> xyzs

    goalExport False expmnt isoexpmnt trups

    goalExport False expmnt isoexpmnt lrnps

    goalExportNamed True expmnt sgdexpmnt $ DirichletSGD <$> csts

    runGnuplot expmnt isoexpmnt defaultGnuplotOptions "multivariate.gpi"
    runGnuplot expmnt sgdexpmnt defaultGnuplotOptions "cross-entropy-descent.gpi"

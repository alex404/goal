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
expmnt = Experiment "probability" "dirichlet"

-- CSV

newtype DirichletSGD = DirichletSGD
    { crossEntropy :: Double }
    deriving (Generic, Show)

instance ToNamedRecord DirichletSGD where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered DirichletSGD where
    headerOrder = goalCSVOrder

-- Training

eps :: Double
eps = -0.01

nsmps :: Int
nsmps = 50

nepchs :: Int
nepchs = 5000

drch0 :: Natural # Dirichlet 3
drch0 = Point $ S.fromTuple (1,1,1)

-- Functions

fitDirichlet
    :: Sample (Dirichlet 3)
    -> [Natural # Dirichlet 3]
fitDirichlet xyzs =
    let backprop drch = joinTangentPair drch $ stochasticCrossEntropyDifferential xyzs drch
     in vanillaGradientSequence backprop eps defaultAdamPursuit drch0

density2d
    :: Natural # Dirichlet 3
    -> (Double,Double)
    -> Double
density2d drch (x,y) =
    let z = 1 - x - y
     in if x + y < 0.995
           then density drch $ S.fromTuple (x,y,z)
           else 0

dirichletIsolines
    :: Natural # Dirichlet 3
    -> ([Isolevel],[[[(Double,Double)]]])
dirichletIsolines drch =
    let f x y = density2d drch (x,y)
     in isolines (mn,mx,100) (mn,mx,100) 20 f
     --   (lvls,lns) = isolines (mn,mx,100) (mn,mx,100) 20 f
     --in (lvls, [ [ln ++ [head ln]] | ln <- head <$> lns ])


-- Main --


main :: IO ()
main = do


    -- Fit --

    xyzs <- realize $ sample nsmps tru

    let drchs = take nepchs $ fitDirichlet xyzs
        csts = stochasticCrossEntropy xyzs <$> drchs

    -- Simple Statistics --

    let mnxs = S.map (/S.sum alphs) alphs
        sdxs = let alph0 = S.sum alphs
                   f alphi = alphi * (alph0 - alphi) / (square alph0 * (alph0 + 1))
                in S.map f alphs

    putStrLn "Dirichlet Means:"
    print mnxs
    putStrLn "Dirichlet Mean Logs:"
    print . coordinates . potentialDifferential $ toNatural tru
    putStrLn "Dirichlet Variances:"
    print sdxs

    let (trulvls,trulns) = dirichletIsolines tru
    let (sgdlvls,sgdlns) = dirichletIsolines $ last drchs

    --let

    goalExportNamed True expmnt Nothing trulvls
    mapM_ (goalExportLines False expmnt Nothing) trulns

    goalExportNamed False expmnt Nothing sgdlvls
    mapM_ (goalExportLines False expmnt Nothing) sgdlns

    goalExport False expmnt Nothing $ S.toList <$> xyzs

    goalExportNamed False expmnt Nothing $ DirichletSGD <$> csts

    runGnuplot expmnt Nothing defaultGnuplotOptions "multivariate.gpi"
    runGnuplot expmnt Nothing defaultGnuplotOptions "cross-entropy-descent.gpi"

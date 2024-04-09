{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Globals ---

--- Model

scattru :: Source # Categorical 2
scattru = fromTuple (0.15, 0.25)

mcattru :: Mean # Categorical 2
mcattru = toMean scattru

nx :: Natural # Categorical 2
nx = 0

nxz :: Natural # Tensor (Categorical 2) (Dirichlet 3)
nxz = fromRows $ S.fromTuple (fromTuple (-1, 1, 0), fromTuple (-1, 0, 1))

lkl :: Natural # Categorical 2 <* Dirichlet 3
lkl = join nx nxz

drch0 :: Natural # Dirichlet 3
drch0 = 0.1

rho :: Natural # Dirichlet 3
rho = fromTuple (-1, 0, 0)

--- Simplices

vrt1x, vrt1y, vrt2x, vrt2y, vrt3x, vrt3y :: Double
vrt1x = 0
vrt1y = 0
vrt2x = 0.5
vrt2y = sqrt 3 / 2
vrt3x = 1
vrt3y = 0

barycentricToCartesian :: Mean # Categorical 2 -> (Double, Double)
barycentricToCartesian mcat =
    let bry = coordinates mcat
        x2 = S.unsafeIndex bry 0
        x3 = S.unsafeIndex bry 1
        x1 = 1 - x2 - x3
     in (vrt1x * x1 + vrt2x * x2 + vrt3x * x3, vrt1y * x1 + vrt2y * x2 + vrt3y * x3)

cartesianToBarycentric :: (Double, Double) -> (Double, Double, Double)
cartesianToBarycentric (x, y) =
    let det = (vrt2y - vrt3y) * (vrt1x - vrt3x) + (vrt3x - vrt2x) * (vrt1y - vrt3y)
        x1 = ((vrt2y - vrt3y) * (x - vrt3x) + (vrt3x - vrt2x) * (y - vrt3y)) / det
        x2 = ((vrt3y - vrt1y) * (x - vrt3x) + (vrt1x - vrt3x) * (y - vrt3y)) / det
        x3 = 1 - x1 - x2
     in (x1, x2, x3)

--- Simulation Parameters

nobs :: Int
nobs = 10

pltres :: Int
pltres = 200

pltrngx, pltrngy :: [Double]
pltrngx = range 0 1 pltres
pltrngy = range 0 vrt2y pltres

pltxs :: [(Double, Double)]
pltxs = [(x, y) | y <- pltrngy, x <- pltrngx]

--- Helper functions

bounded :: Double -> Bool
bounded x = x > 0 && x < 1

density2d :: Natural # Dirichlet 3 -> (Double, Double, Double) -> Double
density2d drch (x1, x2, x3) =
    if all bounded [x1, x2, x3]
        then density drch $ S.fromTuple (x1, x2, x3)
        else 0

batchInference ::
    Natural # Dirichlet 3 ->
    Random (Sample (Categorical 2), Natural # Dirichlet 3)
batchInference drch = do
    xs <- sample nobs scattru
    return (xs, last $ conjugatedRecursiveBayesianInference0 rho lkl drch xs)

averageObservations :: Sample (Categorical 2) -> (Double, Double)
averageObservations =
    barycentricToCartesian . averageSufficientStatistic

-- -- Emission Distribution
--
main :: IO ()
main = do
    (xs1, drch1) <- realize $ batchInference drch0
    (xs2, drch2) <- realize $ batchInference drch1
    (xs3, drch3) <- realize $ batchInference drch2

    let drch0dnss = density2d drch0 . cartesianToBarycentric <$> pltxs
        drch1dnss = density2d drch1 . cartesianToBarycentric <$> pltxs
        drch2dnss = density2d drch2 . cartesianToBarycentric <$> pltxs
        drch3dnss = density2d drch3 . cartesianToBarycentric <$> pltxs

    let observations = [xs1, xs2, xs3]

    let njson =
            toJSON
                [ "true-categorical" .= barycentricToCartesian mcattru
                , "observations" .= observations
                , "average-observations" .= map averageObservations observations
                , "dirichlets" .= map coordinates [drch0, drch1, drch2, drch3]
                , "plot-xs" .= pltrngx
                , "plot-ys" .= pltrngy
                , "vertices" .= [(vrt1x, vrt1y), (vrt2x, vrt2y), (vrt3x, vrt3y)]
                , "dirichlet-densities" .= [drch0dnss, drch1dnss, drch2dnss, drch3dnss]
                ]

    --- Process data
    mvnfl <- resultsFilePath "categorical-inference.json"

    exportJSON mvnfl njson

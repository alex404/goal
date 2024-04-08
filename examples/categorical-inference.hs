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
scattru = fromTuple (0.1, 0.2)

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

--- Simulation Parameters

nobs :: Int
nobs = 10

pltres :: Int
pltres = 100

pltxys :: [(Double, Double)]
pltxys = [(x, y) | x <- range 0 1 pltres, y <- range 0 1 pltres]

--- Helper functions

density2d :: Natural # Dirichlet 3 -> (Double, Double) -> Double
density2d drch (x, y) =
    let z = 1 - x - y
     in if x + y < 0.995
            then density drch $ S.fromTuple (x, y, z)
            else 0

batchInference ::
    Natural # Dirichlet 3 ->
    Random (Sample (Categorical 2), Natural # Dirichlet 3)
batchInference drch = do
    xs <- sample nobs scattru
    return (xs, last $ conjugatedRecursiveBayesianInference0 rho lkl drch xs)

averageObservations :: Sample (Categorical 2) -> Mean # Categorical 2
averageObservations = averageSufficientStatistic

-- -- Emission Distribution
--
main :: IO ()
main = do
    (xs1, drch1) <- realize $ batchInference drch0
    (xs2, drch2) <- realize $ batchInference drch1
    (xs3, drch3) <- realize $ batchInference drch2

    let drch0dnss = density2d drch0 <$> pltxys
        drch1dnss = density2d drch1 <$> pltxys
        drch2dnss = density2d drch2 <$> pltxys
        drch3dnss = density2d drch3 <$> pltxys

    let observations = [xs1, xs2, xs3]

    let njson =
            toJSON
                [ "true-categorical" .= coordinates scattru
                , "observations" .= observations
                , "average-observations" .= map (coordinates . averageObservations) observations
                , "dirichlets" .= map coordinates [drch0, drch1, drch2, drch3]
                , "dirichlet-densities" .= [drch0dnss, drch1dnss, drch2dnss, drch3dnss]
                ]

    --- Process data
    mvnfl <- resultsFilePath "categorical-inference.json"

    exportJSON mvnfl njson

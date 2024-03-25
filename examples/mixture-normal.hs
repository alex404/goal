{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Globals ---

pltsmps :: Int
pltsmps = 100

--- Mixture Model

mu1, mu2, mu3 :: Double
mu1 = 0
mu2 = 2
mu3 = -1

mixnrm1, mixnrm2, mixnrm3 :: Source # Normal
mixnrm1 = fromTuple (mu1, 1)
mixnrm2 = fromTuple (mu2, 0.5)
mixnrm3 = fromTuple (mu3, 2)

wghts :: Source # Categorical 2
wghts = fromTuple (0.2, 0.3)

nwghts :: Natural # Categorical 2
nwghts = toNatural wghts

mix :: Source # Mixture Normal 2
mix = joinSourceMixture (S.fromTuple (mixnrm1, mixnrm2, mixnrm3)) wghts

lkl :: Natural # Normal <* Categorical 2
lkl = fst . split $ toNatural mix

mixobs :: [S.Vector 1 Double]
mixobs = S.singleton <$> [-1, 0.5, 1.5]

mixpstrs :: [Natural # Categorical 2]
mixpstrs = conjugatedBayesRule lkl nwghts <$> mixobs

--- Main ---

main :: IO ()
main = do
    --- Mixtures
    let mixrng = range (-4) 4 pltsmps
        mixrng' = S.singleton <$> mixrng
        mxdns = observableDensity (toNatural mix) <$> mixrng'
        nrmdnss = [densities mixnrm mixrng' | mixnrm <- [mixnrm1, mixnrm2, mixnrm3]]

    --- Export
    let json =
            toJSON
                [ "plot-samples" .= mixrng
                , "mixture-density" .= mxdns
                , "component-densities" .= nrmdnss
                , "weights" .= S.toList (categoricalWeights wghts)
                , "observations" .= mixobs
                , "posteriors" .= (S.toList . categoricalWeights <$> mixpstrs)
                ]

    --- Process data
    flnm <- resultsFilePath "mixture-normal.json"

    exportJSON flnm json

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedStrings #-}

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Graphical
import Goal.Probability

import Goal.Core.Vector.Storable qualified as S

--- Misc

import Data.List qualified as L
import Data.Proxy (Proxy (..))

--- Globals ---

pltsmps :: Int
pltsmps = 100

--- Factor Analysis

-- Definitions
famux, famuy, favrx, favry :: Double
famux = 1
famuy = -1
favrx = 3
favry = 2

fasx :: Source # DiagonalNormal 2
fasx = fromTuple (famux, famuy, favrx, favry)

wmtx :: Source # Tensor (StandardNormal 2) (StandardNormal 1)
wmtx = fromTuple (2, 1)

falkl :: Source # FactorAnalysis 2 1
falkl = join fasx wmtx

faprr :: Source # Normal
faprr = standardNormal

faobs :: S.Vector 2 Double
faobs = S.fromTuple (0, 1)

fapstr :: Natural # Normal
fapstr = conjugatedBayesRule (toNatural falkl) (toNatural faprr) faobs

--- Mixture Model

mixnrm1, mixnrm2, mixnrm3 :: Source # Normal
mixnrm1 = fromTuple (0, 1)
mixnrm2 = fromTuple (2, 0.5)
mixnrm3 = fromTuple (-1, 2)

mixprr :: Source # Categorical 2
mixprr = fromTuple (0.2, 0.3)

mix :: Source # Mixture Normal 2
mix = joinSourceMixture (S.fromTuple (mixnrm1, mixnrm2, mixnrm3)) mixprr

mixlkl :: Natural # Normal <* Categorical 2
mixlkl = fst . split $ toNatural mix

mixobs :: S.Vector 1 Double
mixobs = S.singleton 0

mixpstr :: Natural # Categorical 2
mixpstr = conjugatedBayesRule (toNatural mixlkl) (toNatural mixprr) mixobs

--- Von Mises Population Code

type Neurons n = Replicated n Poisson
type N = 10

ppcprr :: Source # VonMises
ppcprr = fromTuple (pi, 1)

joinVonMisesIndependent ::
    (KnownNat n, LegendreExponentialFamily x) =>
    -- | Gains
    Natural # Neurons n ->
    -- | Von Mises Curves
    S.Vector n (Natural # x) ->
    -- | Population Likelihood
    Natural # Neurons n <* x
joinVonMisesIndependent nz0 nps =
    let mtx = fromRows nps
        nz = nz0 - Point (S.map potential nps)
     in join nz mtx

kp :: Double
kp = 1

mus :: S.Vector N Double
mus = S.generate (\k -> 2 * pi * fromIntegral k / fromIntegral (natVal (Proxy @N)))

tcs :: S.Vector N (Source # VonMises)
tcs = S.map (fromTuple . (,kp)) mus

gns :: Source # Neurons N
gns = Point $ S.replicate 2

ppclkl :: Natural # Neurons 10 <* VonMises
ppclkl = joinVonMisesIndependent (toNatural gns) (S.map toNatural tcs)

--- Main ---

main :: IO ()
main = do
    --- Factor Analysis
    let z0 :: S.Vector 1 Double
        z0 = 0

        obsrngs :: ([Double], [Double])
        obsrngs =
            let sdx = sqrt favrx
                sdy = sqrt favry
                (mux0, muy0) = S.toPair . coordinates . fst . split . toSource $ toNatural falkl >.>* z0
                xrng = range (mux0 - 3 * sdx) (mux0 + 3 * sdx) pltsmps
                yrng = range (muy0 - 3 * sdy) (muy0 + 3 * sdy) pltsmps
             in (xrng, yrng)

        ltntrng :: [Double]
        ltntrng = range (-3) 3 pltsmps

    let (obsxs, obsys) = obsrngs
        lklmvn = toNatural falkl >.>* z0

    let falkldns = do
            x <- obsxs
            y <- obsys
            return (x, y, density lklmvn $ S.fromTuple (x, y))
        faprrdns = [(x, density faprr $ S.singleton x) | x <- ltntrng]
        fapstrdns = [(x, density fapstr $ S.singleton x) | x <- ltntrng]

    --- Mixtures
    let mixrng = range (-4) 4 pltsmps
        mixrng' = S.singleton <$> mixrng
        mxdns = observableDensity (toNatural mix) <$> mixrng'
        nrm1dns = density mixnrm1 <$> mixrng'
        nrm2dns = density mixnrm2 <$> mixrng'
        nrm3dns = density mixnrm3 <$> mixrng'

    --- Population Code
    let stmsmps = range 0 (2 * pi) pltsmps
        nrnss = ppclkl >$>* stmsmps
        frtss = L.transpose $ listCoordinates . toMean <$> nrnss
        ppcprrdns = density ppcprr <$> stmsmps
        stm0 = pi

    ppxsmp <- realize $ samplePoint (ppclkl >.>* stm0)

    let ppcpstr :: Natural # VonMises
        ppcpstr = conjugatedBayesRule ppclkl (toNatural ppcprr) ppxsmp

    --- Export
    let json =
            toJSON
                [ "fa-observable-ranges" .= obsrngs
                , "fa-latent-range" .= ltntrng
                , "fa-likelihood" .= falkldns
                , "fa-prior" .= faprrdns
                , "fa-observation" .= faobs
                , "fa-posterior" .= fapstrdns
                , "mix-rnage" .= mixrng
                , "mix-densities" .= mxdns
                , "mix-component-1-densities" .= nrm1dns
                , "mix-component-2-densities" .= nrm2dns
                , "mix-component-3-densities" .= nrm3dns
                , "mix-prior" .= S.toList (categoricalWeights mixprr)
                , "mix-observation" .= mixobs
                , "mix-posterior" .= S.toList (categoricalWeights mixpstr)
                , "stimulus-samples" .= stmsmps
                , "neuron-responses" .= frtss
                , "von-mises-density" .= ppcprrdns
                ]

    --- Process data
    flnm <- resultsFilePath "conjugated-harmoniums.json"

    exportJSON flnm json

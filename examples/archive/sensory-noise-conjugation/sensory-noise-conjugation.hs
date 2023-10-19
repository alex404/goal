#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L


--- Globals ---

-- Types
type N = 200
type Neurons = Replicated N Poisson

-- Ranges
mnx,mxx :: Double
mnx = 0
mxx = 2*pi

xs :: [Double]
xs = range mnx mxx 200

mny,mxy :: Double
mny = 0
mxy = 2*pi

ys :: [Double]
ys = xs

stm :: Double
stm = pi

-- Functions
--joinVonMisesIndependent
--    :: Natural # Neurons -- ^ Gains
--    -> S.Vector N (Natural # VonMises) -- ^ Von Mises Curves
--    -> Natural # Neurons <* VonMises -- ^ Population Likelihood
--joinVonMisesIndependent nz0 nps =
--    let mtx = fromRows nps
--        nz = nz0 - Point (S.map potential nps)
--     in join nz mtx


-- Sensory Noise
kp :: Double
kp = 10

gyx :: Double -> Source # VonMises
gyx x = fromTuple (x,kp)

-- Population Encoder
--rtc :: Random (Natural # VonMises)
--rtc = uniformInitialize (-1,1)


rfzy :: Random (Natural # Neurons <* VonMises)
rfzy = uniformInitialize (-1,1)

--rfzy :: Random (Natural # Neurons <* VonMises)
--rfzy = do
--    tcs <- S.replicateM rtc
--    return $ joinVonMisesIndependent 0 tcs

sx :: Double -> S.Vector 3 Double
sx x = S.fromTuple (1,cos x, sin x)



--- Main ---


main :: IO ()
main = do

    -- Source Population
    fzy <- realize rfzy

    let stcs :: [Double]
        stcs = potential <$> fzy >$>* ys

    let bts :: S.Vector 3 Double
        bts = S.linearLeastSquares (sx <$> xs) stcs

        stchts :: [Double]
        stchts = S.dotMap bts $ sx <$> ys

    let fys :: [[Double]]
        fys = L.transpose $ listCoordinates . toMean <$> fzy >$>* ys

    z0 <- realize . samplePoint $ fzy >.>* stm

    let unnormalizedLogPosterior1 z x = logDensity (fzy >.>* x) z

    let lgpsts1 = unnormalizedLogPosterior1 z0 <$> xs
        elgpsts1 = exp <$> lgpsts1

    let pstbts1 :: S.Vector 3 Double
        pstbts1 = S.linearLeastSquares (sx <$> xs) lgpsts1

        lgpsthts1 :: [Double]
        lgpsthts1 = S.dotMap pstbts1 $ sx <$> xs
        elgpsthts1 = exp <$> lgpsthts1

        --prt1 = logIntegralExp 1e-4 (\x -> S.dotProduct pstbts1 (sx x)) mny mxy ys

    -- Information limited
    let code0 x z y = logDensity (fzy >.>* y) z + logDensity (gyx x) y
        unnormalizedLogPosterior2 z x = logIntegralExp 1e-4 (code0 x z) mny mxy ys
        --posteriorPartition z = logIntegralExp 1e-2 (unnormalizedLogPosterior z) mny mxy ys
        --logPosterior z x = unnormalizedLogPosterior z x - posteriorPartition z


    let lgpsts2 = unnormalizedLogPosterior2 z0 <$> xs
        elgpsts2 = exp <$> lgpsts2

    let pstbts2 :: S.Vector 3 Double
        pstbts2 = S.linearLeastSquares (sx <$> xs) lgpsts2

        lgpsthts2 :: [Double]
        lgpsthts2 = S.dotMap pstbts2 $ sx <$> xs
        elgpsthts2 = exp <$> lgpsthts2

        --prt2 = logIntegralExp 1e-4 (\x -> S.dotProduct pstbts2 (sx x)) mny mxy ys

    goalExport "." "sensory-noise" . zip xs $ densities (gyx stm) xs
    runGnuplot "." "sensory-noise"

    goalExport "." "conjugation-curve" $ zip3 ys stcs stchts
    putStrLn "(rho0, rho1, rho2):"
    print bts
    runGnuplot "." "conjugation-curve"

    goalExport "." "population-coding" . L.transpose $ ys:fys
    runGnuplot "." "population-coding"

    goalExport "." "source-code" $ L.zip4 xs elgpsts1 elgpsthts1 (zipWith subtract elgpsts1 elgpsthts1)
    goalExport "." "limited-code" $ L.zip4 xs elgpsts2 elgpsthts2 (zipWith subtract elgpsts2 elgpsthts2)
    runGnuplot "." "limited-code"



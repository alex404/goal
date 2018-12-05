{-# LANGUAGE GADTs,FlexibleContexts,TypeOperators,DataKinds #-}

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Plot

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

prjdr :: String
prjdr = "probability/von-mises-mixture"

-- Manifolds --

type Latent = Categorical Int 3
type Observable = (VonMises, VonMises)
type Harmonium' = Harmonium Tensor Observable Latent

-- Mixture Distributions --

mux1,mux2,mux3,muy1,muy2,muy3 :: Double
mux1 = 1.5
muy1 = 1.5
mux2 = 3
muy2 = 3
mux3 = 4.5
muy3 = 4.5

kpx1,kpx2,kpx3,kpy1,kpy2,kpy3 :: Double
kpx1 = 2
kpy1 = 2
kpx2 = 2
kpy2 = 2
kpx3 = 2
kpy3 = 2

vm1,vm2,vm3 :: Source # Observable
vm1 = Point $ S.fromTuple (mux1, kpx1, muy1, kpy1)
vm2 = Point $ S.fromTuple (mux2, kpx2, muy2, kpy2)
vm3 = Point $ S.fromTuple (mux3, kpx3, muy3, kpy3)

vms :: S.Vector 3 (Natural # Observable)
vms = S.fromTuple (toNatural vm1,toNatural vm2,toNatural vm3)

mix1,mix2 :: Double
mix1 = 0.2
mix2 = 0.5

wghts :: Source # Latent
wghts = Point $ S.doubleton mix1 mix2

truhrm :: Natural # Harmonium Tensor Observable Latent
truhrm = buildMixtureModel vms $ toNatural wghts

vm :: Source # VonMises
vm = Point $ S.fromTuple (2,2)

nvm :: Natural # VonMises
nvm = toNatural vm

naff :: Mean #> Natural # Observable <* Latent
nx0 :: Natural # OneHarmonium Latent
(naff,nx0) = splitBottomHarmonium truhrm

nz :: Natural # Observable
nzx :: Mean #> Natural # Tensor Observable Latent
(nz,nzx) = splitAffine naff


--- Main ---


main :: IO ()
main = do

    zxs <- realize $ sample 10000 truhrm

    let mzxs = dualTransition truhrm

    let mzxs' :: Mean # Harmonium Tensor Observable Latent
        mzxs' = sufficientStatisticT zxs

    let udensity :: Double -> Double -> Double
        udensity x y = unnormalizedHarmoniumObservableDensity truhrm (x,y)
        -- udensity x y = exp (mixtureModelLogLikelihood truhrm (x,y)) / (2*pi)
        nrm = fst $ integrate 1e-12 (\x -> fst $ integrate 1e-12 (udensity x) 0 (2*pi)) 0 (2*pi)
        fooDens (x,y) = udensity x y / nrm
        squaredError xy = square $ mixtureDensity0 truhrm xy - fooDens xy
        err = sqrt . average $ squaredError . hHead <$> zxs

    let udensity' x y = sum [ unnormalizedDensity truhrm ((x,y) :+: 0 :+: Null)
                            , unnormalizedDensity truhrm ((x,y) :+: 1 :+: Null)
                            , unnormalizedDensity truhrm ((x,y) :+: 2 :+: Null) ]
        nrm' = fst $ integrate 1e-12 (\x -> fst $ integrate 1e-12 (udensity' x) 0 (2*pi)) 0 (2*pi)
        fooDens' (x,y) = udensity' x y / nrm'
        squaredError' xy = square $ mixtureDensity0 truhrm xy - fooDens' xy
        err' = sqrt . average $ squaredError' . hHead <$> zxs
        squaredError'' xy = square $ mixtureDensity0 truhrm xy - mixtureDensity truhrm xy
        err'' = sqrt . average $ squaredError'' . hHead <$> zxs

        (rho0,rprms) = mixtureLikelihoodRectificationParameters naff

    putStrLn "True:"
    print mzxs
    putStrLn "Estimate:"
    print mzxs'
    putStrLn "Analytic Log-Posterior:"
    print $ potential truhrm
    putStrLn "Rectified Log-Posterior:"
    print $ rho0 + potential (rprms <+> fromOneHarmonium nx0)
    putStrLn "Density Error Final:"
    print err'

    putStrLn "Density Error:"
    print err
    putStrLn "Integral Log-Posterior:"
    print $ log nrm

    putStrLn "Density Error 2:"
    print err'
    putStrLn "Integral Log-Posterior 2:"
    print $ log nrm'

    putStrLn "VM1 Log-Posterior:"
    print $ potential nvm
    putStrLn "VM1 Integral Posterior:"
    print . log . fst $ integrate 1e-12 (unnormalizedDensity nvm) 0 (2*pi)





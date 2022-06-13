#! stack runghc

{-# LANGUAGE TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

kp :: Double
kp = 1

rotationMatrix :: Double -> Natural # Tensor VonMises VonMises
rotationMatrix mu = kp .> fromTuple (cos mu,-sin mu,sin mu,cos mu)

rt :: Double
rt = 1/2 * pi

nxz :: Natural # Tensor VonMises VonMises
nxz = rotationMatrix rt
--nxz = fromTuple (-1,1,-1,2)

nx :: Natural # VonMises
nx = fromTuple (0,5)

nz :: Natural # VonMises
nz = fromTuple (-2,2)

pstr :: Natural # Affine Tensor VonMises VonMises VonMises
pstr = join nz $ transpose nxz

lkl :: Natural # Affine Tensor VonMises VonMises VonMises
lkl = join nx nxz

xs :: [Double]
xs = range 0 (2*pi) 1000

ptns :: [Double]
ptns = potential <$> lkl >$>* xs

rgrs :: Double -> S.Vector 3 Double
rgrs x = S.fromTuple (1,cos x, sin x)

bts :: S.Vector 3 Double
bts = S.linearLeastSquares (rgrs <$> xs) ptns

ptnhts :: [Double]
ptnhts = S.dotMap bts $ rgrs <$> xs

vmhLogPartition :: Double
vmhLogPartition =
    let bm = square $ 1/(2*pi)
        err = 1e-3
        krn x z =
            let sx = sufficientStatistic x
                sz = sufficientStatistic z
             in (*bm) . exp $ nx <.> sx + nz <.> sz + sx <.> (nxz >.> sz)
        krn' x = fst $ integrate err (krn x) 0 (2*pi)
     in log . fst $ integrate err krn' 0 (2*pi)


priorDensities :: [Double]
priorDensities =
    let mxs = sufficientStatistic <$> xs
        nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
        lgprt = vmhLogPartition
     in map (exp . subtract lgprt) . zipWith (+) nrgs $ logBaseMeasure (Proxy @VonMises) <$> xs

--- Main ---


main :: IO ()
main = do

    let dnss = priorDensities
        ldnss = densities (lkl >.> pi) xs
    goalExport "." "conjugation-curve" $ zip3 xs ptns ptnhts
    goalExport "." "prior-density" $ zip xs dnss
    goalExport "." "likelihood-density" $ zip xs ldnss
    putStrLn "(rho0, rho1, rho2):"
    print bts
    runGnuplot "." "conjugation-curve"
    runGnuplot "." "densities"

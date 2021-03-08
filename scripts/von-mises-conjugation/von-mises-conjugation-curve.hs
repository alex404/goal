#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S


--- Globals ---

kp :: Double
kp = 2

nzx :: Natural # Tensor VonMises VonMises
nzx = 1 * fromTuple (2,0,0,2)

nz,nx :: Natural # VonMises
nz = fromTuple (2,2)
nx = fromTuple (0,0)

fzx :: Natural # Affine Tensor VonMises VonMises VonMises
fzx = join nz nzx

vmhrm :: Natural # Harmonium Tensor VonMises VonMises VonMises VonMises
vmhrm = join fzx nx

xs :: [Double]
xs = range 0 (2*pi) 1000

ys :: [Double]
ys = potential <$> fzx >$>* xs

sx :: Double -> S.Vector 3 Double
sx x = S.fromTuple (1,cos x, sin x)

bts :: S.Vector 3 Double
bts = S.linearLeastSquares (sx <$> xs) ys

yhts :: [Double]
yhts = S.dotMap bts $ sx <$> xs

--- Main ---


main :: IO ()
main = do


    goalExport "." "conjugation-curve" $ zip3 xs ys yhts
    putStrLn "(rho0, rho1, rho2):"
    print bts
    runGnuplot "." "conjugation-curve"

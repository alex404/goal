{-# LANGUAGE FlexibleContexts #-}


-- Goal --

import Pendulum

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation
import Goal.Cognition

-- Qualified --

import qualified Data.Vector.Unboxed as V
import qualified Statistics.Correlation as S

--- Program ---


-- Globals --

ntrls = 1000
nz1 = 10
nz2 = 20
nstp = 100
x0 = (1,0)
mlpfl = "multilayer-perceptron"
dtnm = mlpfl ++ "ef1"

-- Functions --

correlation xs0 ys0 =
    S.pearson . V.fromList $ zip xs0 ys0

-- Main --

main :: IO ()
main = do

    bl <- goalDoesFileExist sbdr dtnm
    g <- if bl
              then fromList mg1 . read <$> goalReadFile sbdr dtnm
              else error $ "Script requires a " ++ dtnm ++ " file"

    xchn <- runWithSystemRandom (pendulumChain x0)
    nmlys <- sequence $ replicate ntrls (runWithSystemRandom pendulumResponseMealy)
    let gflt = harmoniumEncodedFilter y01 amtx1 bmtx1 g
        zchns = [xchn >>> nmly >>> gflt | nmly <- nmlys ]
        zs = [ listCoordinates $ streamChain zchn !! nstp | zchn <- zchns ]
        zs1 = (!! nz1) <$> zs
        zs2 = (!! nz2) <$> zs
        crs = [ correlation zs1 $ (!! n) <$> zs | n <- [0..39] ]
    print $ crs
    print . average $ abs <$> crs
    return ()


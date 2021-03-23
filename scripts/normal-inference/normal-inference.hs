#! stack runghc

{-# LANGUAGE DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds #-}


--- Imports ---


import Goal.Geometry
import Goal.Probability
import Goal.Graphical


--- Globals ---

-- Manifolds --

-- Mixture Distributions --

psrx :: Source # Normal
psrx = fromTuple (0,10)

prx :: Natural # Normal
prx = toNatural psrx

-- Emission Distribution

scl,vr :: Double
scl = 1
vr = 5

enzx :: Natural # Tensor NormalMean NormalMean
enzx = singleton $ scl / vr

enrz :: Natural # Normal
enrz = fromTuple (1,-1/(2*vr))

efzx :: Natural # Affine Tensor NormalMean Normal NormalMean
efzx = join enrz enzx

main :: IO ()
main = do

    x0 <- realize $ samplePoint prx
    zs <- realize . sample 30 $ efzx >.>* x0

    let nblfs = conjugatedRecursiveBayesianInference efzx prx zs
        sblfs = toSource <$> nblfs
    print x0
    mapM_ print sblfs


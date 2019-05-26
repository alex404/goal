{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    ScopedTypeVariables,
    DataKinds
    #-}


--- Imports ---


import NeuralData
import NeuralData.VonMises
import NeuralData.Mixture

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

--- Globals ---


nbns :: Int
nbns = 10

xsmps :: [Double]
xsmps = init $ range mnx mxx 101

-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
mixtureStochasticConditionalCrossEntropy' xs ys f =
    let nys = f >$>* xs
     in negate . average $ log <$> zipWith mixtureDensity' nys ys

mixtureDensity' hrm x =
    let (ncmpnts,nwghts) = splitMixture hrm
        dnss = (`density` x) <$> S.toList ncmpnts
        wghts = S.toList $ categoricalWeights nwghts
     in weightedAverage $ zip wghts dnss



--- CLI ---


main :: IO ()
main = do

    let expnm = "synthetic-3k-500n"
        dst = "random"

    mtrulkl0 <- getMixtureLikelihood expnm dst

    zxs0  <- snd <$> getNeuralData expnm dst

    let (k,m,cs) = fromJust mtrulkl0

    let mx = 8

    case someNatVal k of
        SomeNat (Proxy :: Proxy k) -> case someNatVal m
            of SomeNat (Proxy :: Proxy m) -> do

                let zxs :: [(Response k, Double)]
                    zxs = strengthenNeuralData zxs0

                let inteps = -0.05
                    intnepchs = 500
                    intnbtch = 100

                let lkl = last $ fitIPLikelihood inteps intnbtch intnepchs zxs

                let mlkl :: Mean #> Natural # MixtureGLM (Neurons k) m VonMises
                    mlkl = strengthenMixtureLikelihood cs

                let (zs,xs) = unzip zxs

                    mxmdl = mlkl >.>* pi
                    nlmdl = lkl >.>* pi

                let z0 = [0..100]

                let zs' :: [Response k]
                    zs' = fromJust . B.fromList <$> replicateM (fromIntegral k) z0

                let ubnd = stochasticConditionalCrossEntropy xs zs lkl
                    lbnd = mixtureStochasticConditionalCrossEntropy xs zs mlkl
                    lbnd' = mixtureStochasticConditionalCrossEntropy' xs zs mlkl

                putStrLn $ "Cost upper bound: " ++ show ubnd
                putStrLn $ "Cost lower bound: " ++ show lbnd
                putStrLn $ "Cost lower bound prime: " ++ show lbnd'

                let dns = mixtureDensity (mlkl >.>* head xs) $ head zs
                let dns' = mixtureDensity' (mlkl >.>* head xs) $ head zs

                putStrLn $ "Density: " ++ show dns
                putStrLn $ "Density Prime: " ++ show dns'

                let avg = sum $ mixtureDensity mxmdl <$> zs'
                let avg' = sum $ mixtureDensity' mxmdl <$> zs'
                let avg0 = sum $ density nlmdl <$> zs'
                putStrLn $ "Average: " ++ show avg
                putStrLn $ "Average Prime: " ++ show avg'
                putStrLn $ "Average Null: " ++ show avg0

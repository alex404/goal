{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeOperators
    #-}

module NeuralData.Mixture
    ( -- * Mixture
      fitMixtureLikelihood
    , getFittedMixtureLikelihood
    , strengthenMixtureLikelihood
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Types ---


prjnm :: String
prjnm = "neural-data"


--- Inference ---


getFittedMixtureLikelihood
    :: String
    -> String
    -> IO (NatNumber,NatNumber,[Double])
getFittedMixtureLikelihood expnm dst = do
    (k,n,xs) <- read <$> goalReadDataset (Experiment prjnm expnm) dst
    return (k,n,xs)

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---


fitMixtureLikelihood
    :: forall k n r . (KnownNat k, KnownNat n)
    => [(Response k,Double)]
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises) -- ^ Function
fitMixtureLikelihood xzs = do
    let eps = -0.05
        nepchs = 500
    kps <- S.replicateM $ uniformR (0.2,0.6)
    let sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps $ S.range 0 (2*pi)
        wghts0 :: Natural # Categorical Int n
        wghts0 = zero
        wghts = toSource wghts0
        gnss0 = S.replicate . transition . sufficientStatisticT $ fst <$> xzs
    gnss' <- S.replicateM . fmap Point . S.replicateM $ uniformR (0,2)
    let gnss = S.zipWith (<+>) gnss0 gnss'
        ppc0 = vonMisesMixturePopulationEncoder True wghts gnss sps
        (zs,xs) = unzip xzs
        backprop mglm = joinTangentPair mglm $ mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm
    return (vanillaGradientSequence backprop eps defaultAdamPursuit ppc0 !! nepchs)

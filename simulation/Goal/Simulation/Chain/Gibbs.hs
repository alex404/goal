{-# LANGUAGE FlexibleContexts,Arrows #-}

-- | Gibbs sampling for exponential family harmoniums.
module Goal.Simulation.Chain.Gibbs where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Simulation.Circuit
import Goal.Simulation.Chain

--- Types ---

-- | Returns a Markov chain over the random variables in a deep harmonium by Gibbs sampling.
bulkGibbsChain
    :: ( Map Mean Natural f m n, KnownNat k, Generative Natural n, ExponentialFamily m
       , Bilinear f m n, Gibbs (f : fs) (n : m : ms), Manifold (DeepHarmonium fs (m : ms)) )
    => Natural # DeepHarmonium (f : fs) (n : m : ms) -- ^ The deep harmonium
    -> Sample k n -- ^ The initial states of the Gibbs chains
    -> Random s (Chain (Sample k (DeepHarmonium (f : fs) (n : m : ms)))) -- ^ The resulting Gibbs chains
{-# INLINE bulkGibbsChain #-}
bulkGibbsChain hrm z0s = do
    xzs0 <- initialPass hrm z0s
    gstp <- accumulateRandomFunction0 (gibbsPass hrm)
    return . accumulateCircuit xzs0 $ proc ((),xzs) -> do
        xzs' <- gstp -< xzs
        returnA -< (xzs,xzs')

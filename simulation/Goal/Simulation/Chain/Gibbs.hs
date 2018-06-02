{-# LANGUAGE FlexibleContexts,Arrows #-}

-- | Gibbs sampling for exponential family harmoniums.
module Goal.Simulation.Chain.Gibbs
    ( bulkGibbsChain
    , contrastiveDivergence
    ) where


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
    :: ( KnownNat k, Generative Natural m, ExponentialFamily n
       , Map Mean Natural f m n, Bilinear f m n
       , Gibbs (f : fs) (m : n : ms), Manifold (DeepHarmonium fs (n : ms)) )
    => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ The deep harmonium
    -> Sample k (DeepHarmonium (f : fs) (m : n : ms)) -- ^ The initial states of the Gibbs chains
    -> Random s (Chain (Sample k (DeepHarmonium (f : fs) (m : n : ms)))) -- ^ The resulting Gibbs chains
{-# INLINE bulkGibbsChain #-}
bulkGibbsChain hrm xzs0 = do
    gstp <- accumulateRandomFunction0 (gibbsPass hrm)
    return . accumulateCircuit xzs0 $ proc ((),xzs) -> do
        xzs' <- gstp -< xzs
        returnA -< (xzs,xzs')

contrastiveDivergence
    :: ( KnownNat k, 1 <= k, Generative Natural m, ExponentialFamily m, ExponentialFamily n
       , Map Mean Natural f m n, Bilinear f m n, Gibbs '[f] '[m,n] )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample k m -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f m n -- ^ The harmonium
      -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    gchn <- bulkGibbsChain hrm xzs0
    let xzs1 = streamChain gchn !! cdn
    return $ estimateStochasticCrossEntropyDifferential xzs0 xzs1


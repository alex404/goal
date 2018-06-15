{-# LANGUAGE FlexibleContexts,Arrows #-}

-- | Gibbs sampling for exponential family harmoniums.
module Goal.Simulation.Chain.Gibbs
    ( bulkGibbsChain
    , contrastiveDivergence
    , harmoniumInformationProjection
    -- , fitRectificationParameters
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Simulation.Circuit
import Goal.Simulation.Circuit.Optimization
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

informationProjectionDifferential0
    :: forall k f m n r
    . ( ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m
      , Generative Natural n, KnownNat k, 1 <= k, 2 <= k )
    => Proxy k
    -> Natural # Harmonium f m n
    -> Natural # n
    -> Random r (CotangentVector Natural n)
{-# INLINE informationProjectionDifferential0 #-}
informationProjectionDifferential0 _ hrm nx = do
    (xs :: Sample k n) <- sample nx
    return $ harmoniumInformationProjectionDifferential nx xs hrm

informationProjectionCircuit1
    :: ExponentialFamily n
    => Double
    -> Double
    -> Double
    -> Double
    -> Natural # n
    -> Circuit (Natural # n) (CotangentVector Natural n)
    -> Chain (Natural # n)
{-# INLINE informationProjectionCircuit1 #-}
informationProjectionCircuit1 eps bt1 bt2 rg nx0 dffcrc = accumulateCircuit0 nx0 $ proc ((),nx) -> do
    dnx <- dffcrc -< nx
    adamAscent eps bt1 bt2 rg -< breakPoint $ joinTangentPair nx dnx

harmoniumInformationProjection
    :: ( ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m
       , Generative Natural n, KnownNat k, 1 <= k, 2 <= k )
    => Proxy k
    -> Double
    -> Double
    -> Double
    -> Double
    -> Int
    -> Natural # Harmonium f m n
    -> Natural # n
    -> Random s (Natural # n)
{-# INLINE harmoniumInformationProjection #-}
harmoniumInformationProjection prxk eps bt1 bt2 rg nstps hrm nx0 = do
    dffcrc <- accumulateRandomFunction0 $ informationProjectionDifferential0 prxk hrm
    return $ streamChain (informationProjectionCircuit1 eps bt1 bt2 rg nx0 dffcrc) !! nstps

--class FitRectificationParameters (fs :: [* -> * -> *]) (ms :: [*]) where
--    fitRectificationParameters
--        :: Double
--        -> Maybe Int
--        -> Natural # DeepHarmonium fs ms
--        -> Natural # Sum (Tail ms)
--        -> Random s (Natural # Sum (Tail ms))
--
--instance FitRectificationParameters '[] '[m] where
--    {-# INLINE fitRectificationParameters #-}
--    fitRectificationParameters _ _ _ _ = zero
--
--instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
--         , ExponentialFamily n, SampleRectified fs (n : ms), Generative Natural m
--         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
--  => SampleRectified (f : fs) (m : n : ms) where
--    {-# INLINE sampleRectified #-}
--    sampleRectified rprms dhrm = do
--        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
--            (rprm,rprms') = splitSum rprms
--        (ys,xs) <- fmap hUnzip . sampleRectified rprms' $ biasBottom rprm dhrm'
--        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
--        return . hZip zs $ hZip ys xs
--
--

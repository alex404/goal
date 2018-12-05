{-# LANGUAGE
    Arrows,
    RankNTypes,
    TypeOperators,
    DataKinds,
    FlexibleContexts
#-}

-- | A collection of 'Circuit's for computing the differentials in gradient
-- descent algorithms.

module Goal.Simulation.Circuit.Optimization
    ( -- * Gradient Pursuit
      gradientCircuit
    -- * Monte Carlo
    , bulkGibbsChain
    , contrastiveDivergence
    , dualContrastiveDivergence
    , harmoniumInformationProjection
    ) where


--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Simulation.Chain
import Goal.Simulation.Circuit


-- | A 'Circuit' for classic gradient descent.
gradientCircuit
    :: Manifold m
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Circuit (TangentPair c m) (Point c m) -- ^ Gradient Ascent
{-# INLINE gradientCircuit #-}
gradientCircuit eps gp = accumulateFunction (repeat zero,0) $ \pdp (vs,k) ->
    let (p',vs') = gradientPursuitStep eps gp k pdp vs
     in (p',(vs',k+1))

-- | Returns a Markov chain over the random variables in a deep harmonium by Gibbs sampling.
bulkGibbsChain
    :: ( Generative Natural m, ExponentialFamily n, Map Mean Natural f m n, Bilinear f m n
       , Gibbs (f : fs) (m : n : ms), Manifold (DeepHarmonium fs (n : ms)) )
    => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ The deep harmonium
    -> Sample (DeepHarmonium (f : fs) (m : n : ms)) -- ^ The initial states of the Gibbs chains
    -> s ~> Chain (Sample (DeepHarmonium (f : fs) (m : n : ms))) -- ^ The resulting Gibbs chains
{-# INLINE bulkGibbsChain #-}
bulkGibbsChain hrm xzs0 = do
    gstp <- accumulateRandomFunction0 (gibbsPass hrm)
    return . accumulateCircuit xzs0 $ proc ((),xzs) -> do
        xzs' <- gstp -< xzs
        returnA -< (xzs,xzs')

contrastiveDivergence
    :: ( Generative Natural m, ExponentialFamily m, ExponentialFamily n
       , Generative Natural n, Map Mean Natural f m n, Bilinear f m n )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample m -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f m n -- ^ The harmonium
      -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    gchn <- bulkGibbsChain hrm xzs0
    let xzs1 = streamChain gchn !! cdn
    return $ stochasticCrossEntropyDifferential' xzs0 xzs1

dualContrastiveDivergence
    :: forall s f z x
    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
      , Map Mean Natural f x z, Bilinear f z x, Bilinear f x z )
      => Int -- ^ The number of contrastive divergence steps
      -> Int -- ^ The number of samples
      -> Natural # x -- ^ Target marginal
      -> Natural # Harmonium f z x -- ^ The harmonium
      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
dualContrastiveDivergence cdn nsmps prr hrm = do
    xs <- sample nsmps prr
    dhrm' <- contrastiveDivergence cdn xs $ transposeHarmonium hrm
    return $ primalIsomorphism . transposeHarmonium $ dualIsomorphism dhrm'

informationProjectionDifferential0
    :: forall f m n r
    . (ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m, Generative Natural n)
    => Int
    -> Natural # Harmonium f m n
    -> Natural # n
    -> Random r (CotangentVector Natural n)
{-# INLINE informationProjectionDifferential0 #-}
informationProjectionDifferential0 nsmps hrm nx = do
    xs <- sample nsmps nx
    return $ harmoniumInformationProjectionDifferential nx xs hrm

informationProjectionCircuit1
    :: ExponentialFamily n
    => Double -- ^ Learning rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Natural # n
    -> Circuit (Natural # n) (CotangentVector Natural n)
    -> Chain (Natural # n)
{-# INLINE informationProjectionCircuit1 #-}
informationProjectionCircuit1 eps gp nx0 dffcrc = accumulateCircuit0 nx0 $ proc ((),nx) -> do
    dnx <- dffcrc -< nx
    gradientCircuit eps gp -< breakPoint $ joinTangentPair nx dnx

-- | Uses SGD to project an exponential family distribution onto the latent
-- distribution of a harmonium.
harmoniumInformationProjection
    :: (ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m, Generative Natural n )
    => Int -- ^ Sample size
    -> Double -- ^ Learning rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Int -- ^ Number of steps
    -> Natural # Harmonium f m n -- ^ Target harmonium (marginal)
    -> Natural # n -- ^ Initial Point
    -> Random s (Natural # n) -- ^ Projected point
{-# INLINE harmoniumInformationProjection #-}
harmoniumInformationProjection nsmps eps gp nstps hrm nx0 = do
    dffcrc <- accumulateRandomFunction0 $ informationProjectionDifferential0 nsmps hrm
    return $ streamChain (informationProjectionCircuit1 eps gp nx0 dffcrc) !! nstps

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

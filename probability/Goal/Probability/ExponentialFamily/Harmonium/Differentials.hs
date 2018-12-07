{-# LANGUAGE
    Arrows,
    RankNTypes,
    TypeOperators,
    DataKinds,
    FlexibleContexts
#-}

-- | A collection of 'Circuit's for computing the differentials in gradient
-- descent algorithms.

module Goal.Probability.ExponentialFamily.Harmonium.Differentials
    ( -- * Differentials
      stochasticRectifiedHarmoniumDifferential
    , harmoniumInformationProjectionDifferential
    , stochasticMixtureModelDifferential
    , contrastiveDivergence
      -- ** Conditional
    , mixtureStochasticConditionalCrossEntropyDifferential
    , conditionalHarmoniumRectificationDifferential
      -- * Pursuit
    , harmoniumInformationProjection
    ) where


--- Imports ---

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Rectification
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conditional


--- Differentials ---


-- | Estimates the stochastic cross entropy differential of a rectified harmonium with
-- respect to the relative entropy, and given an observation.
--
-- NB: Right now I'm using length on a list where I probably don't need to...
stochasticRectifiedHarmoniumDifferential
    :: ( Map Mean Natural f m n, Bilinear f m n, ExponentialFamily m
       , ExponentialFamily n, Generative Natural m, Generative Natural n )
       => Sample m -- ^ Observations
       -> Natural # n -- ^ Rectification Parameters
       -> Natural # Harmonium f m n -- ^ Harmonium
       -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ Differential
{-# INLINE stochasticRectifiedHarmoniumDifferential #-}
stochasticRectifiedHarmoniumDifferential zs rprms hrm = do
    pzxs <- initialPass hrm zs
    qzxs <- sampleRectifiedHarmonium (length zs) (toSingletonSum rprms) hrm
    return $ stochasticCrossEntropyDifferential' pzxs qzxs

-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: (ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m)
    => Natural # Harmonium f m n -- ^ Harmonium
    -> Natural # n -- ^ Model Distribution
    -> Sample n -- ^ Model Samples
    -> CotangentVector Natural n -- ^ Differential Estimate
{-# INLINE harmoniumInformationProjectionDifferential #-}
harmoniumInformationProjectionDifferential hrm px xs =
    let (affmn,nm0) = splitBottomHarmonium hrm
        (nn,nmn) = splitAffine affmn
        nm = fromOneHarmonium nm0
        mxs = sufficientStatistic <$> xs
        mys0 = nmn >$> mxs
        mys = zipWith (\mx my0 -> mx <.> (px <-> nm) - potential (nn <+> my0)) mxs mys0
        ln = fromIntegral $ length xs
        mxht = averagePoint mxs
        myht = sum mys / ln
        foldfun (mx,my) (k,z0) = (k+1,z0 <+> ((my - myht) .> (mx <-> mxht)))
        cvr = uncurry (/>) . foldr foldfun (-1,zero) $ zip mxs mys
     in primalIsomorphism cvr

-- | The stochastic cross entropy differential of a mixture model.
stochasticMixtureModelDifferential
    :: ( Enum e, Legendre Natural o, KnownNat n, ExponentialFamily o )
      => Sample o -- ^ Observations
      -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
      -> CotangentVector Natural (Harmonium Tensor o (Categorical e n)) -- ^ Differential
{-# INLINE stochasticMixtureModelDifferential #-}
stochasticMixtureModelDifferential zs hrm =
    let pxs = empiricalHarmoniumExpectations zs hrm
        qxs = dualTransition hrm
     in primalIsomorphism $ qxs <-> pxs

contrastiveDivergence
    :: ( Generative Natural m, ExponentialFamily m, ExponentialFamily n
       , Generative Natural n, Map Mean Natural f m n, Bilinear f m n )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample m -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium f m n -- ^ The harmonium
      -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- loopM cdn (gibbsPass hrm) xzs0
    return $ stochasticCrossEntropyDifferential' xzs0 xzs1

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates (this is
-- primarily of internal use).
mixtureStochasticConditionalCrossEntropyDifferential
    :: ( Enum e, ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k )
    => Sample x -- ^ Input mean distributions
    -> Sample z -- ^ Output mean distributions
    -> Mean #> Natural # MixtureGLM z e k x -- ^ Function
    -> CotangentVector (Mean #> Natural) (MixtureGLM z e k x) -- ^ Differential
{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential #-}
mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
    -- This could be better optimized but not throwing out the second result of propagate
    let dmglms = dualIsomorphism
            <$> zipWith stochasticMixtureModelDifferential ((:[]) <$> zs) (mglm >$>* xs)
        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
        f = snd $ splitBottomSubLinear mglm
        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
     in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df

-- | A gradient for rectifying gains which won't allow them to be negative.
conditionalHarmoniumRectificationDifferential
    :: ( ExponentialFamily x, Manifold (f z x), Map Mean Natural f z x, Manifold (DeepHarmonium gs (z : zs))
       , Legendre Natural (DeepHarmonium gs (z : zs)) )
    => Double -- ^ Rectification shift
    -> Natural # x -- ^ Rectification parameters
    -> Sample x -- ^ Sample points
    -> Mean #> Natural # f z x -- ^ linear part of ppc
    -> Natural # DeepHarmonium gs (z : zs) -- ^ Gains
    -> CotangentPair Natural (DeepHarmonium gs (z : zs)) -- ^ Rectified PPC
{-# INLINE conditionalHarmoniumRectificationDifferential #-}
conditionalHarmoniumRectificationDifferential rho0 rprms xsmps tns dhrm =
    let lkl = joinBottomSubLinear dhrm tns
        rcts = rectificationCurve rho0 rprms xsmps
        ndhrmlkls = lkl >$>* xsmps
        mdhrmlkls = dualTransition <$> ndhrmlkls
        ptns = potential <$> ndhrmlkls
     in joinTangentPair dhrm . averagePoint
         $ [ primalIsomorphism $ (ptn - rct) .> mdhrmlkl | (rct,mdhrmlkl,ptn) <- zip3 rcts mdhrmlkls ptns ]

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
harmoniumInformationProjection nsmps eps gp nstps hrm nx0 =
    loopChain nstps . chainCircuit nx0 $ proc nx -> do
        xs <- arrM (sample nsmps) -< nx
        dnx <- arr . uncurry $ harmoniumInformationProjectionDifferential hrm -< (nx,xs)
        gradientCircuit eps gp -< breakPoint $ joinTangentPair nx dnx

--dualContrastiveDivergence
--    :: forall s f z x
--    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
--      , Map Mean Natural f x z, Bilinear f z x, Bilinear f x z )
--      => Int -- ^ The number of contrastive divergence steps
--      -> Int -- ^ The number of samples
--      -> Natural # x -- ^ Target marginal
--      -> Natural # Harmonium f z x -- ^ The harmonium
--      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
--dualContrastiveDivergence cdn nsmps prr hrm = do
--    xs <- sample nsmps prr
--    dhrm' <- contrastiveDivergence cdn xs $ transposeHarmonium hrm
--    return $ primalIsomorphism . transposeHarmonium $ dualIsomorphism dhrm'
--
----class FitRectificationParameters (fs :: [* -> * -> *]) (ms :: [*]) where
----    fitRectificationParameters
----        :: Double
----        -> Maybe Int
----        -> Natural # DeepHarmonium fs ms
----        -> Natural # Sum (Tail ms)
----        -> Random s (Natural # Sum (Tail ms))
----
----instance FitRectificationParameters '[] '[m] where
----    {-# INLINE fitRectificationParameters #-}
----    fitRectificationParameters _ _ _ _ = zero
----
----instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
----         , ExponentialFamily n, SampleRectified fs (n : ms), Generative Natural m
----         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
----  => SampleRectified (f : fs) (m : n : ms) where
----    {-# INLINE sampleRectified #-}
----    sampleRectified rprms dhrm = do
----        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
----            (rprm,rprms') = splitSum rprms
----        (ys,xs) <- fmap hUnzip . sampleRectified rprms' $ biasBottom rprm dhrm'
----        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
----        return . hZip zs $ hZip ys xs
----
----

{-# LANGUAGE Arrows #-}
-- | A collection of algorithms for optimizing harmoniums.

module Goal.Graphical.Learning
    ( -- * Expectation Maximization
      expectationMaximization
    , expectationMaximizationAscent
    , gibbsExpectationMaximization
    -- * Differentials
    , harmoniumInformationProjectionDifferential
    , contrastiveDivergence
    , conditionalHarmoniumConjugationDifferential
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium
import Goal.Graphical.Hybrid
import Goal.Graphical.Inference


--- Differentials ---


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: ( Map Natural f z x, LegendreExponentialFamily z
       , ExponentialFamily x, Generative Natural x)
    => Int
    -> Natural # Harmonium z f x -- ^ Harmonium
    -> Natural # x -- ^ Model Distribution
    -> Random r (Mean # x) -- ^ Differential Estimate
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
    let (affmn,nm0) = splitBottomHarmonium hrm
        (nn,nmn) = splitAffine affmn
        nm = fromOneHarmonium nm0
        mxs = sufficientStatistic <$> xs
        mys0 = nmn >$> mxs
        mys = zipWith (\mx my0 -> mx <.> (px - nm) - potential (nn + my0)) mxs mys0
        ln = fromIntegral $ length xs
        mxht = average mxs
        myht = sum mys / ln
        foldfun (mx,my) (k,z0) = (k+1,z0 + ((my - myht) .> (mx - mxht)))
    return . uncurry (/>) . foldr foldfun (-1,0) $ zip mxs mys

-- | Contrastive divergence on harmoniums (<https://www.mitpressjournals.org/doi/abs/10.1162/089976602760128018?casa_token=x_Twj1HaXcMAAAAA:7-Oq181aubCFwpG-f8Lo1wRKvGnmujzl8zjn9XbeO5nGhfvKCCQjsu4K4pJCkMNYUYWqc2qG7TRXBg Hinton, 2019>).
contrastiveDivergence
    :: ( Generative Natural z, ExponentialFamily z, Generative Natural x
       , ExponentialFamily x, Bilinear f z x, Map Natural f x z
       , Map Natural f z x )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample z -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium z f x -- ^ The harmonium
      -> Random s (Mean # Harmonium z f x) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticRelativeEntropyDifferential xzs0 xzs1

-- | An approximate differntial for conjugating a harmonium likelihood.
conditionalHarmoniumConjugationDifferential
    :: ( Propagate Natural f y z, Manifold (g y x)
       , LegendreExponentialFamily (Harmonium y g x)
       , LegendreExponentialFamily x, ExponentialFamily y, ExponentialFamily z )
    => Double -- ^ Conjugation shift
    -> Natural # z -- ^ Conjugation parameters
    -> Sample z -- ^ Sample points
    -> Natural # ConditionalHarmonium f y g x z
    -> Mean # ConditionalHarmonium f y g x z
conditionalHarmoniumConjugationDifferential rho0 rprms xsmps chrm =
    let rcts = conjugationCurve rho0 rprms xsmps
        mhrms = transition <$> nhrms
        ptns = potential <$> nhrms
        dhrms = [ (ptn - rct) .> mhrm | (rct,mhrm,ptn) <- zip3 rcts mhrms ptns ]
        (dchrm,nhrms) = propagate dhrms (sufficientStatistic <$> xsmps) chrm
     in dchrm


--- Expectation Maximization ---


-- | EM for latent aviarlbe models.
expectationMaximization
    :: ( Transition Mean d x, ExpectationMaximization c x, LatentVariable x o l )
    => [o]
    -> c # x
    -> d # x
expectationMaximization zs hrm = transition $ expectationStep zs hrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
expectationMaximizationAscent
    :: ( LegendreExponentialFamily x, ExpectationMaximization Natural x
       , LatentVariable x o l )
    => Double -> GradientPursuit -> [o] -> (Natural # x) -> [Natural # x]
expectationMaximizationAscent eps gp zs nhrm =
    let mhrm' = expectationStep zs nhrm
     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
gibbsExpectationMaximization
    :: ( Generative Natural z, Generative Natural x, LegendreExponentialFamily x
       , Manifold (Harmonium z f x), Map Natural f x z
       , ExponentialFamily z, Bilinear f z x, Map Natural f z x )
    => Double
    -> Int
    -> Int
    -> GradientPursuit
    -> Sample z -- ^ Observations
    -> Natural # Harmonium z f x -- ^ Current Harmonium
    -> Chain (Random r) (Natural # Harmonium z f x) -- ^ Harmonium Chain
gibbsExpectationMaximization eps cdn nbtch gp zs0 nhrm0 =
    let mhrm0 = expectationStep zs0 nhrm0
     in chainCircuit nhrm0 $ proc nhrm -> do
         zs <- minibatcher nbtch zs0 -< ()
         xzs0 <- arrM (uncurry initialPass) -< (nhrm,zs)
         xzs1 <- arrM (\(x,y) -> iterateM' cdn (gibbsPass x) y) -< (nhrm,xzs0)
         let dff = mhrm0 - averageSufficientStatistic xzs1
         gradientCircuit eps gp -< (nhrm,vanillaGradient dff)


---- | Estimates the stochastic cross entropy differential of a conjugated harmonium with
---- respect to the relative entropy, and given an observation.
--stochasticConjugatedHarmoniumDifferential
--    :: ( Map Natural f z x, Bilinear f z x, ExponentialFamily z
--       , ExponentialFamily x, Generative Natural z, Generative Natural x )
--       => Sample z -- ^ Observations
--       -> Natural # x -- ^ Conjugation Parameters
--       -> Natural # Harmonium f z x -- ^ Harmonium
--       -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ Differential
--stochasticConjugatedHarmoniumDifferential zs rprms hrm = do
--    pzxs <- initialPass hrm zs
--    qzxs <- sampleConjugatedHarmonium (length zs) (toSingletonSum rprms) hrm
--    return $ stochasticCrossEntropyDifferential' pzxs qzxs
--
---- | The stochastic conditional cross-entropy differential, based on target
---- inputs and outputs expressed as distributions in mean coordinates.
--mixtureStochasticConditionalCrossEntropyDifferential
--    :: ( ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k )
--    => Sample x -- ^ Input mean distributions
--    -> Sample z -- ^ Output mean distributions
--    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
--    -> CotangentVector (Mean #> Natural) (MixtureGLM z k x) -- ^ Differential
--mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
--    -- This could be better optimized but not throwing out the second result of propagate
--    let dmglms = dualIsomorphism
--            <$> zipWith stochasticMixtureDifferential ((:[]) <$> zs) (mglm >$>* xs)
--        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
--        f = snd $ splitBottomSubLinear mglm
--        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
--     in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df
--
--
--
----dualContrastiveDivergence
----    :: forall s f z x
----    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
----      , Map Natural f x z, Bilinear f z x, Bilinear f x z )
----      => Int -- ^ The number of contrastive divergence steps
----      -> Int -- ^ The number of samples
----      -> Natural # x -- ^ Target marginal
----      -> Natural # Harmonium f z x -- ^ The harmonium
----      -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ The gradient estimate
----dualContrastiveDivergence cdn nsmps prr hrm = do
----    xs <- sample nsmps prr
----    dhrm' <- contrastiveDivergence cdn xs $ transposeHarmonium hrm
----    return $ primalIsomorphism . transposeHarmonium $ dualIsomorphism dhrm'
----
------class FitConjugationParameters (fs :: [* -> * -> *]) (ms :: [*]) where
------    fitConjugationParameters
------        :: Double
------        -> Maybe Int
------        -> Natural # DeepHarmonium fs ms
------        -> Natural # Sum (Tail ms)
------        -> Random s (Natural # Sum (Tail ms))
------
------instance FitConjugationParameters '[] '[m] where
------    fitConjugationParameters _ _ _ _ = zero
------
------instance ( Manifold (DeepHarmonium fs (n : ms)), Map Natural f z x, Manifold (Sum ms)
------         , ExponentialFamily n, SampleConjugated fs (n : ms), Generative Natural m
------         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
------  => SampleConjugated (f : fs) (m : n : ms) where
------    sampleConjugated rprms dhrm = do
------        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
------            (rprm,rprms') = splitSum rprms
------        (ys,xs) <- fmap hUnzip . sampleConjugated rprms' $ biasBottom rprm dhrm'
------        zs <- samplePoint $ mapReplicatedPoint (pn +) (pf >$>* ys)
------        return . hZip zs $ hZip ys xs
------
------
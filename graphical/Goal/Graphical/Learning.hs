{-# LANGUAGE Arrows #-}
-- | A collection of algorithms for optimizing harmoniums.

module Goal.Graphical.Learning where
    --( -- * Expectation Maximization
    --  expectationMaximization
    --, expectationMaximizationAscent
    --, gibbsExpectationMaximization
    --, latentProcessExpectationMaximization
    --, latentProcessExpectationMaximizationAscent
    ---- * Differentials
    --, harmoniumInformationProjectionDifferential
    --, contrastiveDivergence
    --) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium
import Goal.Graphical.Generative.Dynamic
import Goal.Graphical.Inference


--- Differentials ---


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: ( Map Natural f y x, LegendreExponentialFamily z
       , SamplePoint w ~ SamplePoint x, Translation z y
       , ExponentialFamily x, ExponentialFamily w, Generative Natural w )
    => Int
    -> Natural # Harmonium f y x z w -- ^ Harmonium
    -> Natural # w -- ^ Model Distribution
    -> Random r (Mean # w) -- ^ Differential Estimate
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
    let (lkl,nw) = split hrm
        mys0 = lkl >$>* xs
        mws = sufficientStatistic <$> xs
        mys = zipWith (\mw my0 -> mw <.> (px - nw) - potential my0) mws mys0
        ln = fromIntegral $ length xs
        mwht = average mws
        myht = sum mys / ln
        foldfun (mw,my) (k,z0) = (k+1,z0 + ((my - myht) .> (mw - mwht)))
    return . uncurry (/>) . foldr foldfun (-1,0) $ zip mws mys

-- | Contrastive divergence on harmoniums (<https://www.mitpressjournals.org/doi/abs/10.1162/089976602760128018?casa_token=x_Twj1HaXcMAAAAA:7-Oq181aubCFwpG-f8Lo1wRKvGnmujzl8zjn9XbeO5nGhfvKCCQjsu4K4pJCkMNYUYWqc2qG7TRXBg Hinton, 2019>).
--contrastiveDivergence
--    :: ( Generative Natural z, ExponentialFamily z, Generative Natural x
--       , ExponentialFamily x, Bilinear f z x, Map Natural f x z
--       , Map Natural f z x )
--      => Int -- ^ The number of contrastive divergence steps
--      -> Sample z -- ^ The initial states of the Gibbs chains
--      -> Natural # Harmonium f y x z w -- ^ The harmonium
--      -> Random s (Mean # Harmonium f y x z w) -- ^ The gradient estimate
--contrastiveDivergence cdn zs hrm = do
--    xzs0 <- initialPass hrm zs
--    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
--    return $ stochasticRelativeEntropyDifferential xzs0 xzs1


--- Expectation Maximization ---


-- | EM for latent variable odels.
expectationMaximization
    :: ( DuallyFlatExponentialFamily f, ExpectationMaximization f )
    => Observations f
    -> Natural # f
    -> Natural # f
expectationMaximization zs hrm = transition $ expectationStep zs hrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
expectationMaximizationAscent
    :: (LegendreExponentialFamily f, ExpectationMaximization f)
    => Double -> GradientPursuit -> Observations f -> Natural # f -> [Natural # f]
expectationMaximizationAscent eps gp zs nhrm =
    let mhrm' = expectationStep zs nhrm
     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm

-- | Ascent of the EM objective on harmoniums for when the expectation
-- step can't be computed in closed-form. The convergent harmonium distribution
-- of the output harmonium-list is the result of 1 iteration of the EM
-- algorithm.
--gibbsExpectationMaximization
--    :: ( Generative Natural z, Generative Natural x, LegendreExponentialFamily x
--       , Manifold (Harmonium f z x), Map Natural f x z
--       , ExponentialFamily z, Bilinear f z x, Map Natural f z x )
--    => Double
--    -> Int
--    -> Int
--    -> GradientPursuit
--    -> Sample z -- ^ Observations
--    -> Natural # Harmonium f z x -- ^ Current Harmonium
--    -> Chain (Random r) (Natural # Harmonium f z x) -- ^ Harmonium Chain
--gibbsExpectationMaximization eps cdn nbtch gp zs0 nhrm0 =
--    let mhrm0 = expectationStep zs0 nhrm0
--     in chainCircuit nhrm0 $ proc nhrm -> do
--         zs <- minibatcher nbtch zs0 -< ()
--         xzs0 <- arrM (uncurry initialPass) -< (nhrm,zs)
--         xzs1 <- arrM (\(x,y) -> iterateM' cdn (gibbsPass x) y) -< (nhrm,xzs0)
--         let dff = mhrm0 - averageSufficientStatistic xzs1
--         gradientCircuit eps gp -< (nhrm,vanillaGradient dff)

latentProcessExpectationStep
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , Transition Natural Mean w, Transition Natural Mean (Harmonium g x x w w)
       , Manifold (Harmonium g x x w w)
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Observations (LatentProcess f g y x z w)
    -> Natural # LatentProcess f g y x z w
    -> (Mean # w, Mean # Harmonium f y x z w, Mean # Harmonium g x x w w)
latentProcessExpectationStep zss ltnt =
    let (prr,emsn,trns) = splitLatentProcess ltnt
        (smthss,hrmss) = unzip $ conjugatedSmoothing trns emsn prr <$> zss
        mprr = average $ toMean . head <$> smthss
        mtrns = average $ toMean <$> concat hrmss
        mws = toMean <$> concat smthss
        mzs = sufficientStatistic <$> concat zss
        mys = anchor <$> mzs
        mxs = anchor <$> mws
        memsn = joinHarmonium (average mzs) (mys >$< mxs) (average mws)
     in (mprr,memsn,mtrns)

latentProcessExpectationMaximization
    :: ( ConjugatedLikelihood g x x w w, ConjugatedLikelihood f y x z w
       , Transition Natural Mean w, Transition Natural Mean (Harmonium g x x w w)
       , Transition Mean Natural w
       , Transition Mean Natural (Harmonium f y x z w)
       , Transition Mean Natural (Harmonium g x x w w)
       , Manifold (Harmonium g x x w w)
       , Bilinear g x x, Map Natural f x y, Bilinear f y x
       , SamplePoint y ~ SamplePoint z )
    => Observations (LatentProcess f g y x z w)
    -> Natural # LatentProcess f g y x z w
    -> Natural # LatentProcess f g y x z w
latentProcessExpectationMaximization zss ltnt =
    let (mprr,memsn,mtrns) = latentProcessExpectationStep zss ltnt
        prr' = toNatural mprr
        emsn' = fst . split $ toNatural memsn
        trns' = fst . split $ toNatural mtrns
     in joinLatentProcess prr' emsn' trns'

--latentProcessExpectationMaximizationAscent
--    :: ( ConjugatedLikelihood g x x, ConjugatedLikelihood f z x
--       , Propagate Natural g x x, Propagate Natural f z x
--       , ExponentialFamily z, Bilinear g x x, DuallyFlatExponentialFamily x
--       , Bilinear f z x, Map Natural f x z
--       , LegendreExponentialFamily (Harmonium g x x)
--       , LegendreExponentialFamily (Harmonium f z x) )
--    => Double
--    -> Int
--    -> GradientPursuit
--    -> [Sample z]
--    -> Natural # LatentProcess f g z x
--    -> Natural # LatentProcess f g z x
--latentProcessExpectationMaximizationAscent eps nstps gp zss ltnt =
--    let (mprr,mehrm,mthrm) = latentProcessExpectationStep zss ltnt
--        (nprr,nemsn,ntrns) = splitLatentProcess ltnt
--        nsmth = toNatural . snd $ splitBottomHarmonium mehrm
--        nehrm = joinConjugatedHarmonium nemsn nsmth
--        nthrm = joinConjugatedHarmonium ntrns nsmth
--        nprr' = (!! nstps)
--            $ vanillaGradientSequence (relativeEntropyDifferential mprr) (-eps) gp nprr
--        nemsn' = fst . splitBottomHarmonium . (!! nstps)
--            $ vanillaGradientSequence (relativeEntropyDifferential mehrm) (-eps) gp nehrm
--        ntrns' = fst . splitBottomHarmonium . (!! nstps)
--            $ vanillaGradientSequence (relativeEntropyDifferential mthrm) (-eps) gp nthrm
--     in joinLatentProcess nprr' nemsn' ntrns'

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
-- -- | An approximate differntial for conjugating a harmonium likelihood.
-- conditionalHarmoniumConjugationDifferential
--     :: ( Propagate Natural g z y, Manifold (g z y)
--        , LegendreExponentialFamily (Harmonium g y x)
--        , LegendreExponentialFamily x, ExponentialFamily y, ExponentialFamily z )
--     => Double -- ^ Conjugation shift
--     -> Natural # z -- ^ Conjugation parameters
--     -> Sample z -- ^ Sample points
--     -> Natural # ConditionalHarmonium g f z x y
--     -> Mean # ConditionalHarmonium g f z x y
-- conditionalHarmoniumConjugationDifferential rho0 rprms xsmps chrm =
--     let rcts = conjugationCurve rho0 rprms xsmps
--         mhrms = transition <$> nhrms
--         ptns = potential <$> nhrms
--         dhrms = [ (ptn - rct) .> mhrm | (rct,mhrm,ptn) <- zip3 rcts mhrms ptns ]
--         (dchrm,nhrms) = propagate dhrms (sufficientStatistic <$> xsmps) chrm
--      in dchrm

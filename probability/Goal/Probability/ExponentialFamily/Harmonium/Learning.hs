{-# LANGUAGE Arrows #-}
-- | A collection of algorithms for optimizing harmoniums.

module Goal.Probability.ExponentialFamily.Harmonium.Learning
    ( -- * Differentials
      harmoniumInformationProjectionDifferential
    , contrastiveDivergence
      -- ** Conditional
    , conditionalExpectationMaximizationAscent
    , conditionalHarmoniumConjugationDifferential
    -- * Expectation Maximization
    , expectationMaximization
    , expectationMaximizationAscent
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conditional
import Goal.Probability.ExponentialFamily.Harmonium.Inference

import qualified Data.Vector as V
import System.Random.MWC.Probability hiding (initialize,sample)
import System.Random.MWC.Distributions (uniformShuffle)



--- Differentials ---


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: ( Map Mean Natural f z x, LegendreExponentialFamily z
       , ExponentialFamily x, Generative Natural x)
    => Int
    -> Natural # Harmonium z f x -- ^ Harmonium
    -> Natural # x -- ^ Model Distribution
    -> Random r (Mean # x) -- ^ Differential Estimate
{-# INLINE harmoniumInformationProjectionDifferential #-}
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
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
    return . uncurry (/>) . foldr foldfun (-1,zero) $ zip mxs mys

contrastiveDivergence
    :: ( Generative Natural z, ExponentialFamily z, Generative Natural x
       , ExponentialFamily x, Bilinear f z x, Map Mean Natural f x z, Map Mean Natural f z x )
      => Int -- ^ The number of contrastive divergence steps
      -> Sample z -- ^ The initial states of the Gibbs chains
      -> Natural # Harmonium z f x -- ^ The harmonium
      -> Random s (Mean # Harmonium z f x) -- ^ The gradient estimate
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticCrossEntropyDifferential xzs0 xzs1


--- Expectation Maximization ---


-- | EM implementation for mixture models/categorical harmoniums.
expectationMaximization
    :: ( DuallyFlatExponentialFamily (Harmonium z f x), LegendreExponentialFamily x
       , ExponentialFamily z, Bilinear f z x, Map Mean Natural f x z )
    => Sample z -- ^ Observations
    -> Natural # Harmonium z f x -- ^ Current Harmonium
    -> Natural # Harmonium z f x -- ^ Updated Harmonium
{-# INLINE expectationMaximization #-}
expectationMaximization zs hrm = transition $ harmoniumEmpiricalExpectations zs hrm

expectationMaximizationAscent
    :: ( LegendreExponentialFamily (Harmonium z f x), LegendreExponentialFamily x
       , ExponentialFamily z, Bilinear f z x, Map Mean Natural f x z )
    => Double
    -> GradientPursuit
    -> Sample z -- ^ Observations
    -> Natural # Harmonium z f x -- ^ Current Harmonium
    -> [Natural # Harmonium z f x] -- ^ Updated Harmonium
{-# INLINE expectationMaximizationAscent #-}
expectationMaximizationAscent eps gp zs nhrm =
    let mhrm' = harmoniumEmpiricalExpectations zs nhrm
     in vanillaGradientSequence (crossEntropyDifferential mhrm') (-eps) gp nhrm

-- | NB: Write now this never shuffles the training data, as it probably should.
conditionalExpectationMaximizationAscent
    :: ( ExponentialFamily z, LegendreExponentialFamily y, KnownNat k )
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient Pursuit Algorithm
    -> Int -- ^ Batch size
    -> Int -- ^ Number of iterations
    -> Sample (y,z) -- ^ Conditioning sample
    -> Natural #> ConditionalMixture y k z
    -> Random r (Natural #> ConditionalMixture y k z)
{-# INLINE conditionalExpectationMaximizationAscent #-}
conditionalExpectationMaximizationAscent eps gp nbtch nstps yzs0 chrm0 = do
    let chrmcrc = loopCircuit' chrm0 $ proc (mhrmzs,chrm) -> do
            let (mhrms,zs) = unzip mhrmzs
            let dhrms = zipWith (<->) mhrms $ transition <$> hrmhts
                (dchrm,hrmhts) = propagate dhrms zs chrm
            gradientCircuit (-eps) gp -< (chrm,vanillaGradient dchrm)
    let (ys0,zs0) = unzip yzs0
        mhrms0 = conditionalHarmoniumEmpiricalExpectations ys0 chrm0
        ncycs = 1 + div (length yzs0 - 1) (nstps * nbtch)
    mhrmzs0 <- replicateM ncycs (shuffleList . zip mhrms0 $ sufficientStatistic <$> zs0)
    let mhrmzss = take nstps . breakEvery nbtch $ concat mhrmzs0
    iterateCircuit chrmcrc mhrmzss

-- | A gradient for conjugateing gains which won't allow them to be negative.
conditionalHarmoniumConjugationDifferential
    :: ( ExponentialFamily x, Map Mean Natural f z x
       , LegendreExponentialFamily (DeepHarmonium z fxs) )
    => Double -- ^ Conjugation shift
    -> Natural # x -- ^ Conjugation parameters
    -> Sample x -- ^ Sample points
    -> Natural #> f z x -- ^ linear part of ppc
    -> Natural # DeepHarmonium z fxs -- ^ Gains
    -> Mean # DeepHarmonium z fxs -- ^ Conjugated PPC
{-# INLINE conditionalHarmoniumConjugationDifferential #-}
conditionalHarmoniumConjugationDifferential rho0 rprms xsmps tns dhrm =
    let lkl = joinConditionalDeepHarmonium dhrm tns
        rcts = conjugationCurve rho0 rprms xsmps
        ndhrmlkls = lkl >$>* xsmps
        mdhrmlkls = transition <$> ndhrmlkls
        ptns = potential <$> ndhrmlkls
     in averagePoint [ (ptn - rct) .> mdhrmlkl | (rct,mdhrmlkl,ptn) <- zip3 rcts mdhrmlkls ptns ]

---- | Estimates the stochastic cross entropy differential of a conjugated harmonium with
---- respect to the relative entropy, and given an observation.
--stochasticConjugatedHarmoniumDifferential
--    :: ( Map Mean Natural f z x, Bilinear f z x, ExponentialFamily z
--       , ExponentialFamily x, Generative Natural z, Generative Natural x )
--       => Sample z -- ^ Observations
--       -> Natural # x -- ^ Conjugation Parameters
--       -> Natural # Harmonium f z x -- ^ Harmonium
--       -> Random s (CotangentVector Natural (Harmonium f z x)) -- ^ Differential
--{-# INLINE stochasticConjugatedHarmoniumDifferential #-}
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
--{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential #-}
--mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
--    -- This could be better optimized but not throwing out the second result of propagate
--    let dmglms = dualIsomorphism
--            <$> zipWith stochasticMixtureDifferential ((:[]) <$> zs) (mglm >$>* xs)
--        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
--        f = snd $ splitBottomSubLinear mglm
--        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
--     in primalIsomorphism $ joinBottomSubLinear (averagePoint dmglms) df
--
shuffleList :: [a] -> Random r [a]
shuffleList xs = fmap V.toList . Prob $ uniformShuffle (V.fromList xs)
--
--
----dualContrastiveDivergence
----    :: forall s f z x
----    . ( Generative Natural z, ExponentialFamily z, ExponentialFamily x, Generative Natural x
----      , Map Mean Natural f x z, Bilinear f z x, Bilinear f x z )
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
------    {-# INLINE fitConjugationParameters #-}
------    fitConjugationParameters _ _ _ _ = zero
------
------instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f z x, Manifold (Sum ms)
------         , ExponentialFamily n, SampleConjugated fs (n : ms), Generative Natural m
------         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
------  => SampleConjugated (f : fs) (m : n : ms) where
------    {-# INLINE sampleConjugated #-}
------    sampleConjugated rprms dhrm = do
------        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
------            (rprm,rprms') = splitSum rprms
------        (ys,xs) <- fmap hUnzip . sampleConjugated rprms' $ biasBottom rprm dhrm'
------        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
------        return . hZip zs $ hZip ys xs
------
------

{-# LANGUAGE Arrows #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

-- | A collection of algorithms for optimizing harmoniums.
module Goal.Graphical.Learning (
    -- * Expectation Maximization
    expectationMaximization,
    expectationMaximizationAscent,
    stochasticConjugatedEMAscent,
    gibbsExpectationMaximization,
    latentProcessExpectationMaximization,
    latentProcessExpectationMaximizationAscent,

    -- * Maximum Likelihood
    stochasticConjugatedMLAscent,

    -- * Differentials
    harmoniumInformationProjectionDifferential,
    contrastiveDivergence,
) where

--- Imports ---

-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models
import Goal.Graphical.Models.Dynamic
import Goal.Graphical.Models.Harmonium

--- Differentials ---

{- | The differential of the dual relative entropy. Minimizing this results in
the information projection of the model against the marginal distribution of
the given harmonium. This is more efficient than the generic version.
-}
harmoniumInformationProjectionDifferential ::
    ( LegendreExponentialFamily z
    , SamplePoint w ~ SamplePoint x
    , KnownLinear f y x
    , Translation z y
    , ExponentialFamily x
    , ExponentialFamily w
    , Generative Natural w
    ) =>
    Int ->
    -- | Harmonium
    Natural # AffineHarmonium f y x z w ->
    -- | Model Distribution
    Natural # w ->
    -- | Differential Estimate
    Random (Mean # w)
harmoniumInformationProjectionDifferential n hrm px = do
    xs <- sample n px
    let (lkl, nw) = split hrm
        mys0 = lkl >$>* xs
        mws = sufficientStatistic <$> xs
        mys = zipWith (\mw my0 -> mw <.> (px - nw) - potential my0) mws mys0
        ln = fromIntegral $ length xs
        mwht = average mws
        myht = sum mys / ln
        foldfun (mw, my) (k, z0) = (k + 1, z0 + ((my - myht) .> (mw - mwht)))
    return . uncurry (/>) . foldr foldfun (-1, 0) $ zip mws mys

-- | Contrastive divergence on harmoniums (<https://www.mitpressjournals.org/doi/abs/10.1162/089976602760128018?casa_token=x_Twj1HaXcMAAAAA:7-Oq181aubCFwpG-f8Lo1wRKvGnmujzl8zjn9XbeO5nGhfvKCCQjsu4K4pJCkMNYUYWqc2qG7TRXBg Hinton, 2019>).
contrastiveDivergence ::
    ( Generative Natural z
    , ExponentialFamily z
    , Translation w x
    , Generative Natural w
    , ExponentialFamily y
    , Translation z y
    , KnownLinear f x y
    , KnownLinear f y x
    , LegendreExponentialFamily w
    , SamplePoint y ~ SamplePoint z
    , SamplePoint x ~ SamplePoint w
    , ExponentialFamily x
    ) =>
    -- | The number of contrastive divergence steps
    Int ->
    -- | The initial states of the Gibbs chains
    Sample z ->
    -- | The harmonium
    Natural # AffineHarmonium f y x z w ->
    -- | The gradient estimate
    Random (Mean # AffineHarmonium f y x z w)
contrastiveDivergence cdn zs hrm = do
    xzs0 <- initialPass hrm zs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticRelativeEntropyDifferential xzs0 xzs1

--- Expectation Maximization ---

-- | A single iteration of EM for 'Harmonium' based models.
expectationMaximization ::
    ( ExponentialFamily x
    , LegendreExponentialFamily z
    , KnownLinear f z0 x0
    , KnownLinear f x0 z0
    , Translation x x0
    , Translation z z0
    , DuallyFlatExponentialFamily (AffineHarmonium f x0 z0 x z)
    ) =>
    Sample x ->
    Natural # AffineHarmonium f x0 z0 x z ->
    Natural # AffineHarmonium f x0 z0 x z
expectationMaximization zs hrm = transition $ expectationStep zs hrm

{- | Ascent of the EM objective on harmoniums for when the expectation
step can't be computed in closed-form. The convergent harmonium distribution
of the output harmonium-list is the result of 1 iteration of the EM
algorithm.
-}
expectationMaximizationAscent ::
    ( LegendreExponentialFamily (AffineHarmonium f y x z w)
    , ExponentialFamily z
    , KnownLinear f x y
    , KnownLinear f y x
    , Translation z y
    , Translation w x
    , LegendreExponentialFamily w
    ) =>
    Double ->
    GradientPursuit ->
    Sample z ->
    Natural # AffineHarmonium f y x z w ->
    [Natural # AffineHarmonium f y x z w]
expectationMaximizationAscent eps gp zs nhrm =
    let mhrm' = expectationStep zs nhrm
     in vanillaGradientSequence (relativeEntropyDifferential mhrm') (-eps) gp nhrm

stochasticConjugatedMLAscent ::
    ( Generative Natural z
    , LegendreExponentialFamily z
    , Generative Natural x
    , ConjugatedLikelihood f x0 z0 x z
    ) =>
    Double ->
    GradientPursuit ->
    Sample x ->
    Int ->
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Harmonium Chain
    Chain Random (Natural # AffineHarmonium f x0 z0 x z)
stochasticConjugatedMLAscent eps gp xs0 nbtch nhrm0 = chainCircuit nhrm0 $ proc nhrm -> do
    xs <- minibatcher nbtch xs0 -< ()
    xzs <- arrM (sample nbtch) -< nhrm
    let mhrm' = expectationStep xs nhrm
    let dff = mhrm' - averageSufficientStatistic xzs
    gradientCircuit eps gp -< (nhrm, vanillaGradient dff)

{- | Ascent of the EM objective on harmoniums for when the expectation
step can't be computed in closed-form. The convergent harmonium distribution
of the output harmonium-list is the result of 1 iteration of the EM
algorithm.
-}
stochasticConjugatedEMAscent ::
    ( Generative Natural z
    , LegendreExponentialFamily z
    , Generative Natural x
    , ConjugatedLikelihood f x0 z0 x z
    ) =>
    Double ->
    GradientPursuit ->
    Sample x ->
    Int ->
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Harmonium Chain
    Chain Random (Natural # AffineHarmonium f x0 z0 x z)
stochasticConjugatedEMAscent eps gp xs0 nbtch nhrm0 = chainCircuit nhrm0 $ proc nhrm -> do
    xs <- minibatcher nbtch xs0 -< ()
    xzs <- arrM (sample nbtch) -< nhrm
    let mhrm' = expectationStep xs nhrm0
    let dff = mhrm' - averageSufficientStatistic xzs
    gradientCircuit eps gp -< (nhrm, vanillaGradient dff)

{- | Ascent of the EM objective on harmoniums for when the expectation
step can't be computed in closed-form. The convergent harmonium distribution
of the output harmonium-list is the result of 1 iteration of the EM
algorithm.
-}
gibbsExpectationMaximization ::
    ( ExponentialFamily z
    , Manifold w
    , Translation z y
    , Translation w x
    , SamplePoint y ~ SamplePoint z
    , SamplePoint w ~ SamplePoint x
    , KnownLinear f x y
    , KnownLinear f y x
    , ExponentialFamily y
    , Generative Natural w
    , ExponentialFamily x
    , Generative Natural z
    , Manifold (AffineHarmonium f y x z w)
    , LegendreExponentialFamily w
    ) =>
    Double ->
    Int ->
    Int ->
    GradientPursuit ->
    -- | Observations
    Sample z ->
    -- | Current Harmonium
    Natural # AffineHarmonium f y x z w ->
    -- | Harmonium Chain
    Chain Random (Natural # AffineHarmonium f y x z w)
gibbsExpectationMaximization eps cdn nbtch gp zs0 nhrm0 =
    let mhrm0 = expectationStep zs0 nhrm0
     in chainCircuit nhrm0 $ proc nhrm -> do
            zs <- minibatcher nbtch zs0 -< ()
            xzs0 <- arrM (uncurry initialPass) -< (nhrm, zs)
            xzs1 <- arrM (\(x, y) -> iterateM' cdn (gibbsPass x) y) -< (nhrm, xzs0)
            let dff = mhrm0 - averageSufficientStatistic xzs1
            gradientCircuit eps gp -< (nhrm, vanillaGradient dff)

latentProcessExpectationStep ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    , Transition Natural Mean w
    , Transition Natural Mean (AffineHarmonium g x x w w)
    , Manifold (AffineHarmonium g x x w w)
    , SamplePoint y ~ SamplePoint z
    ) =>
    Observations (LatentProcess f g y x z w) ->
    Natural # LatentProcess f g y x z w ->
    (Mean # w, Mean # AffineHarmonium f y x z w, Mean # AffineHarmonium g x x w w)
latentProcessExpectationStep zss ltnt =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        (smthss, hrmss) = unzip $ conjugatedSmoothing0 prr emsn trns <$> zss
        mprr = average $ toMean . head <$> smthss
        mtrns = average $ toMean <$> concat hrmss
        mws = toMean <$> concat smthss
        mzs = sufficientStatistic <$> concat zss
        mys = anchor <$> mzs
        mxs = anchor <$> mws
        memsn = joinHarmonium (average mzs) (mys >$< mxs) (average mws)
     in (mprr, memsn, mtrns)

-- | Direct expectation maximization for 'LatentProcess'es.
latentProcessExpectationMaximization ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    , Transition Natural Mean w
    , Transition Natural Mean (AffineHarmonium g x x w w)
    , Transition Mean Natural w
    , Transition Mean Natural (AffineHarmonium f y x z w)
    , Transition Mean Natural (AffineHarmonium g x x w w)
    , Manifold (AffineHarmonium g x x w w)
    , SamplePoint y ~ SamplePoint z
    ) =>
    Observations (LatentProcess f g y x z w) ->
    Natural # LatentProcess f g y x z w ->
    Natural # LatentProcess f g y x z w
latentProcessExpectationMaximization zss ltnt =
    let (mprr, memsn, mtrns) = latentProcessExpectationStep zss ltnt
        prr' = toNatural mprr
        emsn' = fst . split $ toNatural memsn
        trns' = fst . split $ toNatural mtrns
     in joinLatentProcess prr' emsn' trns'

{- | Expectation maximization for 'LatentProcess'es approximated through
gradient ascent.
-}
latentProcessExpectationMaximizationAscent ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    , DuallyFlatExponentialFamily w
    , LegendreExponentialFamily (AffineHarmonium f y x z w)
    , LegendreExponentialFamily (AffineHarmonium g x x w w)
    , SamplePoint y ~ SamplePoint z
    ) =>
    Double ->
    Int ->
    GradientPursuit ->
    [Sample z] ->
    Natural # LatentProcess f g y x z w ->
    Natural # LatentProcess f g y x z w
latentProcessExpectationMaximizationAscent eps nstps gp zss ltnt =
    let (mprr, mehrm, mthrm) = latentProcessExpectationStep zss ltnt
        (nprr, nemsn, ntrns) = splitLatentProcess ltnt
        neql0 = toNatural . snd $ split mehrm
        neql1 = toNatural . snd $ split mthrm
        nehrm = joinConjugatedHarmonium nemsn neql0
        nthrm = joinConjugatedHarmonium ntrns neql1
        nprr' =
            (!! nstps)
                $ vanillaGradientSequence (relativeEntropyDifferential mprr) (-eps) gp nprr
        nemsn' =
            fst
                . split
                . (!! nstps)
                $ vanillaGradientSequence (relativeEntropyDifferential mehrm) (-eps) gp nehrm
        ntrns' =
            fst
                . split
                . (!! nstps)
                $ vanillaGradientSequence (relativeEntropyDifferential mthrm) (-eps) gp nthrm
     in joinLatentProcess nprr' nemsn' ntrns'

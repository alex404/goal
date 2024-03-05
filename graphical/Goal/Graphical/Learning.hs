{-# LANGUAGE Arrows #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

-- | A collection of algorithms for optimizing harmoniums.
module Goal.Graphical.Learning (
    -- * Expectation Maximization
    expectationMaximization,
    expectationMaximizationAscent,
    stochasticConjugatedEMAscent,
    gibbsExpectationMaximization,

    -- ** Dynamic
    latentProcessExpectationMaximization,
    latentProcessExpectationMaximizationAscent,

    -- * Maximum Likelihood
    stochasticConjugatedMLAscent,

    -- * Differentials
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

--- I need to fix this. It's cool but I've lost the plot
-}

-- harmoniumInformationProjectionDifferential ::
--     ( LegendreExponentialFamily z
--     , SamplePoint x ~ SamplePoint x0
--     , KnownLinear f z0 x0
--     , LinearSubspace z z0
--     , ExponentialFamily x
--     , ExponentialFamily x0
--     , Generative Natural x
--     ) =>
--     Int ->
--     -- | Harmonium
--     Natural # AffineHarmonium f x0 z0 x z ->
--     -- | Model Distribution
--     Natural # x ->
--     -- | Differential Estimate
--     Random (Mean # x)
-- harmoniumInformationProjectionDifferential n hrm px = do
--     xs <- sample n px
--     let (lkl, nw) = split hrm
--         mys0 = lkl >$>* xs
--         mws = sufficientStatistic <$> xs
--         mys = zipWith (\mw my0 -> mw <.> (px - nw) - potential my0) mws mys0
--         ln = fromIntegral $ length xs
--         mwht = average mws
--         myht = sum mys / ln
--         foldfun (mw, my) (k, z0) = (k + 1, z0 + ((my - myht) .> (mw - mwht)))
--     return . uncurry (/>) . foldr foldfun (-1, 0) $ zip mws mys
--

-- | Contrastive divergence on harmoniums (<https://www.mitpressjournals.org/doi/abs/10.1162/089976602760128018?casa_token=x_Twj1HaXcMAAAAA:7-Oq181aubCFwpG-f8Lo1wRKvGnmujzl8zjn9XbeO5nGhfvKCCQjsu4K4pJCkMNYUYWqc2qG7TRXBg Hinton, 2019>).
contrastiveDivergence ::
    ( KnownAffineHarmonium f x0 z0 x z
    , LegendreExponentialFamily z
    , Generative Natural z
    , Generative Natural x
    ) =>
    -- | The number of contrastive divergence steps
    Int ->
    -- | The initial states of the Gibbs chains
    Sample x ->
    -- | The harmonium
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | The gradient estimate
    Random (Mean # AffineHarmonium f x0 z0 x z)
contrastiveDivergence cdn xs hrm = do
    xzs0 <- initialPass hrm xs
    xzs1 <- iterateM' cdn (gibbsPass hrm) xzs0
    return $ stochasticRelativeEntropyDifferential xzs0 xzs1

--- Expectation Maximization ---

-- | A single iteration of EM for 'Harmonium' based models.
expectationMaximization ::
    ( KnownAffineHarmonium f x0 z0 x z
    , LegendreExponentialFamily z
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
    ( KnownAffineHarmonium f x0 z0 x z
    , LegendreExponentialFamily (AffineHarmonium f x0 z0 x z)
    , LegendreExponentialFamily z
    ) =>
    Double ->
    GradientPursuit ->
    Sample x ->
    Natural # AffineHarmonium f x0 z0 x z ->
    [Natural # AffineHarmonium f x0 z0 x z]
expectationMaximizationAscent eps gp xs nhrm =
    let mhrm' = expectationStep xs nhrm
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
    ( KnownAffineHarmonium f x0 z0 x z
    , LegendreExponentialFamily z
    , Generative Natural x
    , Generative Natural z
    ) =>
    Double ->
    Int ->
    Int ->
    GradientPursuit ->
    -- | Observations
    Sample x ->
    -- | Current Harmonium
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Harmonium Chain
    Chain Random (Natural # AffineHarmonium f x0 z0 x z)
gibbsExpectationMaximization eps cdn nbtch gp xs0 nhrm0 =
    let mhrm0 = expectationStep xs0 nhrm0
     in chainCircuit nhrm0 $ proc nhrm -> do
            xs <- minibatcher nbtch xs0 -< ()
            xzs0 <- arrM (uncurry initialPass) -< (nhrm, xs)
            xzs1 <- arrM (\(x, y) -> iterateM' cdn (gibbsPass x) y) -< (nhrm, xzs0)
            let dff = mhrm0 - averageSufficientStatistic xzs1
            gradientCircuit eps gp -< (nhrm, vanillaGradient dff)

latentProcessExpectationStep ::
    ( KnownLatentProcess f g x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , ConjugatedLikelihood f x0 z0 x z
    , LegendreExponentialFamily z
    , LegendreExponentialFamily (AffineHarmonium f x0 z0 x z)
    , LegendreExponentialFamily (AffineHarmonium g z0 z0 z z)
    ) =>
    Observations (LatentProcess f g x0 z0 x z) ->
    Natural # LatentProcess f g x0 z0 x z ->
    (Mean # z, Mean # AffineHarmonium f x0 z0 x z, Mean # AffineHarmonium g z0 z0 z z)
latentProcessExpectationStep zss ltnt =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        (smthss, hrmss) = unzip $ conjugatedSmoothing0 prr emsn trns <$> zss
        mprr = average $ toMean . head <$> smthss
        mtrns = average $ toMean <$> concat hrmss
        mws = toMean <$> concat smthss
        mzs = sufficientStatistic <$> concat zss
        mys = linearProjection <$> mzs
        mxs = linearProjection <$> mws
        memsn = joinHarmonium (average mzs) (mys >$< mxs) (average mws)
     in (mprr, memsn, mtrns)

latentProcessExpectationMaximization ::
    ( KnownLatentProcess f g x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , ConjugatedLikelihood f x0 z0 x z
    , DuallyFlatExponentialFamily z
    , DuallyFlatExponentialFamily (AffineHarmonium f x0 z0 x z)
    , DuallyFlatExponentialFamily (AffineHarmonium g z0 z0 z z)
    ) =>
    Observations (LatentProcess f g x0 z0 x z) ->
    Natural # LatentProcess f g x0 z0 x z ->
    Natural # LatentProcess f g x0 z0 x z
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
    ( KnownLatentProcess f g x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , ConjugatedLikelihood f x0 z0 x z
    , LegendreExponentialFamily z
    , LegendreExponentialFamily (AffineHarmonium f x0 z0 x z)
    , LegendreExponentialFamily (AffineHarmonium g z0 z0 z z)
    , DuallyFlatExponentialFamily z
    ) =>
    Double ->
    Int ->
    GradientPursuit ->
    Observations (LatentProcess f g x0 z0 x z) ->
    Natural # LatentProcess f g x0 z0 x z ->
    Natural # LatentProcess f g x0 z0 x z
latentProcessExpectationMaximizationAscent eps nstps gp zss ltnt =
    let (mprr, mehrm, mthrm) = latentProcessExpectationStep zss ltnt
        (nprr, nemsn, ntrns) = splitLatentProcess ltnt
        neql0 = toNatural . snd $ split mehrm
        neql1 = toNatural . snd $ split mthrm
        nehrm = joinConjugatedHarmonium nemsn neql0
        nthrm = joinConjugatedHarmonium ntrns neql1
        nprr' =
            (!! nstps) $
                vanillaGradientSequence (relativeEntropyDifferential mprr) (-eps) gp nprr
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

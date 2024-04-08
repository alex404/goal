{-# LANGUAGE Arrows #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | An Exponential Family 'Harmonium' is a product exponential family with a
particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
A 'Mixture' model is a special case of harmonium.
-}
module Goal.Graphical.Models.Harmonium.Approximate (
    -- * General
    conjugationParameterRegression,
    joinConjugatedHarmonium0,
    splitConjugatedHarmonium0,
    conjugatedBayesRule0,
    conjugatedRecursiveBayesianInference0,

    -- * Population Codes
    ProbabilisticPopulationCode,
    GaussianBoltzmannPopulationCode,
    KnownPopulationCode,
    joinPopulationCode,
    samplePPC,
    ppcExpectationBiases,
    ppcExpectationMaximizationAscent,
    ppcStochasticExpectationMaximization,
    ppcStochasticMaximumLikelihood,
    ppcLogLikelihood,
) where

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models.Harmonium

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.List (scanl')
import Data.Proxy (Proxy (..))

--- General ---

conjugationParameterRegression ::
    forall f x0 z x.
    (KnownAffineHarmonium f x0 z x z, Legendre x) =>
    Sample z ->
    Natural # Affine f x0 x z ->
    (Double, Natural # z)
{-# INLINE conjugationParameterRegression #-}
conjugationParameterRegression zs lkl =
    let sz0 :: SamplePoint z -> Mean # z
        sz0 = sufficientStatistic
        sz = S.cons 1 . coordinates . sz0
        ptns = map potential $ lkl >$>* zs
        bts = S.linearLeastSquares (sz <$> zs) ptns
        (rho0, rprms) = S.splitAt bts
     in (S.head rho0, Point rprms)

-- | The conjugation parameters of a conjugated `Harmonium`.
splitConjugatedHarmonium0 ::
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    Natural # AffineHarmonium f x0 z0 x z ->
    (Natural # Affine f x0 x z0, Natural # z)
{-# INLINE splitConjugatedHarmonium0 #-}
splitConjugatedHarmonium0 rprms hrm =
    let (lkl, nw) = split hrm
     in (lkl, nw >+> rprms)

-- | The conjugation parameters of a conjugated `Harmonium`.
joinConjugatedHarmonium0 ::
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    Natural # Affine f x0 x z0 ->
    Natural # z ->
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z
{-# INLINE joinConjugatedHarmonium0 #-}
joinConjugatedHarmonium0 rprms lkl prr = join lkl $ prr >+> (-rprms)

{- | The posterior distribution given a prior and likelihood, where the
likelihood is approximately conjugated.
-}
conjugatedBayesRule0 ::
    forall f x0 z0 x z.
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    Natural # Affine f x0 x z0 ->
    Natural # z ->
    SamplePoint x ->
    Natural # z
conjugatedBayesRule0 rprms lkl prr x =
    let hrm = joinConjugatedHarmonium0 rprms lkl prr
        pstr = fst . split $ transposeHarmonium hrm
        mx :: Mean # x
        mx = sufficientStatistic x
     in pstr >.+> mx

{- | The posterior distribution given a prior and likelihood, where the
likelihood is conjugated.
-}
conjugatedRecursiveBayesianInference0 ::
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    -- | Likelihood
    Natural # Affine f x0 x z0 ->
    -- | Prior
    Natural # z ->
    -- | Observations
    Sample x ->
    -- | Updated prior
    [Natural # z]
conjugatedRecursiveBayesianInference0 rho lkl = scanl' (conjugatedBayesRule0 rho lkl)

--- Population Codes ---

type ProbabilisticPopulationCode n y0 y = AffineHarmonium L.Full (Replicated n Poisson) y0 (Replicated n Poisson) y

type KnownPopulationCode n y0 y =
    ( KnownAffineHarmonium L.Full (Replicated n Poisson) y0 (Replicated n Poisson) y
    )

joinPopulationCode ::
    (KnownNat n, LegendreExponentialFamily x) =>
    -- | Gains
    Natural # Replicated n Poisson ->
    -- | Von Mises Curves
    S.Vector n (Natural # x) ->
    -- | Population Likelihood
    Natural # Replicated n Poisson <* x
{-# INLINE joinPopulationCode #-}
joinPopulationCode nz0 nps =
    let mtx = fromRows nps
        nz = nz0 - Point (S.map potential nps)
     in join nz mtx

--- Hierarchical Boltzman Machine

type GaussianBoltzmannPopulationCode f n m k =
    ProbabilisticPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k)

--- Functions ---

samplePPC ::
    forall n y0 y.
    (KnownPopulationCode n y0 y, Generative Natural y) =>
    Int ->
    Natural # y0 ->
    Natural # ProbabilisticPopulationCode n y0 y ->
    Random (Sample (ProbabilisticPopulationCode n y0 y))
{-# INLINE samplePPC #-}
samplePPC n rprms ppc = do
    let (lkl, gbhrm) = splitConjugatedHarmonium0 rprms ppc
    yzs <- sample n gbhrm
    let myzs :: [Mean # y]
        myzs = sufficientStatistic <$> yzs
    xs <- mapM samplePoint $ lkl >$> (linearProjection <$> myzs)
    return $ zip xs yzs

ppcExpectationBiases ::
    (KnownNat n, ExponentialFamily y) =>
    -- | Model Samples
    Sample (Replicated n Poisson) ->
    -- | Harmonium
    Natural # Replicated n Poisson <* y ->
    -- | Harmonium expected sufficient statistics
    [Natural # y]
{-# INLINE ppcExpectationBiases #-}
ppcExpectationBiases ns pc =
    let pstr0 = snd $ split pc
     in ns *<$< pstr0

ppcExpectationStep ::
    (Transition Natural Mean y, LinearSubspace y y0) =>
    [Natural # y0] ->
    Natural # y ->
    Mean # y
{-# INLINE ppcExpectationStep #-}
ppcExpectationStep ny0s ny =
    average $ toMean . (ny >+>) <$> ny0s

ppcExpectationMaximizationAscent ::
    (LinearSubspace y y0, LegendreExponentialFamily y) =>
    Double ->
    GradientPursuit ->
    [Natural # y0] ->
    Natural # y ->
    [Natural # y]
{-# INLINE ppcExpectationMaximizationAscent #-}
ppcExpectationMaximizationAscent eps gp ny0s ny =
    let my0 = ppcExpectationStep ny0s ny
     in vanillaGradientSequence (relativeEntropyDifferential my0) (-eps) gp ny

ppcStochasticExpectationMaximization ::
    (LinearSubspace y y0, Generative Natural y, LegendreExponentialFamily y) =>
    Double ->
    GradientPursuit ->
    Int ->
    [Natural # y0] ->
    (Natural # y) ->
    Chain Random (Natural # y)
{-# INLINE ppcStochasticExpectationMaximization #-}
ppcStochasticExpectationMaximization eps gp nbtch ny0s ny =
    let my = ppcExpectationStep ny0s ny
     in chainCircuit ny $ proc ny' -> do
            ys <- arrM (sample nbtch) -< ny'
            let dff = my - averageSufficientStatistic ys
            gradientCircuit eps gp -< (ny', vanillaGradient dff)

ppcStochasticMaximumLikelihood ::
    (LinearSubspace y y0, Generative Natural y, ExponentialFamily y) =>
    Double ->
    GradientPursuit ->
    Int ->
    [Natural # y0] ->
    (Natural # y) ->
    Chain Random (Natural # y)
{-# INLINE ppcStochasticMaximumLikelihood #-}
ppcStochasticMaximumLikelihood eps gp nbtch ny0s ny =
    chainCircuit ny $ proc ny1 -> do
        ny0 <- minibatcher 1 ny0s -< ()
        ys0 <- arrM (sample nbtch) -< ny1 >+> head ny0
        ys1 <- arrM (sample nbtch) -< ny1
        let dff = averageSufficientStatistic ys0 - averageSufficientStatistic ys1
        gradientCircuit eps gp -< (ny1, vanillaGradient dff)

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
ppcLogLikelihood ::
    forall n y0 y.
    (KnownPopulationCode n y0 y, Legendre y) =>
    (Double, Natural # y0) ->
    Sample (Replicated n Poisson) ->
    Natural # ProbabilisticPopulationCode n y0 y ->
    Double
{-# INLINE ppcLogLikelihood #-}
ppcLogLikelihood (rho0, rprms) xs ppc =
    let (pstr, nx) = split $ transposeHarmonium ppc
        mxs = sufficientStatistic <$> xs
        nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
        udns = zipWith (+) nrgs $ logBaseMeasure (Proxy @(Replicated n Poisson)) <$> xs
        nz = snd $ split ppc
     in average $ subtract (potential (nz >+> rprms) + rho0) <$> udns

--- Graveyard ---

-- A self contained, sampling based version of PPC EM.
-- ppcExpectationMaximization' ::
--     (KnownPopulationCode n y0 y, LegendreExponentialFamily y) =>
--     Sample (Replicated n Poisson) ->
--     Double ->
--     Int ->
--     GradientPursuit ->
--     Natural # y0 ->
--     Natural # ProbabilisticPopulationCode n y0 y ->
--     Chain Random (Natural # ProbabilisticPopulationCode n y0 y)
-- ppcExpectationMaximization' xs0 eps nsmps gp rprms nppc0 =
--     let mppc0 = expectationStep xs0 nppc0
--      in chainCircuit nppc0 $ proc nppc -> do
--             xyzs <- arrM (samplePPC nsmps rprms) -< nppc
--             let dff0 = mppc0 - averageSufficientStatistic xyzs
--                 dff = join 0 . snd $ split dff0
--             gradientCircuit eps gp -< (nppc, vanillaGradient dff)
--
-- sampleGBPPC ::
--     (KnownPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k)) =>
--     Int ->
--     Natural # MultivariateNormal f m ->
--     Natural # GaussianBoltzmannPopulationCode f n m k ->
--     Random [(S.Vector n Int, (S.Vector m Double, S.Vector k Bool))]
-- sampleGBPPC n rprms gbppc = do
--     let (ppc, gbhrm) = approximateSplitConjugatedHarmonium rprms gbppc
--     yzs <- sample n gbhrm
--     let ys = fst <$> yzs
--     xs <- mapM samplePoint $ ppc >$>* ys
--     return $ zip xs yzs
--
-- gbppcExpectationMaximization ::
--     ( KnownPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k)
--     , KnownCovariance f m
--     , KnownNat k
--     , 1 <= k
--     ) =>
--     Sample (Replicated n Poisson) ->
--     Double ->
--     Int ->
--     GradientPursuit ->
--     Natural # MultivariateNormal f m ->
--     Natural # ProbabilisticPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k) ->
--     Chain Random (Natural # ProbabilisticPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k))
-- gbppcExpectationMaximization xs0 eps nsmps gp rprms nppc0 =
--     let mppc0 = expectationStep xs0 nppc0
--      in chainCircuit nppc0 $ proc nppc -> do
--             xyzs <- arrM (sampleGBPPC nsmps rprms) -< nppc
--             let dff0 = mppc0 - averageSufficientStatistic xyzs
--                 dff = join 0 . snd $ split dff0
--             gradientCircuit eps gp -< (nppc, vanillaGradient dff)
--
--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
-- gbppcLogLikelihood ::
--     forall f n m k.
--     ( KnownPopulationCode n (MultivariateNormal f m) (GaussianBoltzmannHarmonium f m k)
--     , KnownCovariance f m
--     , KnownNat k
--     , 1 <= k
--     ) =>
--     (Double, Natural # MultivariateNormal f m) ->
--     Sample (Replicated n Poisson) ->
--     Natural # GaussianBoltzmannPopulationCode f n m k ->
--     Double
-- gbppcLogLikelihood (rho0, rprms) xs gbppc =
--     let (pstr, nx) = split $ transposeHarmonium gbppc
--         mxs = sufficientStatistic <$> xs
--         nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
--         udns = zipWith (+) nrgs $ logBaseMeasure (Proxy @(Replicated n Poisson)) <$> xs
--         nz = snd $ split gbppc
--      in average $ subtract (potential (nz >+> rprms) + rho0) <$> udns
--

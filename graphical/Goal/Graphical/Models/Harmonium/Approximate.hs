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
    approximateJoinConjugatedHarmonium,
    approximateSplitConjugatedHarmonium,

    -- * Population Codes
    ProbabilisticPopulationCode,
    GaussianBoltzmannPopulationCode,
    KnownPopulationCode,
    joinPopulationCode,
    samplePPC,
    ppcExpectationMaximization,
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

import Data.Proxy (Proxy (..))
import Debug.Trace

--- General ---

conjugationParameterRegression ::
    forall f x0 z x.
    (KnownAffineHarmonium f x0 z x z, Legendre x) =>
    Sample z ->
    Natural # Affine f x0 x z ->
    (Double, Natural # z)
conjugationParameterRegression zs lkl =
    let sz0 :: SamplePoint z -> Mean # z
        sz0 = sufficientStatistic
        sz = S.cons 1 . coordinates . sz0
        ptns = map potential $ lkl >$>* zs
        bts = S.linearLeastSquares (sz <$> zs) ptns
        (rho0, rprms) = S.splitAt bts
     in (S.head rho0, Point rprms)

-- | The conjugation parameters of a conjugated `Harmonium`.
approximateSplitConjugatedHarmonium ::
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    Natural # AffineHarmonium f x0 z0 x z ->
    (Natural # Affine f x0 x z0, Natural # z)
approximateSplitConjugatedHarmonium rprms hrm =
    let (lkl, nw) = split hrm
     in (lkl, nw >+> rprms)

-- | The conjugation parameters of a conjugated `Harmonium`.
approximateJoinConjugatedHarmonium ::
    (KnownAffineHarmonium f x0 z0 x z) =>
    Natural # z0 ->
    Natural # Affine f x0 x z0 ->
    Natural # z ->
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z
approximateJoinConjugatedHarmonium rprms lkl prr = join lkl $ prr >+> (-rprms)

--- Population Codes ---

type ProbabilisticPopulationCode n y0 y = AffineHarmonium L.Full (Replicated n Poisson) y0 (Replicated n Poisson) y

type KnownPopulationCode n y0 y =
    ( KnownAffineHarmonium L.Full (Replicated n Poisson) y0 (Replicated n Poisson) y
    , Generative Natural y
    )

joinPopulationCode ::
    (KnownNat n, LegendreExponentialFamily x) =>
    -- | Gains
    Natural # Replicated n Poisson ->
    -- | Von Mises Curves
    S.Vector n (Natural # x) ->
    -- | Population Likelihood
    Natural # Replicated n Poisson <* x
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
    (KnownPopulationCode n y0 y) =>
    Int ->
    Natural # y0 ->
    Natural # ProbabilisticPopulationCode n y0 y ->
    Random (Sample (ProbabilisticPopulationCode n y0 y))
samplePPC n rprms ppc = do
    let (lkl, gbhrm) = approximateSplitConjugatedHarmonium rprms ppc
    yzs <- sample n gbhrm
    let myzs :: [Mean # y]
        myzs = sufficientStatistic <$> yzs
    xs <- mapM samplePoint $ lkl >$> (linearProjection <$> myzs)
    return $ zip xs yzs

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
ppcExpectationMaximization ::
    (KnownPopulationCode n y0 y, LegendreExponentialFamily y) =>
    Sample (Replicated n Poisson) ->
    Double ->
    Int ->
    GradientPursuit ->
    Natural # y0 ->
    Natural # ProbabilisticPopulationCode n y0 y ->
    Chain Random (Natural # ProbabilisticPopulationCode n y0 y)
ppcExpectationMaximization xs0 eps nsmps gp rprms nppc0 =
    let mppc0 = expectationStep xs0 nppc0
     in chainCircuit nppc0 $ proc nppc -> do
            xyzs <- arrM (samplePPC nsmps rprms) -< nppc
            let dff0 = mppc0 - averageSufficientStatistic xyzs
                dff = join 0 . snd $ split dff0
            gradientCircuit eps gp -< (nppc, vanillaGradient dff)

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
ppcLogLikelihood ::
    forall n y0 y.
    (KnownPopulationCode n y0 y, Legendre y) =>
    (Double, Natural # y0) ->
    Sample (Replicated n Poisson) ->
    Natural # ProbabilisticPopulationCode n y0 y ->
    Double
ppcLogLikelihood (rho0, rprms) xs ppc =
    let (pstr, nx) = split $ transposeHarmonium ppc
        mxs = sufficientStatistic <$> xs
        nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
        udns = zipWith (+) nrgs $ logBaseMeasure (Proxy @(Replicated n Poisson)) <$> xs
        nz = snd $ split ppc
     in average $ subtract (potential (nz >+> rprms) + rho0) <$> udns

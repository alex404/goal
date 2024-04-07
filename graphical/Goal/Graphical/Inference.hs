{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=100 #-}

-- | Infering latent variables in graphical models.
module Goal.Graphical.Inference (
    -- * Inference
    conjugatedBayesRule,

    -- * Recursive
    conjugatedRecursiveBayesianInference,

    -- * Dynamic
    conjugatedPredictionStep,
    conjugatedForwardStep,

    -- * Conjugation
    regressConjugationParameters,
    conjugationCurve,
) where

--- Imports ---

-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models.Harmonium

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

import Data.Kind (Type)
import Data.List

--- Conjugate Prior Exponential Families ---

type family ConjugatePriorFamily x :: Type

-- type PoissonGammaLikelihood = Natural # Affine L.Identity Poisson Poisson Gamma

--- Inference ---

{- | The posterior distribution given a prior and likelihood, where the
likelihood is conjugated.
-}
conjugatedBayesRule ::
    forall f x0 z0 x z.
    (ConjugatedLikelihood f x0 z0 x z) =>
    Natural # Affine f x0 x z0 ->
    Natural # z ->
    SamplePoint x ->
    Natural # z
conjugatedBayesRule lkl prr x =
    let hrm = joinConjugatedHarmonium lkl prr
        pstr = fst . split $ transposeHarmonium hrm
        mx :: Mean # x
        mx = sufficientStatistic x
     in pstr >.+> mx

--- Recursive ---

{- | The posterior distribution given a prior and likelihood, where the
likelihood is conjugated.
-}
conjugatedRecursiveBayesianInference ::
    (ConjugatedLikelihood f x0 z0 x z) =>
    -- | Likelihood
    Natural # Affine f x0 x z0 ->
    -- | Prior
    Natural # z ->
    -- | Observations
    Sample x ->
    -- | Updated prior
    [Natural # z]
conjugatedRecursiveBayesianInference lkl = scanl' (conjugatedBayesRule lkl)

-- Dynamical ---

{- | The predicted distribution given a current distribution and transition
distribution, where the transition distribution is (doubly) conjugated.
-}
conjugatedPredictionStep ::
    (ConjugatedLikelihood f z0 z0 z z) =>
    Natural # Affine f z0 z z0 ->
    Natural # z ->
    Natural # z
conjugatedPredictionStep trns prr =
    snd
        . splitConjugatedHarmonium
        . transposeHarmonium
        $ joinConjugatedHarmonium trns prr

{- | Forward inference based on conjugated models: priors at a previous time are
first predicted into the current time, and then updated with Bayes rule.
-}
conjugatedForwardStep ::
    (ConjugatedLikelihood g z0 z0 z z, ConjugatedLikelihood f x0 z0 x z) =>
    -- | Transition Distribution
    Natural # Affine g z0 z z0 ->
    -- | Emission Distribution
    Natural # Affine f x0 x z0 ->
    -- | Beliefs at time $t-1$
    Natural # z ->
    -- | Observation at time $t$
    SamplePoint x ->
    -- | Beliefs at time $t$
    Natural # z
conjugatedForwardStep trns emsn prr z =
    flip (conjugatedBayesRule emsn) z $ conjugatedPredictionStep trns prr

--- Approximate Conjugation ---

{- | Computes the conjugation curve given a set of conjugation parameters,
at the given set of points.
-}
conjugationCurve ::
    (ExponentialFamily z) =>
    -- | Conjugation shift
    Double ->
    -- | Conjugation parameters
    Natural # z ->
    -- | Samples points
    Sample z ->
    -- | Conjugation curve at sample points
    [Double]
conjugationCurve rho0 rprms mus = (\z -> rprms <.> sufficientStatistic z + rho0) <$> mus

-- Linear Least Squares

{- | Returns the conjugation parameters which best satisfy the conjugation
equation for the given population code according to linear regression.
-}
regressConjugationParameters ::
    (Map Natural f z x, LegendreExponentialFamily z, ExponentialFamily x) =>
    -- | PPC
    Natural # f z x ->
    -- | Sample points
    Sample x ->
    -- | Approximate conjugation parameters
    (Double, Natural # x)
regressConjugationParameters lkl mus =
    let dpnds = potential <$> lkl >$>* mus
        indpnds = independentVariables0 lkl mus
        (rho0, rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

--- Internal ---

independentVariables0 ::
    forall f x z.
    (ExponentialFamily x) =>
    Natural # f z x ->
    Sample x ->
    [S.Vector (Dimension x + 1) Double]
independentVariables0 _ mus =
    let sss :: [Mean # x]
        sss = sufficientStatistic <$> mus
     in (S.singleton 1 S.++) . coordinates <$> sss

---- | The posterior distribution given a prior and likelihood, where the
---- posterior is normalized via numerical integration.
-- numericalRecursiveBayesianInference
--    :: forall f z x .
--        ( Map Natural f x z, Map Natural f z x, Bilinear f z x
--        , LegendreExponentialFamily z, ExponentialFamily x, SamplePoint x ~ Double)
--    => Double -- ^ Integral error bound
--    -> Double -- ^ Sample space lower bound
--    -> Double -- ^ Sample space upper bound
--    -> Sample x -- ^ Centralization samples
--    -> [Natural # Affine f z x] -- ^ Likelihoods
--    -> Sample z -- ^ Observations
--    -> (Double -> Double) -- ^ Prior
--    -> (Double -> Double, Double) -- ^ Posterior Density and Log-Partition Function
-- numericalRecursiveBayesianInference errbnd mnx mxx xsmps lkls zs prr =
--    let logbm = logBaseMeasure (Proxy @ x)
--        logupst0 x lkl z =
--            (z *<.< snd (splitAffine lkl)) <.> sufficientStatistic x - potential (lkl >.>* x)
--        logupst x = sum $ logbm x : log (prr x) : zipWith (logupst0 x) lkls zs
--        logprt = logIntegralExp errbnd logupst mnx mxx xsmps
--        dns x = exp $ logupst x - logprt
--     in (dns,logprt)

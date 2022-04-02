{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    RankNTypes,
    PolyKinds,
    DataKinds,
    TypeOperators,
    FlexibleContexts,
    FlexibleInstances,
    TypeApplications,
    ScopedTypeVariables,
    TypeFamilies
#-}
-- | Infering latent variables in graphical models.
module Goal.Graphical.Inference
    ( -- * Inference
      conjugatedBayesRule
    -- * Recursive
    , conjugatedRecursiveBayesianInference
    -- * Dynamic
    , conjugatedPredictionStep
    , conjugatedForwardStep
    -- * Conjugation
    , regressConjugationParameters
    , conjugationCurve
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models.Harmonium

import qualified Goal.Core.Vector.Storable as S

import Data.List


--- Inference ---


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedBayesRule
    :: forall f y x z w
    . ( Bilinear Natural f y x, ConjugatedLikelihood f y x z w )
    => Natural # Affine f y z x
    -> Natural # w
    -> SamplePoint z
    -> Natural # w
conjugatedBayesRule lkl prr z =
    let pstr = fst . split . transposeHarmonium $ joinConjugatedHarmonium lkl prr
        mz :: Mean # z
        mz = sufficientStatistic z
     in pstr >.+> mz


--- Recursive ---


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedRecursiveBayesianInference
    :: ( Bilinear Natural f y x, ConjugatedLikelihood f y x z w )
    => Natural # Affine f y z x -- ^ Likelihood
    -> Natural # w -- ^ Prior
    -> Sample z -- ^ Observations
    -> [Natural # w] -- ^ Updated prior
conjugatedRecursiveBayesianInference lkl = scanl' (conjugatedBayesRule lkl)


-- Dynamical ---


-- | The predicted distribution given a current distribution and transition
-- distribution, where the transition distribution is (doubly) conjugated.
conjugatedPredictionStep
    :: (ConjugatedLikelihood f x x w w, Bilinear Natural f x x)
    => Natural # Affine f x w x
    -> Natural # w
    -> Natural # w
conjugatedPredictionStep trns prr =
    snd . splitConjugatedHarmonium . transposeHarmonium
        $ joinConjugatedHarmonium trns prr

-- | Forward inference based on conjugated models: priors at a previous time are
-- first predicted into the current time, and then updated with Bayes rule.
conjugatedForwardStep
    :: ( ConjugatedLikelihood g x x w w, Bilinear Natural g x x
       , ConjugatedLikelihood f y x z w, Bilinear Natural f y x
       , Map Natural g x x, Map Natural f x y )
    => Natural # Affine g x w x -- ^ Transition Distribution
    -> Natural # Affine f y z x -- ^ Emission Distribution
    -> Natural # w -- ^ Beliefs at time $t-1$
    -> SamplePoint z -- ^ Observation at time $t$
    -> Natural # w -- ^ Beliefs at time $t$
conjugatedForwardStep trns emsn prr z =
    flip (conjugatedBayesRule emsn) z $ conjugatedPredictionStep trns prr


--- Approximate Conjugation ---


-- | Computes the conjugation curve given a set of conjugation parameters,
-- at the given set of points.
conjugationCurve
    :: ExponentialFamily x
    => Double -- ^ Conjugation shift
    -> Natural # x -- ^ Conjugation parameters
    -> Sample x -- ^ Samples points
    -> [Double] -- ^ Conjugation curve at sample points
conjugationCurve rho0 rprms mus = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- Linear Least Squares

-- | Returns the conjugation parameters which best satisfy the conjugation
-- equation for the given population code according to linear regression.
regressConjugationParameters
    :: (Map Natural f z x, LegendreExponentialFamily z, ExponentialFamily x)
    => Natural # f z x -- ^ PPC
    -> Sample x -- ^ Sample points
    -> (Double, Natural # x) -- ^ Approximate conjugation parameters
regressConjugationParameters lkl mus =
    let dpnds = potential <$> lkl >$>* mus
        indpnds = independentVariables0 lkl mus
        (rho0,rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

--- Internal ---

independentVariables0
    :: forall f x z . ExponentialFamily x
    => Natural # f z x
    -> Sample x
    -> [S.Vector (Dimension x + 1) Double]
independentVariables0 _ mus =
    let sss :: [Mean # x]
        sss = sufficientStatistic <$> mus
     in (S.singleton 1 S.++) . coordinates <$> sss


---- | The posterior distribution given a prior and likelihood, where the
---- posterior is normalized via numerical integration.
--numericalRecursiveBayesianInference
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
--numericalRecursiveBayesianInference errbnd mnx mxx xsmps lkls zs prr =
--    let logbm = logBaseMeasure (Proxy @ x)
--        logupst0 x lkl z =
--            (z *<.< snd (splitAffine lkl)) <.> sufficientStatistic x - potential (lkl >.>* x)
--        logupst x = sum $ logbm x : log (prr x) : zipWith (logupst0 x) lkls zs
--        logprt = logIntegralExp errbnd logupst mnx mxx xsmps
--        dns x = exp $ logupst x - logprt
--     in (dns,logprt)



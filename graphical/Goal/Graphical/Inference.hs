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
-- | Exponential Family Harmoniums and Conjugation.
module Goal.Graphical.Inference
    ( -- * Inference
      conjugatedBayesRule
    -- * Recursive
    , conjugatedRecursiveBayesianInference
    -- * Dynamic
    , conjugatedPredictionStep
    , conjugatedForwardStep
    , conjugatedFiltering
    , conjugatedSmoothing
    -- * Conjugation
    , regressConjugationParameters
    , conjugationCurve
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Conditional
import Goal.Graphical.Generative.Harmonium

import qualified Goal.Core.Vector.Storable as S

import Data.List


--- Inference ---


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedBayesRule
    :: ( Map Natural f x z, ExponentialFamily z, ExponentialFamily x
       , Bilinear f z x, ConjugatedLikelihood f z x )
    => Natural # Affine f z x
    -> Natural # x
    -> SamplePoint z
    -> Natural # x
conjugatedBayesRule lkl prr z =
    transposeHarmonium (joinConjugatedHarmonium lkl prr) >.>* z


--- Recursive ---


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedRecursiveBayesianInference
    :: ( Map Natural f x z, ExponentialFamily z, ExponentialFamily x
       , Bilinear f z x, ConjugatedLikelihood f z x )
    => Natural # Affine f z x -- ^ Likelihood
    -> Natural # x -- ^ Prior
    -> Sample z -- ^ Observations
    -> [Natural # x] -- ^ Updated prior
conjugatedRecursiveBayesianInference lkl = scanl' (conjugatedBayesRule lkl)


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


-- Dynamical ---


conjugatedPredictionStep
    :: (ConjugatedLikelihood f x x, Bilinear f x x)
    => Natural # Affine f x x
    -> Natural # x
    -> Natural # x
conjugatedPredictionStep trns prr =
    snd . splitConjugatedHarmonium . transposeHarmonium
        $ joinConjugatedHarmonium trns prr

conjugatedForwardStep
    :: ( ExponentialFamily z, ExponentialFamily x, ConjugatedLikelihood g z x
       , ConjugatedLikelihood f x x, Bilinear f x x, Bilinear g z x
       , Map Natural g x z)
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> Natural # x
    -> SamplePoint z
    -> Natural # x
conjugatedForwardStep trns emsn prr z =
    flip (conjugatedBayesRule emsn) z $ conjugatedPredictionStep trns prr

conjugatedFiltering
    :: ( ExponentialFamily z, ExponentialFamily x, ConjugatedLikelihood g z x
       , ConjugatedLikelihood f x x, Bilinear f x x, Bilinear g z x
       , Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> Natural # x
    -> Sample z
    -> [Natural # x]
conjugatedFiltering trns emsn prr zs =
    tail $ scanl' (conjugatedForwardStep trns emsn) prr zs

conjugatedBackwardStep
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, Map Natural g x z, Bilinear g z x)
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> Natural # x
    -> Natural # x
    -> SamplePoint z
    -> (Natural # x, Natural # x)
conjugatedBackwardStep trns emsn dff flt z =
    let (tx,txx) = splitAffine trns
        ezx = snd $ splitAffine emsn
        ecnj = snd $ conjugationParameters emsn
        trns' = joinAffine (dff + tx - ecnj + z *<.< ezx) txx
        dff' = snd (conjugationParameters trns') - snd (conjugationParameters trns)
     in (flt + dff', dff')

conjugatedSmoothing
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Bilinear g z x, Map Natural g x z )
    => Natural # x
    -> Natural # Affine f x x
    -> Natural # Affine g z x
    -> Sample z
    -> [Natural # x]
conjugatedSmoothing prr trns emsn zs =
    let flts = conjugatedFiltering trns emsn prr zs
        (flt:flts') = reverse flts
     in fst <$> scanr scanner (flt,0) (zip zs $ reverse flts')
        where scanner (z,flt) (_,dff) =
                conjugatedBackwardStep trns emsn dff flt z


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
-- equation for the given population code.
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

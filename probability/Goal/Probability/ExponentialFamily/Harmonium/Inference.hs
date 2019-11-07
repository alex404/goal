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
module Goal.Probability.ExponentialFamily.Harmonium.Inference
    ( -- * Inference
      (<|<)
    , (<|<*)
    , numericalRecursiveBayesianInference
    -- ** Conjugated
    , conjugatedBayesRule
    , conjugatedRecursiveBayesianInference
    , conjugatedRecursiveBayesianInference'
    -- * Conjugation
    , regressConjugationParameters
    , conjugationCurve
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Inference ---


-- | The given deep harmonium conditioned on a mean distribution over the bottom layer.
(<|<) :: ( Bilinear f z y, Map Mean Natural f y z, Manifold (DeepHarmonium z fxs) )
      => Natural # DeepHarmonium z ('(f,y) : fxs) -- ^ Deep harmonium
      -> Mean # z -- ^ Input means
      -> Natural # DeepHarmonium y fxs -- ^ Conditioned deep harmonium
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (p <.< snd (splitAffine f)) dhrm'

-- | The given deep harmonium conditioned on a sample from its bottom layer.
-- In other words, the posterior of the model given an observation of
-- the observable variable.
(<|<*) :: ( Bilinear f z y, Map Mean Natural f y z
          , Manifold (DeepHarmonium z fxs), ExponentialFamily z )
      => Natural # DeepHarmonium z ('(f,y) : fxs) -- ^ Deep harmonium
      -> SamplePoint z -- ^ Input means
      -> Natural # DeepHarmonium y fxs -- ^ Conditioned deep harmonium
{-# INLINE (<|<*) #-}
(<|<*) dhrm x = dhrm <|< sufficientStatistic x

-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedBayesRule
    :: (Map Mean Natural f y z, Bilinear f z y, ExponentialFamily z)
    => Natural # y -- ^ Conjugation Parameters
    -> Natural #> Affine f z y -- ^ Likelihood
    -> SamplePoint z -- ^ Observation
    -> Natural # DeepHarmonium y fxs -- ^ Prior
    -> Natural # DeepHarmonium y fxs -- ^ Updated prior
{-# INLINE conjugatedBayesRule #-}
conjugatedBayesRule rprms lkl z =
    biasBottom (z *<.< snd (splitAffine lkl) - rprms)

-- | The posterior distribution given a prior and likelihood, where the
-- posterior is normalized via numerical integration.
numericalRecursiveBayesianInference
    :: forall f z x .
        ( Map Mean Natural f x z, Map Mean Natural f z x, Bilinear f z x
        , LegendreExponentialFamily z, ExponentialFamily x, SamplePoint x ~ Double)
    => Double -- ^ Integral error bound
    -> Double -- ^ Sample space lower bound
    -> Double -- ^ Sample space upper bound
    -> Sample x -- ^ Centralization samples
    -> [Natural #> Affine f z x] -- ^ Likelihoods
    -> Sample z -- ^ Observations
    -> (Double -> Double) -- ^ Prior
    -> (Double -> Double, Double) -- ^ Posterior Density and Log-Partition Function
{-# INLINE numericalRecursiveBayesianInference #-}
numericalRecursiveBayesianInference errbnd mnx mxx xsmps lkls zs prr =
    let logbm = log . baseMeasure (Proxy @ x)
        logupst0 x lkl z =
            (z *<.< snd (splitAffine lkl)) <.> sufficientStatistic x - potential (lkl >.>* x)
        logupst x = sum $ logbm x : log (prr x) : zipWith (logupst0 x) lkls zs
        logprt = logIntegralExp errbnd logupst mnx mxx xsmps
        dns x = exp $ logupst x - logprt
     in (dns,logprt)

-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedRecursiveBayesianInference'
    :: (Map Mean Natural f x z, Bilinear f z x, ExponentialFamily z)
    => Natural # x -- ^ Conjugation Parameters
    -> Natural #> Affine f z x -- ^ Likelihood
    -> Sample z -- ^ Observations
    -> Natural # x -- ^ Prior
    -> Natural # x -- ^ Posterior
{-# INLINE conjugatedRecursiveBayesianInference' #-}
conjugatedRecursiveBayesianInference' rprms lkl zs prr =
    let pstr0 = foldr (+) 0 $ (subtract rprms) <$> zs *<$< snd (splitAffine lkl)
     in pstr0 + prr


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is conjugated.
conjugatedRecursiveBayesianInference
    :: (Map Mean Natural f y z, Bilinear f z y, ExponentialFamily z)
    => [Natural # y] -- ^ Conjugation Parameters
    -> [Natural #> Affine f z y] -- ^ Likelihood
    -> Sample z -- ^ Observations
    -> Natural # DeepHarmonium y fxs -- ^ Prior
    -> Natural # DeepHarmonium y fxs -- ^ Updated prior
{-# INLINE conjugatedRecursiveBayesianInference #-}
conjugatedRecursiveBayesianInference rprmss lkls zs prr =
    foldl' (\pstr' (rprms,lkl,z) -> conjugatedBayesRule rprms lkl z pstr') prr (zip3 rprmss lkls zs)



-- | Computes the conjugation curve given a set of conjugation parameters,
-- at the given set of points.
conjugationCurve
    :: ExponentialFamily x
    => Double -- ^ Conjugation shift
    -> Natural # x -- ^ Conjugation parameters
    -> Sample x -- ^ Samples points
    -> [Double] -- ^ Conjugation curve at sample points
{-# INLINE conjugationCurve #-}
conjugationCurve rho0 rprms mus = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- Linear Least Squares

-- | Returns the conjugation parameters which best satisfy the conjugation
-- equation for the given population code.
regressConjugationParameters
    :: (Map Mean Natural f z x, LegendreExponentialFamily z, ExponentialFamily x)
    => Natural #> f z x -- ^ PPC
    -> Sample x -- ^ Sample points
    -> (Double, Natural # x) -- ^ Approximate conjugation parameters
{-# INLINE regressConjugationParameters #-}
regressConjugationParameters lkl mus =
    let dpnds = potential <$> lkl >$>* mus
        indpnds = independentVariables0 lkl mus
        (rho0,rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

--- Internal ---

independentVariables0
    :: forall f x z . ExponentialFamily x
    => Natural #> f z x
    -> Sample x
    -> [S.Vector (Dimension x + 1) Double]
independentVariables0 _ mus =
    let sss :: [Mean # x]
        sss = sufficientStatistic <$> mus
     in (S.singleton 1 S.++) . coordinates <$> sss




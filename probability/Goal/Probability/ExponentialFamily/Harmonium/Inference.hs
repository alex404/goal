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
-- | Exponential Family Harmoniums and Rectification.
module Goal.Probability.ExponentialFamily.Harmonium.Inference
    ( -- * Inference
      (<|<)
    , (<|<*)
    , numericalRecursiveBayesianInference'
    -- ** Rectified
    , rectifiedBayesRule
    , rectifiedRecursiveBayesianInference
    , rectifiedRecursiveBayesianInference'
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium

import Data.List (foldl')


--- Inference ---


-- | The given deep harmonium conditioned on a mean distribution over the bottom layer.
(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> Mean # m -- ^ Input means
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Conditioned deep harmonium
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (p <.< snd (splitAffine f)) dhrm'

-- | The given deep harmonium conditioned on a sample from its bottom layer.
-- This can be interpreted as the posterior of the model given an observation of
-- the bottom layer.
(<|<*) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)), ExponentialFamily m )
      => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> SamplePoint m -- ^ Observations
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Posterior
{-# INLINE (<|<*) #-}
(<|<*) dhrm x = dhrm <|< sufficientStatistic x

-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is rectified.
rectifiedBayesRule
    :: (Map Mean Natural f z y, Bilinear f z y, ExponentialFamily z)
    => Natural # y -- ^ Rectification Parameters
    -> Mean #> Natural # Affine f z y -- ^ Likelihood
    -> SamplePoint z -- ^ Observation
    -> Natural # DeepHarmonium gs (y : xs) -- ^ Prior
    -> Natural # DeepHarmonium gs (y : xs) -- ^ Updated prior
{-# INLINE rectifiedBayesRule #-}
rectifiedBayesRule rprms lkl z prr =
    biasBottom (z *<.< snd (splitAffine lkl) <-> rprms) prr

-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is rectified.
numericalRecursiveBayesianInference'
    :: forall f z x . ( Map Mean Natural f z x, Bilinear f z x, Legendre Natural z
                      , ExponentialFamily z, ExponentialFamily x, SamplePoint x ~ Double)
    => Double -- ^ Integral error bound
    -> Double -- ^ Sample space lower bound
    -> Double -- ^ Sample space upper bound
    -> Int -- ^ Number of centralization samples
    -> Mean #> Natural # Affine f z x -- ^ Likelihoods
    -> Sample z -- ^ Observations
    -> (Double -> Double) -- ^ Prior
    -> (Double -> Double, Double) -- ^ Posterior Density and Log-Partition Function
{-# INLINE numericalRecursiveBayesianInference' #-}
numericalRecursiveBayesianInference' errbnd mnx mxx ncntrs lkl zs prr =
    let lnr = snd (splitAffine lkl)
        logbm = log . baseMeasure (Proxy @ x)
        logupst0 x z =
            (z *<.< lnr) <.> sufficientStatistic x - potential (lkl >.>* x) - logbm x
        logupst x = sum $ log (prr x) : (logupst0 x <$> zs)
        logprt = logIntegralExp errbnd logupst mnx mxx (range mnx mxx ncntrs)
        dns x = exp $ logupst x - logprt
     in (dns,logprt)


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is rectified.
rectifiedRecursiveBayesianInference'
    :: (Map Mean Natural f z x, Bilinear f z x, ExponentialFamily z)
    => Natural # x -- ^ Rectification Parameters
    -> Mean #> Natural # Affine f z x -- ^ Likelihood
    -> Sample z -- ^ Observations
    -> Natural # x -- ^ Prior
    -> Natural # x -- ^ Posterior
{-# INLINE rectifiedRecursiveBayesianInference' #-}
rectifiedRecursiveBayesianInference' rprms lkl zs prr =
    let pstr0 = foldr (<+>) zero $ (<-> rprms) <$> zs *<$< snd (splitAffine lkl)
     in pstr0 <+> prr


-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is rectified.
rectifiedRecursiveBayesianInference
    :: (Map Mean Natural f z y, Bilinear f z y, ExponentialFamily z)
    => [Natural # y] -- ^ Rectification Parameters
    -> [Mean #> Natural # Affine f z y] -- ^ Likelihood
    -> Sample z -- ^ Observations
    -> Natural # DeepHarmonium gs (y : xs) -- ^ Prior
    -> Natural # DeepHarmonium gs (y : xs) -- ^ Updated prior
{-# INLINE rectifiedRecursiveBayesianInference #-}
rectifiedRecursiveBayesianInference rprmss lkls zs prr =
    foldl' (\pstr' (rprms,lkl,z) -> rectifiedBayesRule rprms lkl z pstr') prr (zip3 rprmss lkls zs)




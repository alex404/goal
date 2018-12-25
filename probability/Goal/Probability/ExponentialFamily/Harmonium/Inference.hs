{-# LANGUAGE
    RankNTypes,
    PolyKinds,
    DataKinds,
    TypeOperators,
    MultiParamTypeClasses,
    FlexibleContexts,
    FlexibleInstances,
    TypeFamilies,
    TypeApplications,
    ScopedTypeVariables,
    UndecidableInstances
#-}
-- | Exponential Family Harmoniums and Rectification.
module Goal.Probability.ExponentialFamily.Harmonium.Inference
    ( -- * Inference
      (<|<)
    , (<|<*)
    ) where

--- Imports ---


-- Goal --

import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium


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

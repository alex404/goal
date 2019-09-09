{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    ( -- * Types
      ConditionalModel
    , ConditionalDeepHarmonium
    , ConditionalHarmonium
    , ConditionalMixture
    -- * Construction
    , joinConditionalDeepHarmonium
    , splitConditionalDeepHarmonium
    -- * Analysis
    , conditionalHarmoniumEmpiricalExpectations
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Storable as S
import Goal.Probability.ExponentialFamily.Harmonium


--- Types ---


-- | A generic conditional model.
data ConditionalModel (f :: Type -> Type -> Type) z x

-- | A conditional 'DeepHarmonium', where the observbale biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalDeepHarmonium y (gxs :: [(Type -> Type -> Type,Type)]) f z
  = ConditionalModel f (DeepHarmonium y gxs) z

-- | A conditional 'Harmonium', where the observbale biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium y g x f z = ConditionalDeepHarmonium y '[ '(g,x)] f z

-- | A conditional 'Mixture', where the observbale biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture y k f z = ConditionalHarmonium y Tensor (Categorical k) f z -- ^ Function

-- | Splits a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalDeepHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c #> ConditionalDeepHarmonium y gxs f z -- ^ Conditional Harmonium
    -> (c # DeepHarmonium y gxs, c #> f y z) -- ^ Matrix function and upper part
{-# INLINE splitConditionalDeepHarmonium #-}
splitConditionalDeepHarmonium dhrm =
    let (dhrmcs,fcs) = S.splitAt $ coordinates dhrm
     in (Point dhrmcs,Point fcs)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalDeepHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c # DeepHarmonium y gxs
    -> c #> f y z
    -> c #> ConditionalDeepHarmonium y gxs f z -- ^ Conditional Harmonium
{-# INLINE joinConditionalDeepHarmonium #-}
joinConditionalDeepHarmonium (Point dcs) (Point fcs) = Point $ dcs S.++ fcs


-- | Empirical expectations of a conditional harmonium.
conditionalHarmoniumEmpiricalExpectations
    :: ( ExponentialFamily y, Bilinear g y x, Map Mean Natural g x y
       , Manifold (f y z), LegendreExponentialFamily x )
    => Sample y -- ^ Model Samples
    -> Natural #> ConditionalHarmonium y g x f z -- ^ Harmonium
    -> [Mean # Harmonium y g x] -- ^ Harmonium expected sufficient statistics
{-# INLINE conditionalHarmoniumEmpiricalExpectations #-}
conditionalHarmoniumEmpiricalExpectations ys hrm =
    let aff = fst . splitBottomHarmonium . transposeHarmonium . fst $ splitConditionalDeepHarmonium hrm
     in harmoniumEmpiricalExpectations0 (sufficientStatistic <$> ys) aff


--- Instances ---


instance (Map Mean Natural f y z, Manifold (DeepHarmonium y gxs))
    => Manifold (ConditionalDeepHarmonium y gxs f z) where
        type Dimension (ConditionalDeepHarmonium y gxs f z)
          = Dimension (DeepHarmonium y gxs) + Dimension (f y z)

instance ( Map Mean Natural f y z, Manifold (DeepHarmonium y gxs) )
     => Map Mean Natural (ConditionalModel f) (DeepHarmonium y gxs) z where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

instance (Propagate Mean Natural f y z, Manifold (DeepHarmonium y gxs))
  => Propagate Mean Natural (ConditionalModel f) (DeepHarmonium y gxs) z where
        {-# INLINE propagate #-}
        propagate dhrms dzs chrm =
            let dys = getBottomBias <$> dhrms
                (hrm,f) = splitConditionalDeepHarmonium chrm
                (df,hrmhts) = propagate dys dzs f
             in (joinConditionalDeepHarmonium (average dhrms) df, flip biasBottom hrm <$> hrmhts)


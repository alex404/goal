{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    ( -- * Types
      ConditionalDeepHarmonium
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


data SubLinear (f :: Type -> Type -> Type) z x

type ConditionalDeepHarmonium y (gxs :: [(Type -> Type -> Type,Type)]) f z
  = SubLinear f (DeepHarmonium y gxs) z

type ConditionalHarmonium y g x f z = ConditionalDeepHarmonium y '[ '(g,x)] f z

type ConditionalMixture y k z = ConditionalHarmonium y Tensor (Categorical k) Tensor z -- ^ Function


-- | Splits the top layer off of a harmonium.
splitConditionalDeepHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c #> ConditionalDeepHarmonium y gxs f z -- ^ Conditional Harmonium
    -> (c # DeepHarmonium y gxs, c #> f y z) -- ^ Matrix function and upper part
{-# INLINE splitConditionalDeepHarmonium #-}
splitConditionalDeepHarmonium dhrm =
    let (dhrmcs,fcs) = S.splitAt $ coordinates dhrm
     in (Point dhrmcs,Point fcs)

-- | Splits the top layer off of a harmonium.
joinConditionalDeepHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c # DeepHarmonium y gxs
    -> c #> f y z
    -> c #> ConditionalDeepHarmonium y gxs f z -- ^ Conditional Harmonium
{-# INLINE joinConditionalDeepHarmonium #-}
joinConditionalDeepHarmonium (Point dcs) (Point fcs) = Point $ dcs S.++ fcs


-- | This might be inefficient due to the use of average point instead of a
-- slightly more complicated foldr.
conditionalHarmoniumEmpiricalExpectations
    :: ( ExponentialFamily y, Bilinear f y x, Map Mean Natural f x y
       , Manifold (f y z), LegendreExponentialFamily x )
    => Sample y -- ^ Model Samples
    -> Natural #> ConditionalHarmonium y f x f z -- ^ Harmonium
    -> [Mean # Harmonium y f x] -- ^ Harmonium expected sufficient statistics
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
     => Map Mean Natural (SubLinear f) (DeepHarmonium y gxs) z where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

instance (Propagate Mean Natural f y z, Manifold (DeepHarmonium y gxs))
  => Propagate Mean Natural (SubLinear f) (DeepHarmonium y gxs) z where
        {-# INLINE propagate #-}
        propagate dhrms dzs chrm =
            let dys = getBottomBias <$> dhrms
                (hrm,f) = splitConditionalDeepHarmonium chrm
                (df,hrmhts) = propagate dys dzs f
             in (joinConditionalDeepHarmonium (averagePoint dhrms) df, flip biasBottom hrm <$> hrmhts)


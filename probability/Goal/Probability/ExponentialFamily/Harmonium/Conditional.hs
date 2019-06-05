{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    ( -- * Types
      ConditionalHarmonium
    , ConditionalMixture
    -- * Construction
    , joinConditionalHarmonium
    , splitConditionalHarmonium
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Storable as S
import Goal.Probability.ExponentialFamily.Harmonium


--- Types ---


data SubLinear (f :: Type -> Type -> Type) z x

type ConditionalHarmonium y (gxs :: [(Type -> Type -> Type,Type)]) f z
  = SubLinear f (DeepHarmonium y gxs) z

type ConditionalMixture y k z = ConditionalHarmonium y '[ '(Tensor,Categorical k)] Tensor z -- ^ Function


-- | Splits the top layer off of a harmonium.
splitConditionalHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => Dual c #> c # ConditionalHarmonium y gxs f z -- ^ Conditional Harmonium
    -> (c # DeepHarmonium y gxs, Dual c #> c # f y z) -- ^ Matrix function and upper part
{-# INLINE splitConditionalHarmonium #-}
splitConditionalHarmonium dhrm =
    let (dhrmcs,fcs) = S.splitAt $ coordinates dhrm
     in (Point dhrmcs,Point fcs)

-- | Splits the top layer off of a harmonium.
joinConditionalHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c # DeepHarmonium y gxs
    -> Dual c #> c # f y z
    -> Dual c #> c # ConditionalHarmonium y gxs f z -- ^ Conditional Harmonium
{-# INLINE joinConditionalHarmonium #-}
joinConditionalHarmonium (Point dcs) (Point fcs) = Point $ dcs S.++ fcs


--- Instances ---


instance (Map Mean Natural f y z, Manifold (DeepHarmonium y gxs))
    => Manifold (ConditionalHarmonium y gxs f z) where
        type Dimension (ConditionalHarmonium y gxs f z)
          = Dimension (DeepHarmonium y gxs) + Dimension (f y z)

instance ( Map Mean Natural f y z, Manifold (DeepHarmonium y gxs) )
     => Map Mean Natural (SubLinear f) (DeepHarmonium y gxs) z where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitConditionalHarmonium pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitConditionalHarmonium pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

instance (Propagate Mean Natural f y z, Manifold (DeepHarmonium y gxs))
  => Propagate Mean Natural (SubLinear f) (DeepHarmonium y gxs) z where
        {-# INLINE propagate #-}
        propagate dhrms dzs chrm =
            let dys = getBottomBias <$> dhrms
                (hrm,f) = splitConditionalHarmonium chrm
                (df,hrmhts) = propagate dys dzs f
             in (joinConditionalHarmonium (averagePoint dhrms) df, flip biasBottom hrm <$> hrmhts)


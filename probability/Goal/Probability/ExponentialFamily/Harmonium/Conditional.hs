{-# LANGUAGE
    KindSignatures,
    DataKinds,
    TypeOperators,
    FlexibleInstances,
    TypeFamilies,
    MultiParamTypeClasses,
    FlexibleContexts,
    UndecidableInstances
#-}

module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    ( -- * Types
    SubLinear
    , ConditionalHarmonium
    , MixtureGLM
    -- * Construction
    , joinBottomSubLinear
    , splitBottomSubLinear
    -- * Computation
    , mixtureStochasticConditionalCrossEntropy
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

type ConditionalHarmonium f gs ns x = SubLinear f (DeepHarmonium gs ns) x

type MixtureGLM z e k x =
    ConditionalHarmonium Tensor '[Tensor] [z, Categorical e k] x -- ^ Function

-- | Splits the top layer off of a harmonium.
splitBottomSubLinear
    :: (Manifold x, Manifold (f z x), Manifold (DeepHarmonium gs (z : zs)))
    => Dual c #> c # ConditionalHarmonium f gs (z : zs) x -- ^ Conditional Harmonium
    -> (c # DeepHarmonium gs (z : zs), Dual c #> c # f z x) -- ^ Matrix function and upper part
{-# INLINE splitBottomSubLinear #-}
splitBottomSubLinear dhrm =
    let (dhrmcs,fcs) = S.splitAt $ coordinates dhrm
     in (Point dhrmcs,Point fcs)

-- | Splits the top layer off of a harmonium.
joinBottomSubLinear
    :: (Manifold x, Manifold (f z x), Manifold (DeepHarmonium gs (z : zs)))
    => c # DeepHarmonium gs (z : zs) -- ^ Matrix and upper part
    -> Dual c #> c # f z x
    -> Dual c #> c # ConditionalHarmonium f gs (z : zs) x -- ^ Conditional Harmonium
{-# INLINE joinBottomSubLinear #-}
joinBottomSubLinear (Point dcs) (Point fcs) = Point $ dcs S.++ fcs

-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
mixtureStochasticConditionalCrossEntropy
    :: ( Enum e, ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k, AbsolutelyContinuous Natural z )
    => Sample x -- ^ Input sample
    -> Sample z -- ^ Output sample
    -> Mean #> Natural # MixtureGLM z e k x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE mixtureStochasticConditionalCrossEntropy #-}
mixtureStochasticConditionalCrossEntropy xs ys f =
    let nys = f >$>* xs
     in average $ negate . log <$> zipWith mixtureDensity nys ys


--- Instances ---


instance (Manifold (f n m), Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns)))
    => Manifold (SubLinear f (DeepHarmonium gs (n : ns)) m) where
        type Dimension (SubLinear f (DeepHarmonium gs (n : ns)) m)
          = Dimension (DeepHarmonium gs (n : ns)) + Dimension (f n m)

instance ( Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns)) )
     => Map Mean Natural (SubLinear f) (DeepHarmonium gs (n : ns)) m where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitBottomSubLinear pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitBottomSubLinear pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

{-# LANGUAGE
    KindSignatures,
    DataKinds,
    TypeOperators,
    FlexibleInstances,
    TypeFamilies,
    MultiParamTypeClasses,
    FlexibleContexts,
    ScopedTypeVariables,
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
    , conditionalMixtureRelativeEntropyUpperBound
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

type MixtureGLM z k x =
    ConditionalHarmonium Tensor '[Tensor] [z, Categorical Int k] x -- ^ Function

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
    :: ( ExponentialFamily z, ExponentialFamily x
       , Legendre Natural z, KnownNat k, AbsolutelyContinuous Natural z )
    => Sample x -- ^ Input sample
    -> Sample z -- ^ Output sample
    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE mixtureStochasticConditionalCrossEntropy #-}
mixtureStochasticConditionalCrossEntropy xs ys f =
    let nys = f >$>* xs
     in negate . average $ zipWith logMixtureDensity nys ys

---- | The stochastic cross-entropy of one distribution relative to another, and conditioned
---- on some third variable.
--mixtureStochasticConditionalCrossEntropy
--    :: forall z x k . ( ExponentialFamily z, ExponentialFamily x
--       , Legendre Natural z, KnownNat k, AbsolutelyContinuous Natural z )
--    => Sample x -- ^ Input sample
--    -> Sample z -- ^ Output sample
--    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
--    -> Double -- ^ conditional cross entropy estimate
--{-# INLINE mixtureStochasticConditionalCrossEntropy #-}
--mixtureStochasticConditionalCrossEntropy xs zs f =
--    let (hrm,nzx) = splitBottomSubLinear f
--        (nzc,nc0) = splitBottomHarmonium hrm
--        nzs = nzx >$>* xs
--        nc = fromOneHarmonium nc0
--     in negate . average $ do
--         (z,ncs' <- zip3 zs
--         return $ sum [ log (baseMeasure (Proxy @ z) z)
--                      , sufficientStatistic ox <.> no
--                      , potential (nl <+> ox *<.< nlo)
--                      , negate $ potential (nl <+> rprms) + rho0 ]
--

     --in average $ negate <$> zipWith logMixtureDensity nys ys
--let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
--    let (f,nl0) = splitBottomHarmonium hrm
--        (no,nlo) = splitAffine f
--        nl = fromOneHarmonium nl0
--     in log (baseMeasure (Proxy @ z) ox) + sum
--            [ sufficientStatistic ox <.> no
--            , potential (nl <+> ox *<.< nlo)
--            , negate $ potential (nl <+> rprms) + rho0 ]


-- | The stochastic cross entropy differential of a mixture model.
conditionalMixtureRelativeEntropyUpperBound
    :: forall k z x . ( ClosedFormExponentialFamily z, KnownNat k, ExponentialFamily x )
    => Sample x -- ^ Categorical harmonium
    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
    -> Mean #> Natural # MixtureGLM z k x -- ^ Function
    -> Double -- ^ Upper bound
{-# INLINE conditionalMixtureRelativeEntropyUpperBound #-}
conditionalMixtureRelativeEntropyUpperBound xsmps plkl qlkl =
    average $ zipWith mixtureRelativeEntropyUpperBound (plkl >$>* xsmps) (qlkl >$>* xsmps)

--- Instances ---


instance (Map Mean Natural f y x, Manifold (DeepHarmonium gs (y : zs)))
    => Manifold (SubLinear f (DeepHarmonium gs (y : zs)) x) where
        type Dimension (SubLinear f (DeepHarmonium gs (y : zs)) x)
          = Dimension (DeepHarmonium gs (y : zs)) + Dimension (f y x)

instance ( Map Mean Natural f y x, Manifold (DeepHarmonium gs (y : zs)) )
     => Map Mean Natural (SubLinear f) (DeepHarmonium gs (y : zs)) x where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitBottomSubLinear pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitBottomSubLinear pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

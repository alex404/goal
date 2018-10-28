{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.ExponentialFamily.Harmonium.Conditional where


--- Types ---

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Rectification

import qualified Goal.Core.Vector.Storable as S
import Goal.Probability.ExponentialFamily.Harmonium

data SubLinear (f :: * -> * -> *) h m

type ConditionalHarmonium f m gs ns = SubLinear f (DeepHarmonium gs ns) m

type MixtureGLM m n e k =
    ConditionalHarmonium Tensor m '[Tensor] [n, Categorical e k] -- ^ Function

-- | Splits the top layer off of a harmonium.
splitBottomSubLinear
    :: (Manifold m, Manifold (f n m), Manifold (DeepHarmonium gs (n : ns)))
    => Dual c #> c # ConditionalHarmonium f m gs (n : ns) -- ^ Conditional Harmonium
    -> (Dual c #> c # f n m, c # DeepHarmonium gs (n : ns)) -- ^ Matrix function and upper part
{-# INLINE splitBottomSubLinear #-}
splitBottomSubLinear dhrm =
    let (mtxcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point mtxcs, Point dcs)

-- | Splits the top layer off of a harmonium.
joinBottomSubLinear
    :: (Manifold m, Manifold (f n m), Manifold (DeepHarmonium gs (n : ns)))
    => Dual c #> c # f n m
    -> c # DeepHarmonium gs (n : ns) -- ^ Matrix and upper part
    -> Dual c #> c # ConditionalHarmonium f m gs (n : ns) -- ^ Conditional Harmonium
{-# INLINE joinBottomSubLinear #-}
joinBottomSubLinear (Point mtxcs) (Point dcs) = Point $ mtxcs S.++ dcs

-- | The stochastic cross-entropy of one distribution relative to another, and conditioned
-- on some third variable.
mixtureStochasticConditionalCrossEntropy
    :: ( Enum e, ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k, 1 <= k )
    => Sample x -- ^ Input sample
    -> Sample z -- ^ Output sample
    -> Mean #> Natural # MixtureGLM x z e k -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE mixtureStochasticConditionalCrossEntropy #-}
mixtureStochasticConditionalCrossEntropy xs ys f =
    let nys = f >$>* xs
     in average $ negate . log <$> zipWith mixtureDensity nys ys

-- | The stochastic conditional cross-entropy differential, based on target
-- inputs and outputs expressed as distributions in mean coordinates (this is
-- primarily of internal use).
mixtureStochasticConditionalCrossEntropyDifferential
    :: ( Enum e, ExponentialFamily z, ExponentialFamily x, Legendre Natural z, KnownNat k, 1 <= k )
    => Sample x -- ^ Input mean distributions
    -> Sample z -- ^ Output mean distributions
    -> Mean #> Natural # MixtureGLM x z e k -- ^ Function
    -> CotangentVector (Mean #> Natural) (MixtureGLM x z e k) -- ^ Differential
{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential #-}
mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm =
    -- This could be better optimized but not throwing out the second result of propagate
    let dmglms = dualIsomorphism
            <$> zipWith stochasticMixtureModelDifferential ((:[]) <$> zs) (mglm >$>* xs)
        dzs = [ fst . splitAffine . fst $ splitBottomHarmonium dmglm | dmglm <- dmglms ]
        f = fst $ splitBottomSubLinear mglm
        df = fst $ propagate dzs (sufficientStatistic <$> xs) f
     in primalIsomorphism . joinBottomSubLinear df $ averagePoint dmglms

-- | A gradient for rectifying gains which won't allow them to be negative.
conditionalHarmoniumRectificationDifferential
    :: ( ExponentialFamily x, Manifold (f z x), Map Mean Natural f z x, Manifold (DeepHarmonium gs (z : zs))
       , Legendre Natural (DeepHarmonium gs (z : zs)) )
    => Double -- ^ Rectification shift
    -> Natural # x -- ^ Rectification parameters
    -> Sample x -- ^ Sample points
    -> Mean #> Natural # f z x -- ^ linear part of ppc
    -> Natural # DeepHarmonium gs (z : zs) -- ^ Gains
    -> CotangentPair Natural (DeepHarmonium gs (z : zs)) -- ^ Rectified PPC
{-# INLINE conditionalHarmoniumRectificationDifferential #-}
conditionalHarmoniumRectificationDifferential rho0 rprms xsmps tns dhrm =
    let lkl = joinBottomSubLinear tns dhrm
        rcts = rectificationCurve rho0 rprms xsmps
        ndhrmlkls = lkl >$>* xsmps
        mdhrmlkls = dualTransition <$> ndhrmlkls
        ptns = potential <$> ndhrmlkls
     in joinTangentPair dhrm . averagePoint
         $ [ primalIsomorphism $ (ptn - rct) .> mdhrmlkl | (rct,mdhrmlkl,ptn) <- zip3 rcts mdhrmlkls ptns ]


--- Instances ---


instance (Manifold (f n m), Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns)))
    => Manifold (SubLinear f (DeepHarmonium gs (n : ns)) m) where
        type Dimension (SubLinear f (DeepHarmonium gs (n : ns)) m)
          = Dimension (DeepHarmonium gs (n : ns)) + Dimension (f n m)

instance ( Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns)) )
     => Map Mean Natural (SubLinear f) (DeepHarmonium gs (n : ns)) m where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (pq,dhrm) = splitBottomSubLinear pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (pq,dhrm) = splitBottomSubLinear pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

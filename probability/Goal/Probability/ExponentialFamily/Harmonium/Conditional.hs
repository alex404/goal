{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.ExponentialFamily.Harmonium.Conditional where


--- Types ---

import Goal.Core
import Goal.Geometry
import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily

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

---- | The stochastic cross-entropy of one distribution relative to another, and conditioned
---- on some third variable.
--mixtureStochasticConditionalCrossEntropy
--    :: (ExponentialFamily n, ExponentialFamily m, Legendre Natural m, 1 <= k)
--    => Sample n -- ^ Input sample
--    -> Sample m -- ^ Output sample
--    -> Mean #> Natural # MixtureGLM m n e k -- ^ Function
--    -> Double -- ^ conditional cross entropy estimate
--{-# INLINE mixtureStochasticConditionalCrossEntropy #-}
--mixtureStochasticConditionalCrossEntropy xs ys f =
--    let ykshts = f >$>* xs
--    average . zipWith stochasticCrossEntropy ((:[]) <$> ys) $ f >$>* xs
--
---- | The stochastic conditional cross-entropy differential, based on target
---- inputs and outputs expressed as distributions in mean coordinates (this is
---- primarily of internal use).
--mixtureStochasticConditionalCrossEntropyDifferential0
--    :: (Propagate Mean Natural f m n, ExponentialFamily n, Legendre Natural m)
--    => [Mean # n] -- ^ Input mean distributions
--    -> [Mean # m] -- ^ Output mean distributions
--    -> Mean #> Natural # f m n -- ^ Function
--    -> CotangentVector (Mean #> Natural) (f m n) -- ^ Differential
--{-# INLINE mixtureStochasticConditionalCrossEntropyDifferential0 #-}
--mixtureStochasticConditionalCrossEntropyDifferential0 xs ys f =
--    let (df,yhts) = propagate mys xs f
--        mys = dualIsomorphism <$> zipWith (<->) (potentialDifferential <$> yhts) (primalIsomorphism <$> ys)
--     in primalIsomorphism df


--- Instances ---


instance (Manifold (f n m), Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns)))
    => Manifold (SubLinear f (DeepHarmonium gs (n : ns)) m) where
        type Dimension (SubLinear f (DeepHarmonium gs (n : ns)) m)
          = Dimension (DeepHarmonium gs (n : ns)) + Dimension (f n m)

instance ( Map Mean Natural f n m, Manifold (DeepHarmonium gs (n : ns))
         , Dimension n <= Dimension (DeepHarmonium gs (n : ns)) )
     => Map Mean Natural (SubLinear f) (DeepHarmonium gs (n : ns)) m where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (pq,dhrm) = splitBottomSubLinear pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (pq,dhrm) = splitBottomSubLinear pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

module Goal.Probability.ExponentialFamily.Harmonium.Rectification where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S


--- Functions ---


categoricalHarmoniumRectificationParameters
    :: (KnownNat k, 1 <= k, Enum e, ClosedFormExponentialFamily z)
    => Point Natural (Harmonium Tensor z (Categorical e k))
    -> (Double, Point Natural (Categorical e k))
{-# INLINE categoricalHarmoniumRectificationParameters #-}
categoricalHarmoniumRectificationParameters hrm =
    let (nz,nxz,_) = splitHeadHarmonium hrm
        rho0 = potential nz
        rprms = S.map (\nxzi -> subtract rho0 . potential $ nz <+> Point nxzi) $ S.toRows (toMatrix nxz)
     in (rho0, Point rprms)

sampleRectifiedHarmonium
    :: ( Bilinear f m n, ExponentialFamily m, Generative Natural m, Generative Natural n
       , Hierarchical '[f] '[n,m], KnownNat k )
    => Natural # m
    -> Natural # Harmonium f m n
    -> Random s (Sample k (Harmonium f m n))
{-# INLINE sampleRectifiedHarmonium #-}
sampleRectifiedHarmonium rprms hrm = do
    let (lx,_,_) = splitHeadHarmonium hrm
    xs <- sample (lx <+> rprms)
    downwardPass hrm xs

--sampleCategoricalHarmonium
--    :: ( KnownNat k, Enum e, KnownNat n, 1 <= n
--       , ClosedFormExponentialFamily o, Generative Natural o, Hierarchical '[Tensor] '[Categorical e n, o] )
--      => Point Natural (Categorical e n <*> o)
--      -> Random s (Sample k (Categorical e n <*> o))
--{-# INLINE sampleCategoricalHarmonium #-}
--sampleCategoricalHarmonium hrm = do
--    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
--    sampleRectifiedHarmonium rx hrm
--
--estimateRectifiedHarmoniumDifferentials
--    :: ( Bilinear f m n, ExponentialFamily (Harmonium f m n), KnownNat k, Hierarchical '[f] [m,n]
--       , ExponentialFamily m, ExponentialFamily n, Generative Natural m, Generative Natural n, 1 <= k )
--      => Sample k n
--      -> Natural # m -- ^ Rectification Parameters
--      -> Natural # Harmonium f m n
--      -> Random s (CotangentVector Natural (Harmonium f m n))
--{-# INLINE estimateRectifiedHarmoniumDifferentials #-}
--estimateRectifiedHarmoniumDifferentials zs rprms hrm = do
--    pxzs <- upwardPass hrm zs
--    qxzs <- sampleRectifiedHarmonium rprms hrm
--    return $ estimateStochasticCrossEntropyDifferential pxzs qxzs
--
--estimateCategoricalHarmoniumDifferentials
--    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Hierarchical '[Tensor] '[Categorical e n, o]
--       , ClosedFormExponentialFamily o, Generative Natural o, KnownNat n )
--      => Sample k o
--      -> Point Natural (Categorical e n <*> o)
--      -> Random s (CotangentVector Natural (Categorical e n <*> o))
--{-# INLINE estimateCategoricalHarmoniumDifferentials #-}
--estimateCategoricalHarmoniumDifferentials zs hrm = do
--    pxzs <- upwardPass hrm zs
--    qxzs <- sampleCategoricalHarmonium hrm
--    return $ estimateStochasticCrossEntropyDifferential pxzs qxzs
--
--rectifiedHarmoniumNegativeLogLikelihood
--    :: ( Bilinear f m n, ExponentialFamily (Harmonium f m n), Map Mean Natural f m n
--       , ClosedFormExponentialFamily m, ClosedFormExponentialFamily n )
--      => (Double, Natural # m) -- ^ Rectification Parameters
--      -> Natural # Harmonium f m n
--      -> SamplePoint n
--      -> Double
--{-# INLINE rectifiedHarmoniumNegativeLogLikelihood #-}
--rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
--    let (nl,nlo,no) = splitHeadHarmonium hrm
--     in negate $ sufficientStatistic ox <.> fromOneHarmonium no + potential (nl <+> nlo >.>* ox) - potential (nl <+> rprms) - rho0
--
--categoricalHarmoniumNegativeLogLikelihood
--    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o )
--    => Point Natural (Categorical e k <*> o)
--    -> SamplePoint o
--    -> Double
--{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
--categoricalHarmoniumNegativeLogLikelihood hrm =
--    rectifiedHarmoniumNegativeLogLikelihood (categoricalHarmoniumRectificationParameters hrm) hrm
----
----
----{-
----rectifiedHarmoniumRelativeEntropyDifferentials
----    :: (ExponentialFamily x, SourceGenerative Natural x, ClosedFormExponentialFamily z, SourceGenerative Natural z)
----    => [Sample z]
----    -> Natural :#: x
----    -> Natural :#: Harmonium x z
----    -> RandST s (Differentials :#: Tangent Natural (Harmonium x z))
----rectifiedHarmoniumRelativeEntropyDifferentials zs' rx hrm = do
----
----    xzs <- sampleRectifiedHarmonium (length zs') rx hrm
----    let Harmonium x z = manifold hrm
----    let (xs,zs) = unzip xzs
----    let mzs' = sufficientStatistic z <$> zs'
----
----    xs' <- mapM standardGenerate $ conditionalLatentDistribution hrm >$> mzs'
----
----    let mxs = sufficientStatistic x <$> xs
----        mxs' = sufficientStatistic x <$> xs'
----        mzs = sufficientStatistic z <$> zs
----
----    return $ harmoniumGradientCalculator mxs mxs' mzs mzs' hrm
----
----
------ | Bayes' rule given a rectified harmonium generative model.
----rectifiedBayesRule
----    :: (ExponentialFamily x, ExponentialFamily z)
----    => Function Mixture Natural :#: Affine z x -- ^ Likelihood
----    -> Natural :#: x
----    -> Sample z -- ^ Observation
----    -> Natural :#: x -- ^ Prior
----    -> Natural :#: x -- ^ Posterior
----rectifiedBayesRule lklhd rprms z p0 =
----    let mtx = matrixTranspose . snd $ splitAffine lklhd
----     in mtx >.>* z <+> p0 <+> rprms
----     -}
----

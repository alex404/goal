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
    let (nz,nxz,_) = splitBottomHarmonium hrm
        rho0 = potential nz
        rprms = S.map (\nxzi -> subtract rho0 . potential $ nz <+> Point nxzi) $ S.toRows (toMatrix nxz)
     in (rho0, Point rprms)


sampleCategoricalHarmonium
    :: ( KnownNat k, Enum e, KnownNat n, 1 <= n, ClosedFormExponentialFamily o
      , Generative Natural o, Hierarchical '[Tensor] '[o, Categorical e n] )
      => Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (Sample k (Harmonium Tensor o (Categorical e n)))
{-# INLINE sampleCategoricalHarmonium #-}
sampleCategoricalHarmonium hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    sampleRectified (toSingletonSum rx) hrm

estimateRectifiedHarmoniumDifferentials
    :: ( Map Mean Natural f m n, Bilinear f m n, ExponentialFamily (Harmonium f n m)
       , KnownNat k, Hierarchical '[f] [n,m] , ExponentialFamily m, ExponentialFamily n
       , Generative Natural m, Generative Natural n, 1 <= k )
      => Sample k n
      -> Natural # m -- ^ Rectification Parameters
      -> Natural # Harmonium f n m
      -> Random s (CotangentVector Natural (Harmonium f n m))
{-# INLINE estimateRectifiedHarmoniumDifferentials #-}
estimateRectifiedHarmoniumDifferentials zs rprms hrm = do
    pzxs <- initialPass hrm zs
    qzxs <- sampleRectified (toSingletonSum rprms) hrm
    return $ estimateStochasticCrossEntropyDifferential pzxs qzxs

estimateCategoricalHarmoniumDifferentials
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Hierarchical '[Tensor] '[o, Categorical e n]
       , ClosedFormExponentialFamily o, Generative Natural o, KnownNat n )
      => Sample k o
      -> Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (CotangentVector Natural (Harmonium Tensor o (Categorical e n)))
{-# INLINE estimateCategoricalHarmoniumDifferentials #-}
estimateCategoricalHarmoniumDifferentials zs hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    estimateRectifiedHarmoniumDifferentials zs rx hrm

rectifiedHarmoniumNegativeLogLikelihood
    :: ( Bilinear f m n, ExponentialFamily (Harmonium f n m), Map Mean Natural f m n
       , ClosedFormExponentialFamily m, ClosedFormExponentialFamily n )
      => (Double, Natural # m) -- ^ Rectification Parameters
      -> Natural # Harmonium f n m
      -> SamplePoint n
      -> Double
{-# INLINE rectifiedHarmoniumNegativeLogLikelihood #-}
rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
    let (no,nlo,nl0) = splitBottomHarmonium hrm
        nl = fromOneHarmonium nl0
     in negate $ sufficientStatistic ox <.> no + potential (nl <+> nlo >.>* ox) - potential (nl <+> rprms) - rho0

categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o )
    => Point Natural (Harmonium Tensor o (Categorical e k))
    -> SamplePoint o
    -> Double
{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
categoricalHarmoniumNegativeLogLikelihood hrm =
    rectifiedHarmoniumNegativeLogLikelihood (categoricalHarmoniumRectificationParameters hrm) hrm

--- Classes ---

class Rectified fs ms where
    sampleRectified
        :: KnownNat l
        => Natural # Sum (Tail ms)
        -> Natural # DeepHarmonium fs ms
        -> Random s (Sample l (DeepHarmonium fs ms))

instance Generative Natural m => Rectified '[] '[m] where
    {-# INLINE sampleRectified #-}
    sampleRectified _ = sample

instance ( Hierarchical fs (m : ms), Bilinear f m n, Manifold (Sum ms), ExponentialFamily m
         , Rectified fs (m : ms), Generative Natural n)
  => Rectified (f : fs) (n : m : ms) where
    {-# INLINE sampleRectified #-}
    sampleRectified rprms dhrm = do
        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleRectified rprms' $ biasBottom rprm dhrm'
        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (ys *<$< pf)
        return . hZip zs $ hZip ys xs


--
--
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

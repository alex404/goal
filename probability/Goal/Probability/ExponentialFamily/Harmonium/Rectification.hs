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
import qualified Goal.Core.Vector.Boxed as B


--- Functions ---


categoricalHarmoniumRectificationParameters
    :: (KnownNat k, 1 <= k, Enum e, ClosedFormExponentialFamily z)
    => Point Natural (Categorical e k <*> z)
    -> (Double, Point Natural (Categorical e k))
{-# INLINE categoricalHarmoniumRectificationParameters #-}
categoricalHarmoniumRectificationParameters hrm =
    let (_,nxz,nz0) = splitHeadHarmonium hrm
        nz = fromOneHarmonium nz0
        rho0 = potential nz
        rprms = S.map (\nxzi -> subtract rho0 . potential $ nz <+> Point nxzi) $ S.toRows (toMatrix nxz)
     in (rho0, Point rprms)

sampleRectifiedHarmonium
    :: ( Bilinear f m n, ExponentialFamily m, Generative Natural m, Generative Natural n, Gibbs '[f] '[m,n], KnownNat k )
    => Natural # m
    -> Natural # Harmonium f m n
    -> Random s (Sample k (Harmonium f m n))
{-# INLINE sampleRectifiedHarmonium #-}
sampleRectifiedHarmonium rprms hrm = do
    let (lx,_,_) = splitHeadHarmonium hrm
    xs <- sample (lx <+> rprms)
    zs <- xs *>|>* hrm
    return $ B.zipWith (:+:) xs zs

--sampleCategoricalHarmonium
--    :: ( KnownNat k, Enum e
--       , KnownNat n, 1 <= n, ClosedFormExponentialFamily o, Generative Natural o )
--      => Point Natural (Categorical e n <*> o)
--      -> Random s (Sample k (Categorical e n <*> o))
--{-# INLINE sampleCategoricalHarmonium #-}
--sampleCategoricalHarmonium hrm = do
--    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
--    sampleRectifiedHarmonium rx hrm
--
--estimateRectifiedHarmoniumDifferentials
--    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f), ExponentialFamily (Codomain f)
--       , ExponentialFamily (Domain f), Generative Natural (Codomain f)
--       , Generative Natural (Domain f), KnownNat k, 1 <= k )
--      => Sample k (Domain f)
--      -> Point Natural (Codomain f) -- ^ Rectification Parameters
--      -> Point Natural (Harmonium f)
--      -> Random s (CotangentVector Natural (Harmonium f))
--{-# INLINE estimateRectifiedHarmoniumDifferentials #-}
--estimateRectifiedHarmoniumDifferentials pzs rprms hrm = do
--    pxs <- samplePoint $ conditionalLatentDistributions hrm pzs
--    let pxzs = B.zip pxs pzs
--    qxzs <- sampleRectifiedHarmonium rprms hrm
--    return $ estimateStochasticCrossEntropyDifferential pxzs qxzs
--
--estimateCategoricalHarmoniumDifferentials
--    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e
--       , ClosedFormExponentialFamily o, Generative Natural o, KnownNat n )
--      => Sample n o
--      -> Point Natural (Categorical e k <*> o)
--      -> Random s (CotangentVector Natural (Categorical e k <*> o))
--{-# INLINE estimateCategoricalHarmoniumDifferentials #-}
--estimateCategoricalHarmoniumDifferentials pos hrm = do
--    pls <- samplePoint $ conditionalLatentDistributions hrm pos
--    let plos = B.zip pls pos
--    qlos <- sampleCategoricalHarmonium hrm
--    return $ estimateStochasticCrossEntropyDifferential plos qlos
--
--rectifiedHarmoniumNegativeLogLikelihood
--    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f)
--       , ClosedFormExponentialFamily (Codomain f), ClosedFormExponentialFamily (Domain f) )
--      => (Double, Point Natural (Codomain f)) -- ^ Rectification Parameters
--      -> Point Natural (Harmonium f)
--      -> SamplePoint (Domain f)
--      -> Double
--{-# INLINE rectifiedHarmoniumNegativeLogLikelihood #-}
--rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
--    let (nl,no,nlo) = splitHarmonium hrm
--     in negate $ sufficientStatistic ox <.> no + potential (nl <+> nlo >.>* ox) - potential (nl <+> rprms) - rho0
--
--categoricalHarmoniumNegativeLogLikelihood
--    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o )
--    => Point Natural (Categorical e k <*> o)
--    -> SamplePoint o
--    -> Double
--{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
--categoricalHarmoniumNegativeLogLikelihood hrm =
--    rectifiedHarmoniumNegativeLogLikelihood (categoricalHarmoniumRectificationParameters hrm) hrm
--
--
--{-
--rectifiedHarmoniumRelativeEntropyDifferentials
--    :: (ExponentialFamily x, SourceGenerative Natural x, ClosedFormExponentialFamily z, SourceGenerative Natural z)
--    => [Sample z]
--    -> Natural :#: x
--    -> Natural :#: Harmonium x z
--    -> RandST s (Differentials :#: Tangent Natural (Harmonium x z))
--rectifiedHarmoniumRelativeEntropyDifferentials zs' rx hrm = do
--
--    xzs <- sampleRectifiedHarmonium (length zs') rx hrm
--    let Harmonium x z = manifold hrm
--    let (xs,zs) = unzip xzs
--    let mzs' = sufficientStatistic z <$> zs'
--
--    xs' <- mapM standardGenerate $ conditionalLatentDistribution hrm >$> mzs'
--
--    let mxs = sufficientStatistic x <$> xs
--        mxs' = sufficientStatistic x <$> xs'
--        mzs = sufficientStatistic z <$> zs
--
--    return $ harmoniumGradientCalculator mxs mxs' mzs mzs' hrm
--
--
---- | Bayes' rule given a rectified harmonium generative model.
--rectifiedBayesRule
--    :: (ExponentialFamily x, ExponentialFamily z)
--    => Function Mixture Natural :#: Affine z x -- ^ Likelihood
--    -> Natural :#: x
--    -> Sample z -- ^ Observation
--    -> Natural :#: x -- ^ Prior
--    -> Natural :#: x -- ^ Posterior
--rectifiedBayesRule lklhd rprms z p0 =
--    let mtx = matrixTranspose . snd $ splitAffine lklhd
--     in mtx >.>* z <+> p0 <+> rprms
--     -}
--

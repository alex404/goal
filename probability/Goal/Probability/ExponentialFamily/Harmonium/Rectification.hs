-- | Tools for working with rectified deep harmoniums. Rectified deep harmoiums
-- are a subset of harmoniums which (approximately) satisfy a particular
-- equation. This affords trivial sampling and inference for the given model.
module Goal.Probability.ExponentialFamily.Harmonium.Rectification
    ( -- * Rectification
      Rectified (sampleRectified)
    , estimateRectifiedHarmoniumDifferentials
    , rectifiedHarmoniumNegativeLogLikelihood
      -- ** Categorical Harmoniums
    , buildCategoricalHarmonium
    , mixtureDensity
    , sampleCategoricalHarmonium
    , estimateCategoricalHarmoniumDifferentials
    , categoricalHarmoniumRectificationParameters
    , categoricalHarmoniumNegativeLogLikelihood
      -- * Miscellaneous
    , harmoniumInformationProjectionDifferential
    ) where


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


-- | Estimates the differential of a rectified harmonium with respect to the
-- relative entropy, and given an observation.
estimateRectifiedHarmoniumDifferentials
    :: ( Map Mean Natural f m n, Bilinear f m n, ExponentialFamily (Harmonium f m n)
       , KnownNat k, Manifold (Harmonium f m n) , ExponentialFamily m, ExponentialFamily n
       , Generative Natural m, Generative Natural n, 1 <= k )
      => Sample k m
      -> Natural # n -- ^ Rectification Parameters
      -> Natural # Harmonium f m n
      -> Random s (CotangentVector Natural (Harmonium f m n))
{-# INLINE estimateRectifiedHarmoniumDifferentials #-}
estimateRectifiedHarmoniumDifferentials zs rprms hrm = do
    pzxs <- initialPass hrm zs
    qzxs <- sampleRectified (toSingletonSum rprms) hrm
    return $ estimateStochasticCrossEntropyDifferential pzxs qzxs

-- | Computes the negative log-likelihood of a sample point of a rectified harmonium.
rectifiedHarmoniumNegativeLogLikelihood
    :: ( Bilinear f m n, ExponentialFamily (Harmonium f m n), Map Mean Natural f m n
       , ClosedFormExponentialFamily m, ClosedFormExponentialFamily n )
      => (Double, Natural # n) -- ^ Rectification Parameters
      -> Natural # Harmonium f m n
      -> SamplePoint m
      -> Double
{-# INLINE rectifiedHarmoniumNegativeLogLikelihood #-}
rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
    let (no,nlo,nl0) = splitBottomHarmonium hrm
        nl = fromOneHarmonium nl0
     in negate $ sufficientStatistic ox <.> no + potential (nl <+> ox *<.< nlo) - potential (nl <+> rprms) - rho0


-- Misc --

-- | The differential of the dual relative entropy.
harmoniumInformationProjectionDifferential
    :: (KnownNat k, 1 <= k, 2 <= k, ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m)
    => Natural # n -- ^ Model Distribution
    -> Sample k n -- ^ Model Samples
    -> Natural # Harmonium f m n -- ^ Harmonium
    -> CotangentVector Natural n -- ^ Differential Estimate
{-# INLINE harmoniumInformationProjectionDifferential #-}
harmoniumInformationProjectionDifferential px xs hrm =
    let (nn,nmn,nm0) = splitBottomHarmonium hrm
        nm = fromOneHarmonium nm0
        mxs0 = sufficientStatistic xs
        mys0 = splitReplicated $ nmn >$> mxs0
        mxs = splitReplicated mxs0
        mys = S.zipWith (\mx my0 -> mx <.> (px <-> nm) - potential (nn <+> my0)) mxs mys0
        ln = fromIntegral $ length xs
        mxht = averagePoint mxs
        myht = S.sum mys / ln
        cvr = (ln - 1) /> S.zipFold (\z0 mx my -> z0 <+> ((my - myht) .> (mx <-> mxht))) zero mxs mys
     in primalIsomorphism cvr


-- Categorical Harmoniums --

mixtureDensity
    :: (KnownNat k, 1 <= k, Enum e, ClosedFormExponentialFamily z, Transition Source Natural z)
    => Natural # Harmonium Tensor z (Categorical e k)
    -> SamplePoint z
    -> Double
mixtureDensity hrm x =
    let (nz,nzx,nx) = splitBottomHarmonium hrm
        wghts = coordinates . toMean $ snd (categoricalHarmoniumRectificationParameters hrm) <+> fromOneHarmonium nx
        dxs0 = S.map (flip density x . (<+> nz) . Point) . S.toColumns $ toMatrix nzx
        dx1 = density nz x * (1 - S.sum wghts)
     in dx1 + S.sum (S.zipWith (*) wghts dxs0)

-- | A convenience function for building a mixture model.
buildCategoricalHarmonium
    :: forall k e z
    . (KnownNat k, 1 <= k, Enum e, ClosedFormExponentialFamily z, Transition Source Natural z)
    => Source # z -- -- ^ Component Bias
    -> S.Vector k (Source # z) -- ^ Mixture components
    -> Mean # Categorical e k
    -> Natural # Harmonium Tensor z (Categorical e k)
{-# INLINE buildCategoricalHarmonium #-}
buildCategoricalHarmonium sz szs mx =
    let nz0 = toNatural sz
        nz' :: S.Vector 1 (Natural # z)
        (nzs0,nz') = S.splitAt $ S.map toNatural szs
        nz'' = S.head nz'
        nz = nz0 <+> nz''
        nzs = S.map (<-> nz'') nzs0
        nzx = fromMatrix . S.fromColumns $ S.map coordinates nzs
        nx' = snd . categoricalHarmoniumRectificationParameters $ joinBottomHarmonium nz nzx zero
        nx = toOneHarmonium $ toNatural mx <-> nx'
     in joinBottomHarmonium nz nzx nx

-- | Computes the rectification parameters of a harmonium with a categorical latent variable.
categoricalHarmoniumRectificationParameters
    :: (KnownNat k, 1 <= k, Enum e, ClosedFormExponentialFamily z)
    => Point Natural (Harmonium Tensor z (Categorical e k))
    -> (Double, Point Natural (Categorical e k))
{-# INLINE categoricalHarmoniumRectificationParameters #-}
categoricalHarmoniumRectificationParameters hrm =
    let (nz,nzx,_) = splitBottomHarmonium hrm
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz <+> Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleCategoricalHarmonium
    :: ( KnownNat k, Enum e, KnownNat n, 1 <= n, ClosedFormExponentialFamily o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
      => Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (Sample k (Harmonium Tensor o (Categorical e n)))
{-# INLINE sampleCategoricalHarmonium #-}
sampleCategoricalHarmonium hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    sampleRectified (toSingletonSum rx) hrm

-- | Estimates the differential of a categorical harmonium with respect to the
-- relative entropy, and given an observation.
estimateCategoricalHarmoniumDifferentials
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Manifold (Harmonium Tensor o (Categorical e n))
       , ClosedFormExponentialFamily o, Generative Natural o, KnownNat n )
      => Sample k o
      -> Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (CotangentVector Natural (Harmonium Tensor o (Categorical e n)))
{-# INLINE estimateCategoricalHarmoniumDifferentials #-}
estimateCategoricalHarmoniumDifferentials zs hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    estimateRectifiedHarmoniumDifferentials zs rx hrm


-- | Computes the negative log-likelihood of a sample point of a categorical harmonium.
categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o )
    => Point Natural (Harmonium Tensor o (Categorical e k))
    -> SamplePoint o
    -> Double
{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
categoricalHarmoniumNegativeLogLikelihood hrm =
    rectifiedHarmoniumNegativeLogLikelihood (categoricalHarmoniumRectificationParameters hrm) hrm


--- Classes ---

-- | A rectified distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class Rectified fs ms where
    sampleRectified
        :: KnownNat l
        => Natural # Sum (Tail ms)
        -> Natural # DeepHarmonium fs ms
        -> Random s (Sample l (DeepHarmonium fs ms))

instance Generative Natural m => Rectified '[] '[m] where
    {-# INLINE sampleRectified #-}
    sampleRectified _ = sample

instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
         , ExponentialFamily n, Rectified fs (n : ms), Generative Natural m
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
  => Rectified (f : fs) (m : n : ms) where
    {-# INLINE sampleRectified #-}
    sampleRectified rprms dhrm = do
        let (pn,pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleRectified rprms' $ biasBottom rprm dhrm'
        zs <- samplePoint $ mapReplicatedPoint (pn <+>) (pf >$>* ys)
        return . hZip zs $ hZip ys xs

instance ( Enum e, KnownNat n, 1 <= n, ClosedFormExponentialFamily o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
  => Generative Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE samplePoint #-}
      samplePoint hrm = do
          (smp :: Sample 1 (Harmonium Tensor o (Categorical e n))) <- sampleCategoricalHarmonium hrm
          return $ B.head smp

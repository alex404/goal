{-# LANGUAGE UndecidableInstances #-}
-- | Inference from and fitting parameters to data in the context of Deep Harmoniums.
module Goal.Probability.ExponentialFamily.Harmonium.Inference
    (
    -- * Inference
      (<|<)
    , (<|<*)
    , rectifiedBayesRule
    -- * Training
    , estimateRectifiedHarmoniumDifferentials
    , estimateCategoricalHarmoniumDifferentials
    , categoricalHarmoniumExpectationMaximization
    -- * Rectification
    , harmoniumInformationProjectionDifferential
     -- * Log-Likelihoods
    , rectifiedHarmoniumNegativeLogLikelihood
    , categoricalHarmoniumNegativeLogLikelihood
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Types ---


-- | The given deep harmonium conditioned on its bottom layer.
(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms))
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Mean # m
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (p <.< snd (splitAffine f)) dhrm'

-- | The given deep harmonium conditioned on a sample from its bottom layer.
-- This can be interpreted as the posterior of the model given an observation of
-- the bottom layer.
(<|<*) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)), ExponentialFamily m
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> SamplePoint m
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (<|<*) #-}
(<|<*) dhrm x = dhrm <|< sufficientStatistic x

rectifiedBayesRule
    :: ( Manifold (DeepHarmonium fs (n : ms)), Bilinear f m n
       , ExponentialFamily m, Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # n -- ^ Rectification Parameters
      -> Mean ~> Natural # Affine f m n -- ^ Likelihood
      -> SamplePoint m -- ^ Observation
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Prior
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Updated prior
{-# INLINE rectifiedBayesRule #-}
rectifiedBayesRule rprms lkl x dhrm =
    let dhrm' = joinBottomHarmonium lkl $ biasBottom ((-1) .> rprms) dhrm
     in dhrm' <|<* x

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
       , Legendre Natural m, Legendre Natural n, ExponentialFamily m, ExponentialFamily n )
      => (Double, Natural # n) -- ^ Rectification Parameters
      -> Natural # Harmonium f m n
      -> SamplePoint m
      -> Double
{-# INLINE rectifiedHarmoniumNegativeLogLikelihood #-}
rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
    let (f,nl0) = splitBottomHarmonium hrm
        (no,nlo) = splitAffine f
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
    let (affmn,nm0) = splitBottomHarmonium hrm
        (nn,nmn) = splitAffine affmn
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

-- | EM implementation for categorical harmoniums.
categoricalHarmoniumExpectationMaximization
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Manifold (Harmonium Tensor z (Categorical e n))
       , Legendre Natural z, Generative Natural z, KnownNat n, ExponentialFamily z
       , Transition Mean Natural z )
      => Sample k z
      -> Point Natural (Harmonium Tensor z (Categorical e n))
      -> Point Natural (Harmonium Tensor z (Categorical e n))
{-# INLINE categoricalHarmoniumExpectationMaximization #-}
categoricalHarmoniumExpectationMaximization zs hrm =
    let aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        muss = splitReplicated $ aff >$>* zs
        mus = averagePoint muss
        szs = splitReplicated $ sufficientStatistic zs
        (cmpnts0,nrms) = S.zipFold folder (S.replicate zero, S.replicate 0) muss szs
        cmpnts = S.zipWith (/>) nrms cmpnts0
     in buildCategoricalHarmonium zero cmpnts mus
    where folder (cmpnts,nrms) (Point cs) sz =
              let ws = cs S.++ S.singleton (1 - S.sum cs)
                  cmpnts' = S.map (.> sz) ws
               in (S.zipWith (<+>) cmpnts cmpnts', S.add nrms ws)

-- | Estimates the differential of a categorical harmonium with respect to the
-- relative entropy, and given an observation.
estimateCategoricalHarmoniumDifferentials
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Manifold (Harmonium Tensor o (Categorical e n))
       , Legendre Natural o, Generative Natural o, KnownNat n, ExponentialFamily o )
      => Sample k o
      -> Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (CotangentVector Natural (Harmonium Tensor o (Categorical e n)))
{-# INLINE estimateCategoricalHarmoniumDifferentials #-}
estimateCategoricalHarmoniumDifferentials zs hrm = do
    let rx = snd . categoricalLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
    estimateRectifiedHarmoniumDifferentials zs rx hrm


-- | Computes the negative log-likelihood of a sample point of a categorical harmonium.
categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, Legendre Natural o, ExponentialFamily o )
    => Point Natural (Harmonium Tensor o (Categorical e k))
    -> SamplePoint o
    -> Double
{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
categoricalHarmoniumNegativeLogLikelihood hrm =
    let rh0rx = categoricalLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
     in rectifiedHarmoniumNegativeLogLikelihood rh0rx hrm

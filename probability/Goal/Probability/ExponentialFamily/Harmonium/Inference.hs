{-# LANGUAGE UndecidableInstances #-}
-- | Inference and parameter fitting for Harmoniums.
module Goal.Probability.ExponentialFamily.Harmonium.Inference
    (
    -- * Inference
      (<|<)
    , (<|<*)
    , rectifiedBayesRule
    -- * Training
    , estimateRectifiedHarmoniumDifferentials
    , empiricalHarmoniumExpectations
    , stochasticCategoricalHarmoniumDifferentials
    , categoricalHarmoniumExpectationMaximization
    , deepCategoricalHarmoniumExpectationStep
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


-- | The given deep harmonium conditioned on a mean distribution over the bottom layer.
(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms))
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> Mean # m -- ^ Input means
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Conditioned deep harmonium
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (p <.< snd (splitAffine f)) dhrm'

-- | The given deep harmonium conditioned on a sample from its bottom layer.
-- This can be interpreted as the posterior of the model given an observation of
-- the bottom layer.
(<|<*) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)), ExponentialFamily m
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> SamplePoint m -- ^ Observations
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Posterior
{-# INLINE (<|<*) #-}
(<|<*) dhrm x = dhrm <|< sufficientStatistic x

-- | The posterior distribution given a prior and likelihood, where the
-- likelihood is rectified.
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

-- | Estimates the stochastic cross entropy differential of a rectified harmonium with
-- respect to the relative entropy, and given an observation.
estimateRectifiedHarmoniumDifferentials
    :: ( Map Mean Natural f m n, Bilinear f m n, ExponentialFamily (Harmonium f m n)
       , KnownNat k, Manifold (Harmonium f m n) , ExponentialFamily m, ExponentialFamily n
       , Generative Natural m, Generative Natural n, 1 <= k )
      => Sample k m -- ^ Observations
      -> Natural # n -- ^ Rectification Parameters
      -> Natural # Harmonium f m n -- ^ Harmonium
      -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ Differentials
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

-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium.
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

-- | The stochastic cross entropy differential of a categorical harmonium.
stochasticCategoricalHarmoniumDifferentials
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Manifold (Harmonium Tensor o (Categorical e n))
       , Legendre Natural o, Generative Natural o, KnownNat n, ExponentialFamily o )
      => Sample k o -- ^ Observations
      -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
      -> CotangentVector Natural (Harmonium Tensor o (Categorical e n)) -- ^ Differentials
{-# INLINE stochasticCategoricalHarmoniumDifferentials #-}
stochasticCategoricalHarmoniumDifferentials zs hrm =
    let pxs = empiricalHarmoniumExpectations zs hrm
        qxs = dualTransition hrm
     in primalIsomorphism $ qxs <-> pxs

empiricalHarmoniumExpectations
    :: ( KnownNat k, 1 <= k, ExponentialFamily m, Bilinear f n m
       , Bilinear f m n, Map Mean Natural f n m, Legendre Natural n)
    => Sample k m -- ^ Model Samples
    -> Natural # Harmonium f m n -- ^ Harmonium
    -> Mean # Harmonium f m n -- ^ Harmonium expected sufficient statistics
{-# INLINE empiricalHarmoniumExpectations #-}
empiricalHarmoniumExpectations zs hrm =
    let mzs = splitReplicated $ sufficientStatistic zs
        aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        mxs = S.map dualTransition . splitReplicated $ aff >$>* zs
        mzx = averagePoint $ S.zipWith (>.<) mzs mxs
        maff = joinAffine (averagePoint mzs) mzx
     in joinBottomHarmonium maff . toOneHarmonium $ averagePoint mxs

-- | EM implementation for categorical harmoniums.
categoricalHarmoniumExpectationMaximization
    :: ( KnownNat k, 1 <= k, 1 <= n, Enum e, Manifold (Harmonium Tensor z (Categorical e n))
       , Legendre Natural z, KnownNat n, ExponentialFamily z, Transition Mean Natural z )
      => Sample k z -- ^ Observations
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Current Harmonium
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Updated Harmonium
{-# INLINE categoricalHarmoniumExpectationMaximization #-}
categoricalHarmoniumExpectationMaximization zs hrm =
    let zs' = hSingleton <$> zs
        (cats,mzs) = deepCategoricalHarmoniumExpectationStep zs' $ transposeHarmonium hrm
     in buildCategoricalHarmonium (S.map (toNatural . fromOneHarmonium) mzs) cats

-- | EM implementation for categorical harmoniums.
deepCategoricalHarmoniumExpectationStep
    :: ( KnownNat k, 1 <= k, KnownNat n, 1 <= n, Enum e, ExponentialFamily x
       , ExponentialFamily (DeepHarmonium fs (x : zs)) )
      => Sample k (DeepHarmonium fs (x ': zs)) -- ^ Observations
      -> Natural # DeepHarmonium (Tensor ': fs) (Categorical e n ': x ': zs) -- ^ Current Harmonium
      -> (Natural # Categorical e n, S.Vector n (Mean # DeepHarmonium fs (x ': zs)))
{-# INLINE deepCategoricalHarmoniumExpectationStep #-}
deepCategoricalHarmoniumExpectationStep xzs dhrm =
    let aff = fst $ splitBottomHarmonium dhrm
        muss = splitReplicated . toMean $ aff >$>* fmap hHead xzs
        sxzs = splitReplicated $ sufficientStatistic xzs
        (cmpnts0,nrms) = S.zipFold folder (S.replicate zero, S.replicate 0) muss sxzs
     in (toNatural $ averagePoint muss, S.zipWith (/>) nrms cmpnts0)
    where folder (cmpnts,nrms) (Point cs) sxz =
              let ws = cs S.++ S.singleton (1 - S.sum cs)
                  cmpnts' = S.map (.> sxz) ws
               in (S.zipWith (<+>) cmpnts cmpnts', S.add nrms ws)

-- | Computes the negative log-likelihood of a sample point of a categorical harmonium.
categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE categoricalHarmoniumNegativeLogLikelihood #-}
categoricalHarmoniumNegativeLogLikelihood hrm =
    let rh0rx = categoricalLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
     in rectifiedHarmoniumNegativeLogLikelihood rh0rx hrm

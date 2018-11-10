{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums and Rectification.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , DeepHarmonium
    , unnormalizedHarmoniumObservableDensity
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    , fromOneMixture
    , toOneMixture
    -- ** Construction
    , biasBottom
    , getBottomBias
    , splitBottomHarmonium
    , joinBottomHarmonium
    -- ** Transposition
    , TransposeHarmonium (transposeHarmonium)
    -- ** Sampling
    , Gibbs (upwardPass,initialPass)
    , gibbsPass
    -- ** Inference
    , (<|<)
    , (<|<*)
    , empiricalHarmoniumExpectations
    , harmoniumInformationProjectionDifferential
    -- * Rectified Harmoniums
    , stochasticRectifiedHarmoniumDifferential
    , marginalizeRectifiedHarmonium
    , SampleRectified (sampleRectifiedHarmonium)
    -- * Mixture Models
    , buildMixtureModel
    , splitMixtureModel
    , mixtureDensity
    , mixtureDensity0
    -- ** Statistics
    , stochasticMixtureModelDifferential
    , mixtureExpectationMaximization
    , deepMixtureModelExpectationStep
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Rectification
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Internal as I
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H


--- Types ---


-- | A hierarchical generative model defined by exponential families. Note that
-- the first elements of ms is the bottom layer of the hierachy, and each
-- subsequent element is the next layer up.
data DeepHarmonium (fs :: [* -> * -> *]) (ms :: [*])

-- | A trivial 1-layer harmonium.
type OneHarmonium m = DeepHarmonium '[] '[m]

-- | A 2-layer harmonium.
type Harmonium f m n = DeepHarmonium '[f] [m,n]


--- Functions ---


-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromOneMixture :: c # Harmonium Tensor z (Categorical e 1) -> c # z
{-# INLINE fromOneMixture #-}
fromOneMixture = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toOneMixture :: c # z -> c # Harmonium Tensor z (Categorical e 1)
{-# INLINE toOneMixture #-}
toOneMixture = breakPoint

-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromOneHarmonium :: c # OneHarmonium m -> c # m
{-# INLINE fromOneHarmonium #-}
fromOneHarmonium = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toOneHarmonium :: c # m -> c # OneHarmonium m
{-# INLINE toOneHarmonium #-}
toOneHarmonium = breakPoint

-- | Adds a layer defined by an affine function to the bottom of a deep harmonium.
joinBottomHarmonium
    :: (Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => Dual c #> c # Affine f m n -- ^ Affine function
    -> c # DeepHarmonium fs (n : ms) -- ^ Upper part of the deep harmonium
    -> c # DeepHarmonium (f : fs) (m : n : ms) -- ^ Combined deep harmonium
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold m, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep Harmonium
    -> (Dual c #> c # Affine f m n, c # DeepHarmonium fs (n : ms)) -- ^ Affine function and upper part
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms)) )
    => c # m -- ^ Bias step
    -> c # DeepHarmonium fs (m : ms) -- ^ Deep Harmonium
    -> c # DeepHarmonium fs (m : ms) -- ^ Biased deep harmonium
{-# INLINE biasBottom #-}
biasBottom pm0 dhrm =
    let dz = natValInt (Proxy :: Proxy (Dimension m))
        (pmcs,css') = VS.splitAt dz . S.fromSized $ coordinates dhrm
        pmcs' = H.add pmcs $ S.fromSized (coordinates pm0)
     in Point . I.Vector $ pmcs' VS.++ css'

-- | Get the bias of the bottom layer of the given 'DeepHarmonium'.
getBottomBias
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms)) )
    => c # DeepHarmonium fs (m : ms) -- ^ Deep harmonium
    -> c # m -- ^ Bottom layer bias
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let dz = natValInt (Proxy :: Proxy (Dimension m))
     in Point . I.Vector . VS.take dz . S.fromSized $ coordinates dhrm


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where

    -- | Given a 'DeepHarmonium' and an element of its sample space, partially
    -- updates the sample by resampling from the bottom to the top layer, but
    -- without updating the bottom layer itself.
    upwardPass
        :: Natural # DeepHarmonium fs ms -- ^ Deep harmonium
        -> Sample (DeepHarmonium fs ms) -- ^ Initial sample
        -> Random s (Sample (DeepHarmonium fs ms)) -- ^ Partial Gibbs resample

    -- | Generates an element of the sample spaec of a deep harmonium based by
    -- starting from a sample point from the bottom layer, and doing a naive
    -- upward sampling. This does not generate a true sample from the deep
    -- harmonium.
    initialPass
        :: Natural # DeepHarmonium fs ms -- Deep harmonium
        -> Sample (Head ms) -- ^ Bottom layer sample
        -> Random s (Sample (DeepHarmonium fs ms)) -- ^ Initial deep harmonium sample

-- | Harmonium transpotion. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
class Manifold (DeepHarmonium fs ms) => TransposeHarmonium fs ms where
    transposeHarmonium :: Primal c => c # DeepHarmonium fs ms -> c # DeepHarmonium (Reverse3 fs) (Reverse ms)


-- | A single pass of Gibbs sampling. Infinite iteration of this function yields
-- a sample from the given 'DeepHarmonium'.
gibbsPass :: ( Manifold (DeepHarmonium fs (n : ms)), Gibbs (f : fs) (m : n : ms)
             , Map Mean Natural f m n, Generative Natural m, ExponentialFamily n )
  => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep Hamonium
  -> Sample (DeepHarmonium (f : fs) (m : n : ms)) -- ^ Initial Sample
  -> s ~> Sample (DeepHarmonium (f : fs) (m : n : ms)) -- ^ Gibbs resample
{-# INLINE gibbsPass #-}
gibbsPass dhrm zyxs = do
    let yxs = snd $ hUnzip zyxs
        ys = fst $ hUnzip yxs
        f = fst $ splitBottomHarmonium dhrm
    zs <- mapM samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip zs yxs


--- Rectification ---


-- | A rectified distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleRectified fs ms where
    -- | A true sample from a rectified harmonium.
    sampleRectifiedHarmonium
        :: Int -- ^ Sample Size
        -> Natural # Sum (Tail ms) -- ^ Rectification parameters
        -> Natural # DeepHarmonium fs ms -- ^ Deep harmonium
        -> Random s (Sample (DeepHarmonium fs ms)) -- ^ Deep harmonium sample

-- | Marginalize the bottom layer out of a deep harmonium.
marginalizeRectifiedHarmonium
    :: ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
       , ExponentialFamily m )
      => Natural # Sum (n : ms) -- ^ Rectification Parameters
      -> Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> (Natural # Sum ms, Natural # DeepHarmonium fs (n : ms)) -- ^ Marginalized deep harmonium
{-# INLINE marginalizeRectifiedHarmonium #-}
marginalizeRectifiedHarmonium rprms dhrm =
    let dhrm' = snd $ splitBottomHarmonium dhrm
        (rprm,rprms') = splitSum rprms
     in (rprms', biasBottom rprm dhrm')

-- Mixture Models --

-- | The observable density of a categorical harmonium.
mixtureDensity0
    :: (KnownNat k, 1 <= k, Enum e, Legendre Natural z, AbsolutelyContinuous Natural z)
    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
    -> SamplePoint z -- ^ Observation
    -> Double -- ^ Probablity density of the observation
{-# INLINE mixtureDensity0 #-}
mixtureDensity0 hrm x =
    let (affzx,nx) = splitBottomHarmonium hrm
        nz = fst $ splitAffine affzx
        wghts = listCoordinates . toMean
            $ snd (mixtureLikelihoodRectificationParameters affzx) <+> fromOneHarmonium nx
        dxs0 = (`density` x) <$> affzx >$>* pointSampleSpace (fromOneHarmonium nx)
        dx1 = density nz x * (1 - sum wghts)
     in dx1 + sum (zipWith (*) wghts dxs0)

-- | A convenience function for building a categorical harmonium/mixture model.
buildMixtureModel
    :: forall k e z . ( KnownNat k, 1 <= k, Enum e, Legendre Natural z )
    => S.Vector k (Natural # z) -- ^ Mixture components
    -> Natural # Categorical e k -- ^ Weights
    -> Natural # Harmonium Tensor z (Categorical e k) -- ^ Mixture Model
{-# INLINE buildMixtureModel #-}
buildMixtureModel nzs0 nx0 =
    let nz0 :: S.Vector 1 (Natural # z)
        (nzs0',nz0) = S.splitAt nzs0
        nz = S.head nz0
        nzs = S.map (<-> nz) nzs0'
        nzx = fromMatrix . S.fromColumns $ S.map coordinates nzs
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodRectificationParameters affzx
        nx = toOneHarmonium $ nx0 <-> rprms
     in joinBottomHarmonium affzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitMixtureModel
    :: forall k e z . ( KnownNat k, 1 <= k, Enum e, Legendre Natural z )
    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
    -> (S.Vector k (Natural # z), Natural # Categorical e k) -- ^ (components, weights)
{-# INLINE splitMixtureModel #-}
splitMixtureModel hrm =
    let (affzx,nx) = splitBottomHarmonium hrm
        rprms = snd $ mixtureLikelihoodRectificationParameters affzx
        nx0 = fromOneHarmonium nx <+> rprms
        (nz,nzx) = splitAffine affzx
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (<+> nz) nzs
        nz0 = S.singleton nz
     in (nzs0' S.++ nz0,nx0)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleMixtureModel
    :: ( Enum e, KnownNat n, 1 <= n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
       => Int
       -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
       -> Random s (Sample (Harmonium Tensor o (Categorical e n))) -- ^ Sample
{-# INLINE sampleMixtureModel #-}
sampleMixtureModel k hrm = do
    let rx = snd . mixtureLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
    sampleRectifiedHarmonium k (toSingletonSum rx) hrm


--- Internal Functions ---


harmoniumBaseMeasure
    :: ExponentialFamily m
    => Proxy m
    -> Proxy (OneHarmonium m)
    -> SamplePoint (OneHarmonium m)
    -> Double
{-# INLINE harmoniumBaseMeasure #-}
harmoniumBaseMeasure prxyl _ (x :+: Null) =
     baseMeasure prxyl x

deepHarmoniumBaseMeasure
    :: (ExponentialFamily m, ExponentialFamily (DeepHarmonium fs ms))
    => Proxy m
    -> Proxy (DeepHarmonium fs ms)
    -> Proxy (DeepHarmonium (f : fs) (m : ms))
    -> SamplePoint (DeepHarmonium (f : fs) (m : ms))
    -> Double
{-# INLINE deepHarmoniumBaseMeasure #-}
deepHarmoniumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

mixtureExpectations
    :: ( Enum e, KnownNat k, 1 <= k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k)
    -> Mean # Harmonium Tensor o (Categorical e k)
{-# INLINE mixtureExpectations #-}
mixtureExpectations hrm =
    let (nzs,nx) = splitMixtureModel hrm
        mzs0 = S.map dualTransition nzs
        mx = dualTransition nx
        pis0 = coordinates mx
        pi' = 1 - S.sum pis0
        pis = pis0 S.++ S.singleton pi'
        mzs0' = S.zipWith (.>) pis mzs0
        mzs = S.take mzs0'
        mz = S.foldr1 (<+>) mzs0'
        mzx = fromMatrix . S.fromColumns $ S.map coordinates mzs
     in joinBottomHarmonium (joinAffine mz mzx) $ toOneHarmonium mx

-- | The given deep harmonium conditioned on a mean distribution over the bottom layer.
(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)) )
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
(<|<*) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)), ExponentialFamily m )
      => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep harmonium
      -> SamplePoint m -- ^ Observations
      -> Natural # DeepHarmonium fs (n : ms) -- ^ Posterior
{-# INLINE (<|<*) #-}
(<|<*) dhrm x = dhrm <|< sufficientStatistic x

-- | Estimates the stochastic cross entropy differential of a rectified harmonium with
-- respect to the relative entropy, and given an observation.
--
-- NB: Right now I'm using length on a list where I probably don't need to...
stochasticRectifiedHarmoniumDifferential
    :: ( Map Mean Natural f m n, Bilinear f m n, ExponentialFamily (Harmonium f m n)
       , Manifold (Harmonium f m n) , ExponentialFamily m, ExponentialFamily n
       , Generative Natural m, Generative Natural n )
       => Sample m -- ^ Observations
       -> Natural # n -- ^ Rectification Parameters
       -> Natural # Harmonium f m n -- ^ Harmonium
       -> Random s (CotangentVector Natural (Harmonium f m n)) -- ^ Differential
{-# INLINE stochasticRectifiedHarmoniumDifferential #-}
stochasticRectifiedHarmoniumDifferential zs rprms hrm = do
    pzxs <- initialPass hrm zs
    qzxs <- sampleRectifiedHarmonium (length zs) (toSingletonSum rprms) hrm
    return $ stochasticCrossEntropyDifferential' pzxs qzxs

-- | Computes the negative log-likelihood of a sample point of a rectified harmonium.
rectifiedHarmoniumDensity
    :: forall f m n . ( Bilinear f m n, ExponentialFamily (Harmonium f m n), Map Mean Natural f m n
       , Legendre Natural m, Legendre Natural n, ExponentialFamily m, ExponentialFamily n )
      => (Double, Natural # n) -- ^ Rectification Parameters
      -> Natural # Harmonium f m n
      -> SamplePoint m
      -> Double
{-# INLINE rectifiedHarmoniumDensity #-}
rectifiedHarmoniumDensity (rho0,rprms) hrm ox =
    let (f,nl0) = splitBottomHarmonium hrm
        (no,nlo) = splitAffine f
        nl = fromOneHarmonium nl0
     in (* baseMeasure (Proxy :: Proxy m) ox) . exp $ sum
            [ sufficientStatistic ox <.> no
            , potential (nl <+> ox *<.< nlo)
            , negate $ potential (nl <+> rprms) + rho0 ]


-- | Computes the negative log-likelihood of a sample point of a mixture model.
mixtureDensity
    :: ( Enum e, KnownNat k, 1 <= k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE mixtureDensity #-}
mixtureDensity hrm =
    let rh0rx = mixtureLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
     in rectifiedHarmoniumDensity rh0rx hrm


-- Misc --

unnormalizedHarmoniumObservableDensity
    :: forall f z x . (ExponentialFamily z, Legendre Natural x, Bilinear f z x)
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableDensity hrm z =
    let (affzx,nx0) = splitBottomHarmonium hrm
        (nz,nzx) = splitAffine affzx
        nx = fromOneHarmonium nx0
        mz = sufficientStatistic z
     in baseMeasure (Proxy :: Proxy z) z * exp (nz <.> mz + potential (nx <+> mz <.< nzx))


-- | The differential of the dual relative entropy. Minimizing this results in
-- the information projection of the model against the marginal distribution of
-- the given harmonium. This is more efficient than the generic version.
harmoniumInformationProjectionDifferential
    :: (ExponentialFamily n, Map Mean Natural f m n, Legendre Natural m)
    => Natural # n -- ^ Model Distribution
    -> Sample n -- ^ Model Samples
    -> Natural # Harmonium f m n -- ^ Harmonium
    -> CotangentVector Natural n -- ^ Differential Estimate
{-# INLINE harmoniumInformationProjectionDifferential #-}
harmoniumInformationProjectionDifferential px xs hrm =
    let (affmn,nm0) = splitBottomHarmonium hrm
        (nn,nmn) = splitAffine affmn
        nm = fromOneHarmonium nm0
        mxs = sufficientStatistic <$> xs
        mys0 = nmn >$> mxs
        mys = zipWith (\mx my0 -> mx <.> (px <-> nm) - potential (nn <+> my0)) mxs mys0
        ln = fromIntegral $ length xs
        mxht = averagePoint mxs
        myht = sum mys / ln
        foldfun (mx,my) (k,z0) = (k+1,z0 <+> ((my - myht) .> (mx <-> mxht)))
        cvr = uncurry (/>) . foldr foldfun (-1,zero) $ zip mxs mys
     in primalIsomorphism cvr

-- | The stochastic cross entropy differential of a mixture model.
stochasticMixtureModelDifferential
    :: ( 1 <= n, Enum e, Manifold (Harmonium Tensor o (Categorical e n))
       , Legendre Natural o, KnownNat n, ExponentialFamily o )
      => Sample o -- ^ Observations
      -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
      -> CotangentVector Natural (Harmonium Tensor o (Categorical e n)) -- ^ Differential
{-# INLINE stochasticMixtureModelDifferential #-}
stochasticMixtureModelDifferential zs hrm =
    let pxs = empiricalHarmoniumExpectations zs hrm
        qxs = dualTransition hrm
     in primalIsomorphism $ qxs <-> pxs

empiricalHarmoniumExpectations
    :: ( ExponentialFamily m, Bilinear f n m
       , Bilinear f m n, Map Mean Natural f n m, Legendre Natural n)
    => Sample m -- ^ Model Samples
    -> Natural # Harmonium f m n -- ^ Harmonium
    -> Mean # Harmonium f m n -- ^ Harmonium expected sufficient statistics
{-# INLINE empiricalHarmoniumExpectations #-}
empiricalHarmoniumExpectations zs hrm =
    let mzs = sufficientStatistic <$> zs
        aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        mxs = dualTransition <$> aff >$>* zs
        mzx = averagePoint $ zipWith (>.<) mzs mxs
        maff = joinAffine (averagePoint mzs) mzx
     in joinBottomHarmonium maff . toOneHarmonium $ averagePoint mxs

-- | EM implementation for mixture models/categorical harmoniums.
mixtureExpectationMaximization
    :: ( 1 <= n, Enum e, Manifold (Harmonium Tensor z (Categorical e n))
       , Legendre Natural z, KnownNat n, ExponentialFamily z, Transition Mean Natural z )
      => Sample z -- ^ Observations
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Current Harmonium
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Updated Harmonium
{-# INLINE mixtureExpectationMaximization #-}
mixtureExpectationMaximization zs hrm =
    let zs' = hSingleton <$> zs
        (cats,mzs) = deepMixtureModelExpectationStep zs' $ transposeHarmonium hrm
     in buildMixtureModel (S.map (toNatural . fromOneHarmonium) mzs) cats

-- | E-step implementation for deep mixture models/categorical harmoniums. Note
-- that for the sake of type signatures, this acts on transposed harmoniums
-- (i.e. the categorical variables are at the bottom of the hierarchy).
deepMixtureModelExpectationStep
    :: ( KnownNat n, 1 <= n, Enum e, ExponentialFamily x
       , ExponentialFamily (DeepHarmonium fs (x : zs)) )
      => Sample (DeepHarmonium fs (x ': zs)) -- ^ Observations
      -> Natural # DeepHarmonium (Tensor ': fs) (Categorical e n ': x ': zs) -- ^ Current Harmonium
      -> (Natural # Categorical e n, S.Vector n (Mean # DeepHarmonium fs (x ': zs)))
{-# INLINE deepMixtureModelExpectationStep #-}
deepMixtureModelExpectationStep xzs dhrm =
    let aff = fst $ splitBottomHarmonium dhrm
        muss = toMean <$> aff >$>* fmap hHead xzs
        sxzs = sufficientStatistic <$> xzs
        (cmpnts0,nrms) = foldr folder (S.replicate zero, S.replicate 0) $ zip muss sxzs
     in (toNatural $ averagePoint muss, S.zipWith (/>) nrms cmpnts0)
    where folder (Point cs,sxz) (cmpnts,nrms) =
              let ws = cs S.++ S.singleton (1 - S.sum cs)
                  cmpnts' = S.map (.> sxz) ws
               in (S.zipWith (<+>) cmpnts cmpnts', S.add nrms ws)


----- Instances ---


instance Manifold m => Manifold (DeepHarmonium fs '[m]) where
    type Dimension (DeepHarmonium fs '[m]) = Dimension m

instance (Manifold m, Manifold n, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
  => Manifold (DeepHarmonium (f : fs) (m : n : ms)) where
      type Dimension (DeepHarmonium (f : fs) (m : n : ms))
        = Dimension m + Dimension (f m n) + Dimension (DeepHarmonium fs (n : ms))

instance Manifold (DeepHarmonium fs ms) => Statistical (DeepHarmonium fs ms) where
    type SamplePoint (DeepHarmonium fs ms) = HList (SamplePoints ms)

instance Generative c m => Generative c (OneHarmonium m) where
    {-# INLINE samplePoint #-}
    samplePoint = fmap (:+: Null) . samplePoint . fromOneHarmonium

instance ExponentialFamily m => ExponentialFamily (OneHarmonium m) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (x :+: Null) =
          toOneHarmonium $ sufficientStatistic x
      {-# INLINE baseMeasure #-}
      baseMeasure = harmoniumBaseMeasure Proxy

instance ( ExponentialFamily n, ExponentialFamily m
         , Bilinear f m n, ExponentialFamily (DeepHarmonium fs (n : ms)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (m : n : ms)) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              pm = sufficientStatistic xm
              pn = sufficientStatistic xn
           in joinBottomHarmonium (joinAffine pm $ pm >.< pn) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n )
  => Gibbs '[f] '[m,n] where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zxs = initialPass dhrm . fst $ hUnzip zxs
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              xp = fromOneHarmonium dhrm'
          xs <- mapM (samplePoint . (<+>) xp) $ zs *<$< f
          return . hZip zs . hZip xs $ repeat Null

instance ( Bilinear f m n, Map Mean Natural g n o, Manifold (DeepHarmonium fs (o : ms))
         , ExponentialFamily m, ExponentialFamily o, Generative Natural n, Gibbs (g : fs) (n : o : ms) )
  => Gibbs (f : g : fs) (m : n : o : ms) where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zyxs = do
          let (zs,yxs) = hUnzip zyxs
              (xs,xs') = hUnzip . snd $ hUnzip yxs
              (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              (g,_) = splitBottomHarmonium dhrm'
          ys <- mapM samplePoint $ zipWith (<+>) (g >$>* xs) (zs *<$< f)
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              yp = fst . splitAffine . fst $ splitBottomHarmonium dhrm'
          ys <- mapM (samplePoint . (<+> yp)) $ zs *<$< f
          yxs' <- initialPass dhrm' ys
          return $ hZip zs yxs'

instance Manifold m => TransposeHarmonium '[] '[m] where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium = id

instance (Bilinear f m n, Bilinear f n m, TransposeHarmonium fs (n : ms))
  => TransposeHarmonium (f : fs) (m : n : ms) where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium dhrm =
        let (aff,dhrm') = splitBottomHarmonium dhrm
            (pm,pmtx) = splitAffine aff
            dhrm'' = transposeHarmonium dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

instance Generative Natural m => SampleRectified '[] '[m] where
    {-# INLINE sampleRectifiedHarmonium #-}
    sampleRectifiedHarmonium k _ = sample k

instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
         , ExponentialFamily n, SampleRectified fs (n : ms), Generative Natural m )
  => SampleRectified (f : fs) (m : n : ms) where
    {-# INLINE sampleRectifiedHarmonium #-}
    sampleRectifiedHarmonium k rprms dhrm = do
        let (pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleRectifiedHarmonium k rprms' $ biasBottom rprm dhrm'
        zs <- mapM samplePoint $ pf >$>* ys
        return . hZip zs $ hZip ys xs

instance ( Enum e, KnownNat n, 1 <= n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
  => Generative Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE samplePoint #-}
      samplePoint hrm = head <$> sampleMixtureModel 1 hrm

instance ( Enum e, KnownNat n, 1 <= n, Legendre Natural o, ExponentialFamily o
  , Manifold (Harmonium Tensor o (Categorical e n)))
  => Legendre Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE potential #-}
      potential hrm =
          let (lkl,nx0) = splitBottomHarmonium hrm
              (rho0,rprms) = mixtureLikelihoodRectificationParameters lkl
              nx = fromOneHarmonium nx0
           in potential (nx <+> rprms) + rho0

           --in logSumExp [ sufficientStatistic i <.> nx + potential (lkl >.>* i) | i <- pointSampleSpace nx ]
      potentialDifferential = primalIsomorphism . mixtureExpectations

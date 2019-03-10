{-# LANGUAGE
    RankNTypes,
    PolyKinds,
    DataKinds,
    TypeOperators,
    MultiParamTypeClasses,
    FlexibleContexts,
    FlexibleInstances,
    TypeFamilies,
    TypeApplications,
    ScopedTypeVariables,
    UndecidableInstances
#-}
-- | Exponential Family Harmoniums and Conjugation.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , DeepHarmonium
    , unnormalizedHarmoniumObservableDensity
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    , fromNullMixture
    , toNullMixture
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
    , harmoniumEmpiricalExpectations
    -- * Conjugated Harmoniums
    , conjugatedHarmoniumDensity
    , logConjugatedHarmoniumDensity
    , marginalizeConjugatedHarmonium
    , SampleConjugated (sampleConjugatedHarmonium)
    -- * Mixture Models
    , joinMixtureModel
    , splitMixtureModel
    , mixtureDensity
    , logMixtureDensity
    , mixtureRelativeEntropyUpperBound
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium.Conjugation
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Internal as I
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H


--- Types ---


-- | A hierarchical generative model defined by exponential families. Note that
-- the first elements of ms is the bottom layer of the hierachy, and each
-- subsequent element is the next layer up.
data DeepHarmonium (fs :: [Type -> Type -> Type]) (xs :: [Type])

-- | A trivial 1-layer harmonium.
type OneHarmonium x = DeepHarmonium '[] '[x]

-- | A 2-layer harmonium.
type Harmonium f z x = DeepHarmonium '[f] [z,x]


--- Functions ---


-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromNullMixture :: c # Harmonium Tensor z (Categorical e 0) -> c # z
{-# INLINE fromNullMixture #-}
fromNullMixture = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toNullMixture :: c # z -> c # Harmonium Tensor z (Categorical e 0)
{-# INLINE toNullMixture #-}
toNullMixture = breakPoint

-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromOneHarmonium :: c # OneHarmonium x -> c # x
{-# INLINE fromOneHarmonium #-}
fromOneHarmonium = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toOneHarmonium :: c # x -> c # OneHarmonium x
{-# INLINE toOneHarmonium #-}
toOneHarmonium = breakPoint

-- | Adds a layer defined by an affine function to the bottom of a deep harmonium.
joinBottomHarmonium
    :: Dual c #> c # Affine f z y -- ^ Affine function
    -> c # DeepHarmonium fs (y : xs) -- ^ Upper part of the deep harmonium
    -> c # DeepHarmonium (f : fs) (z : y : xs) -- ^ Combined deep harmonium
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold z, Manifold (f z y))
    => c # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep Harmonium
    -> (Dual c #> c # Affine f z y, c # DeepHarmonium fs (y : xs)) -- ^ Affine function and upper part
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs z xs c . Manifold z
    => c # z -- ^ Bias step
    -> c # DeepHarmonium fs (z : xs) -- ^ Deep Harmonium
    -> c # DeepHarmonium fs (z : xs) -- ^ Biased deep harmonium
{-# INLINE biasBottom #-}
biasBottom pm0 dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
        (pmcs,css') = VS.splitAt dz . S.fromSized $ coordinates dhrm
        pmcs' = H.add pmcs $ S.fromSized (coordinates pm0)
     in Point . I.Vector $ pmcs' VS.++ css'

-- | Get the bias of the bottom layer of the given 'DeepHarmonium'.
getBottomBias
    :: forall fs z xs c . Manifold z
    => c # DeepHarmonium fs (z : xs) -- ^ Deep harmonium
    -> c # z -- ^ Bottom layer bias
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
     in Point . I.Vector . VS.take dz . S.fromSized $ coordinates dhrm


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [Type -> Type -> Type]) (xs :: [Type]) where

    -- | Given a 'DeepHarmonium' and an element of its sample space, partially
    -- updates the sample by resampling from the bottom to the top layer, but
    -- without updating the bottom layer itself.
    upwardPass
        :: Natural # DeepHarmonium fs xs -- ^ Deep harmonium
        -> Sample (DeepHarmonium fs xs) -- ^ Initial sample
        -> Random s (Sample (DeepHarmonium fs xs)) -- ^ Partial Gibbs resample

    -- | Generates an element of the sample spaec of a deep harmonium based by
    -- starting from a sample point from the bottom layer, and doing a naive
    -- upward sampling. This does not generate a true sample from the deep
    -- harmonium.
    initialPass
        :: Natural # DeepHarmonium fs xs -- Deep harmonium
        -> Sample (Head xs) -- ^ Bottom layer sample
        -> Random s (Sample (DeepHarmonium fs xs)) -- ^ Initial deep harmonium sample

-- | Harmonium transpotion. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
class Manifold (DeepHarmonium fs xs) => TransposeHarmonium fs xs where
    transposeHarmonium :: Primal c => c # DeepHarmonium fs xs -> c # DeepHarmonium (Reverse fs) (Reverse xs)


-- | A single pass of Gibbs sampling. Infinite iteration of this function yields
-- a sample from the given 'DeepHarmonium'.
gibbsPass :: ( Manifold (DeepHarmonium fs (y : xs)), Gibbs (f : fs) (z : y : xs)
             , Map Mean Natural f z y, Generative Natural z, ExponentialFamily y )
  => Natural # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep Hamonium
  -> Sample (DeepHarmonium (f : fs) (z : y : xs)) -- ^ Initial Sample
  -> Random s (Sample (DeepHarmonium (f : fs) (z : y : xs))) -- ^ Gibbs resample
{-# INLINE gibbsPass #-}
gibbsPass dhrm zyxs = do
    let yxs = snd $ hUnzip zyxs
        ys = fst $ hUnzip yxs
        f = fst $ splitBottomHarmonium dhrm
    zs <- mapM samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip zs yxs


--- Conjugation ---


-- | A conjugated distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleConjugated fs xs where
    -- | A true sample from a conjugated harmonium.
    sampleConjugatedHarmonium
        :: Int -- ^ Sample Size
        -> Natural # Sum (Tail xs) -- ^ Conjugation parameters
        -> Natural # DeepHarmonium fs xs -- ^ Deep harmonium
        -> Random s (Sample (DeepHarmonium fs xs)) -- ^ Deep harmonium sample

-- | Marginalize the bottom layer out of a deep harmonium.
marginalizeConjugatedHarmonium
    :: ( Manifold (DeepHarmonium fs (y : xs)), Map Mean Natural f z y, Manifold (Sum xs) )
      => Natural # Sum (y : xs) -- ^ Conjugation Parameters
      -> Natural # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep harmonium
      -> (Natural # Sum xs, Natural # DeepHarmonium fs (y : xs)) -- ^ Marginalized deep harmonium
{-# INLINE marginalizeConjugatedHarmonium #-}
marginalizeConjugatedHarmonium rprms dhrm =
    let dhrm' = snd $ splitBottomHarmonium dhrm
        (rprm,rprms') = splitSum rprms
     in (rprms', biasBottom rprm dhrm')

-- Mixture Models --

-- | The stochastic cross entropy differential of a mixture model.
mixtureRelativeEntropyUpperBound
    :: forall z e n . ( ClosedFormExponentialFamily z, Enum e, KnownNat n )
      => Natural # Harmonium Tensor z (Categorical e n) -- ^ Categorical harmonium
      -> Natural # Harmonium Tensor z (Categorical e n) -- ^ Categorical harmonium
      -> Double -- ^ Upper bound
{-# INLINE mixtureRelativeEntropyUpperBound #-}
mixtureRelativeEntropyUpperBound phrm qhrm =
    let pzc = fst $ splitBottomHarmonium phrm
        npc = snd $ splitMixtureModel phrm
        spc = toSource npc
        qzc = fst $ splitBottomHarmonium qhrm
        qc = snd $ splitMixtureModel qhrm
        wghts = (1 - S.sum (coordinates spc)) : listCoordinates spc
        smps = sampleSpace $ Proxy @ (Categorical e n)
        pzs = pzc >$>* smps
        qzs = qzc >$>* smps
        dvg0 = weightedAverage (zip wghts . zipWith divergence pzs $ dualTransition <$> qzs)
     in divergence npc (transition qc) + dvg0


-- | A convenience function for building a categorical harmonium/mixture model.
joinMixtureModel
    :: forall k e z . ( KnownNat k, Enum e, Legendre Natural z )
    => S.Vector (k+1) (Natural # z) -- ^ Mixture components
    -> Natural # Categorical e k -- ^ Weights
    -> Natural # Harmonium Tensor z (Categorical e k) -- ^ Mixture Model
{-# INLINE joinMixtureModel #-}
joinMixtureModel nzs0 nx0 =
    let nz0 :: S.Vector 1 (Natural # z)
        (nz0,nzs0') = S.splitAt nzs0
        nz = S.head nz0
        nzs = S.map (<-> nz) nzs0'
        nzx = fromMatrix . S.fromColumns $ S.map coordinates nzs
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx = toOneHarmonium $ nx0 <-> rprms
     in joinBottomHarmonium affzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitMixtureModel
    :: forall k e z . ( KnownNat k, Enum e, Legendre Natural z )
    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # z), Natural # Categorical e k) -- ^ (components, weights)
{-# INLINE splitMixtureModel #-}
splitMixtureModel hrm =
    let (affzx,nx) = splitBottomHarmonium hrm
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx0 = fromOneHarmonium nx <+> rprms
        (nz,nzx) = splitAffine affzx
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (<+> nz) nzs
     in (S.cons nz nzs0',nx0)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleMixtureModel
    :: ( Enum e, KnownNat n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
       => Int
       -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
       -> Random s (Sample (Harmonium Tensor o (Categorical e n))) -- ^ Sample
{-# INLINE sampleMixtureModel #-}
sampleMixtureModel k hrm = do
    let rx = snd . mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
    sampleConjugatedHarmonium k (toSingletonSum rx) hrm


--- Internal Functions ---


harmoniumBaseMeasure
    :: ExponentialFamily x
    => Proxy x
    -> Proxy (OneHarmonium x)
    -> SamplePoint (OneHarmonium x)
    -> Double
{-# INLINE harmoniumBaseMeasure #-}
harmoniumBaseMeasure prxyl _ (x :+: Null) =
     baseMeasure prxyl x

deepHarmoniumBaseMeasure
    :: (ExponentialFamily x, ExponentialFamily (DeepHarmonium fs xs))
    => Proxy x
    -> Proxy (DeepHarmonium fs xs)
    -> Proxy (DeepHarmonium (f : fs) (x : xs))
    -> SamplePoint (DeepHarmonium (f : fs) (x : xs))
    -> Double
{-# INLINE deepHarmoniumBaseMeasure #-}
deepHarmoniumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

mixtureExpectations
    :: ( Enum e, KnownNat k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k)
    -> Mean # Harmonium Tensor o (Categorical e k)
{-# INLINE mixtureExpectations #-}
mixtureExpectations hrm =
    let (nzs,nx) = splitMixtureModel hrm
        mzs0 = S.map dualTransition nzs
        mx = dualTransition nx
        pis0 = coordinates mx
        pi' = 1 - S.sum pis0
        pis = S.cons pi' pis0
        mzs0' = S.zipWith (.>) pis mzs0
        mzs = S.tail mzs0'
        mz = S.foldr1 (<+>) mzs0'
        mzx = fromMatrix . S.fromColumns $ S.map coordinates mzs
     in joinBottomHarmonium (joinAffine mz mzx) $ toOneHarmonium mx

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
conjugatedHarmoniumDensity
    :: forall f z x . ( Bilinear f z x, ExponentialFamily (Harmonium f z x), Map Mean Natural f z x
       , Legendre Natural z, Legendre Natural x, ExponentialFamily z, ExponentialFamily x )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f z x
      -> SamplePoint z
      -> Double
{-# INLINE conjugatedHarmoniumDensity #-}
conjugatedHarmoniumDensity (rho0,rprms) hrm ox =
    let (f,nl0) = splitBottomHarmonium hrm
        (no,nlo) = splitAffine f
        nl = fromOneHarmonium nl0
     in (* baseMeasure (Proxy @ z) ox) . exp $ sum
            [ sufficientStatistic ox <.> no
            , potential (nl <+> ox *<.< nlo)
            , negate $ potential (nl <+> rprms) + rho0 ]

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedHarmoniumDensity
    :: forall f z x . ( Bilinear f z x, ExponentialFamily (Harmonium f z x), Map Mean Natural f z x
       , Legendre Natural x, Legendre Natural z, ExponentialFamily x, ExponentialFamily z )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f z x
      -> SamplePoint z
      -> Double
{-# INLINE logConjugatedHarmoniumDensity #-}
logConjugatedHarmoniumDensity (rho0,rprms) hrm ox =
    let (f,nl0) = splitBottomHarmonium hrm
        (no,nlo) = splitAffine f
        nl = fromOneHarmonium nl0
     in log (baseMeasure (Proxy @ z) ox) + sum
            [ sufficientStatistic ox <.> no
            , potential (nl <+> ox *<.< nlo)
            , negate $ potential (nl <+> rprms) + rho0 ]


-- | Computes the negative log-likelihood of a sample point of a mixture model.
mixtureDensity
    :: ( Enum e, KnownNat k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE mixtureDensity #-}
mixtureDensity hrm =
    let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
     in conjugatedHarmoniumDensity rh0rx hrm

-- | Computes the negative log-likelihood of a sample point of a mixture model.
logMixtureDensity
    :: ( Enum e, KnownNat k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE logMixtureDensity #-}
logMixtureDensity hrm =
    let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
     in logConjugatedHarmoniumDensity rh0rx hrm


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
     in baseMeasure (Proxy @ z) z * exp (nz <.> mz + potential (nx <+> mz <.< nzx))

harmoniumEmpiricalExpectations
    :: ( ExponentialFamily z, Map Mean Natural f x z
       , Bilinear f z x, Map Mean Natural f z x, Legendre Natural x)
    => Sample z -- ^ Model Samples
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Mean # Harmonium f z x -- ^ Harmonium expected sufficient statistics
{-# INLINE harmoniumEmpiricalExpectations #-}
harmoniumEmpiricalExpectations zs hrm =
    let mzs = sufficientStatistic <$> zs
        aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        mxs = dualTransition <$> aff >$>* zs
        mzx = averagePoint $ zipWith (>.<) mzs mxs
        maff = joinAffine (averagePoint mzs) mzx
     in joinBottomHarmonium maff . toOneHarmonium $ averagePoint mxs


--- Instances ---


instance Manifold m => Manifold (DeepHarmonium fs '[m]) where
    type Dimension (DeepHarmonium fs '[m]) = Dimension m

instance (Manifold z, Manifold y, Manifold (f z y), Manifold (DeepHarmonium fs (y : xs)))
  => Manifold (DeepHarmonium (f : fs) (z : y : xs)) where
      type Dimension (DeepHarmonium (f : fs) (z : y : xs))
        = Dimension z + Dimension (f z y) + Dimension (DeepHarmonium fs (y : xs))

instance Manifold (DeepHarmonium fs xs) => Statistical (DeepHarmonium fs xs) where
    type SamplePoint (DeepHarmonium fs xs) = HList (SamplePoints xs)

instance Generative c x => Generative c (OneHarmonium x) where
    {-# INLINE samplePoint #-}
    samplePoint = fmap (:+: Null) . samplePoint . fromOneHarmonium

instance ExponentialFamily x => ExponentialFamily (OneHarmonium x) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (x :+: Null) =
          toOneHarmonium $ sufficientStatistic x
      {-# INLINE baseMeasure #-}
      baseMeasure = harmoniumBaseMeasure Proxy

instance ( ExponentialFamily z, ExponentialFamily y, Bilinear f z y
         , ExponentialFamily (DeepHarmonium fs (y : xs)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (z : y : xs)) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              pm = sufficientStatistic xm
              pn = sufficientStatistic xn
           in joinBottomHarmonium (joinAffine pm $ pm >.< pn) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f z x, ExponentialFamily z, Generative Natural x )
  => Gibbs '[f] '[z,x] where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zxs = initialPass dhrm . fst $ hUnzip zxs
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              xp = fromOneHarmonium dhrm'
          xs <- mapM (samplePoint . (<+>) xp) $ zs *<$< f
          return . hZip zs . hZip xs $ repeat Null

instance ( Bilinear f z y, Map Mean Natural g y x, Manifold (DeepHarmonium fs (x : xs))
         , ExponentialFamily z, ExponentialFamily x, Generative Natural y, Gibbs (g : fs) (y : x : xs) )
  => Gibbs (f : g : fs) (z : y : x : xs) where
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

instance Manifold x => TransposeHarmonium '[] '[x] where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium = id

instance (Bilinear f z y, Bilinear f z y, TransposeHarmonium fs (y : xs))
  => TransposeHarmonium (f : fs) (z : y : xs) where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium dhrm =
        let (aff,dhrm') = splitBottomHarmonium dhrm
            (pm,pmtx) = splitAffine aff
            dhrm'' = transposeHarmonium dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

instance Generative Natural x => SampleConjugated '[] '[x] where
    {-# INLINE sampleConjugatedHarmonium #-}
    sampleConjugatedHarmonium k _ = sample k

instance ( Manifold (DeepHarmonium fs (y : xs)), Map Mean Natural f z y, Manifold (Sum xs)
         , ExponentialFamily y, SampleConjugated fs (y : xs), Generative Natural z )
  => SampleConjugated (f : fs) (z : y : xs) where
    {-# INLINE sampleConjugatedHarmonium #-}
    sampleConjugatedHarmonium k rprms dhrm = do
        let (pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleConjugatedHarmonium k rprms' $ biasBottom rprm dhrm'
        zs <- mapM samplePoint $ pf >$>* ys
        return . hZip zs $ hZip ys xs

instance ( Enum e, KnownNat n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
  => Generative Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE samplePoint #-}
      samplePoint hrm = head <$> sampleMixtureModel 1 hrm

instance ( Enum e, KnownNat n, Legendre Natural o, ExponentialFamily o
  , Manifold (Harmonium Tensor o (Categorical e n)))
  => Legendre Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE potential #-}
      potential hrm =
          let (lkl,nx0) = splitBottomHarmonium hrm
              (rho0,rprms) = mixtureLikelihoodConjugationParameters lkl
              nx = fromOneHarmonium nx0
           in potential (nx <+> rprms) + rho0
      potentialDifferential = primalIsomorphism . mixtureExpectations

--- Graveyard

---- | The observable density of a categorical harmonium.
--mixtureDensity0
--    :: (KnownNat k, Enum e, Legendre Natural z, AbsolutelyContinuous Natural z)
--    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
--    -> SamplePoint z -- ^ Observation
--    -> Double -- ^ Probablity density of the observation
--{-# INLINE mixtureDensity0 #-}
--mixtureDensity0 hrm x =
--    let (affzx,nx) = splitBottomHarmonium hrm
--        nz = fst $ splitAffine affzx
--        wghts = listCoordinates . toMean
--            $ snd (mixtureLikelihoodConjugationParameters affzx) <+> fromOneHarmonium nx
--        dxs0 = (`density` x) <$> affzx >$>* pointSampleSpace (fromOneHarmonium nx)
--        dx1 = density nz x * (1 - sum wghts)
--     in dx1 + sum (zipWith (*) wghts dxs0)
--
--

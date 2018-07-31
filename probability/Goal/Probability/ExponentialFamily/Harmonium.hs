{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums and Rectification.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , DeepHarmonium
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    -- ** Construction
    , biasBottom
    , getBottomBias
    , splitBottomHarmonium
    , joinBottomHarmonium
    -- * Sampling
    , Gibbs (upwardPass,initialPass)
    , gibbsPass
    -- ** Transposition
    , TransposeHarmonium (transposeHarmonium)
     -- * Rectification
    , marginalizeRectifiedHarmonium
    , SampleRectified (sampleRectified)
    -- * Mixture Models
    , buildMixtureModel
    , splitMixtureModel
    , mixtureDensity
    , categoricalLikelihoodRectificationParameters
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic.Internal as I


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
    => Dual c ~> c # Affine f m n -- ^ Affine function
    -> c # DeepHarmonium fs (n : ms) -- ^ Upper part of the deep harmonium
    -> c # DeepHarmonium (f : fs) (m : n : ms) -- ^ Combined deep harmonium
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold m, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep Harmonium
    -> (Dual c ~> c # Affine f m n, c # DeepHarmonium fs (n : ms)) -- ^ Affine function and upper part
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms))
      , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
    => c # m -- ^ Bias step
    -> c # DeepHarmonium fs (m : ms) -- ^ Deep Harmonium
    -> c # DeepHarmonium fs (m : ms) -- ^ Biased deep harmonium
{-# INLINE biasBottom #-}
biasBottom pm' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double
        (pmcs,css') = S.splitAt $ coordinates dhrm
        pm = pm' <+> Point pmcs
     in Point $ coordinates pm S.++ css'

-- | Get the bias of the bottom layer of the given 'DeepHarmonium'.
getBottomBias
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms))
      , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
    => c # DeepHarmonium fs (m : ms) -- ^ Deep harmonium
    -> c # m -- ^ Bottom layer bias
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let (pmcs,_ :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double)
          = S.splitAt $ coordinates dhrm
       in Point pmcs


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where

    -- | Given a 'DeepHarmonium' and an element of its sample space, partially
    -- updates the sample by resampling from the bottom to the top layer, but
    -- without updating the bottom layer itself.
    -- harmonium.
    upwardPass
        :: KnownNat l
        => Natural # DeepHarmonium fs ms -- ^ Deep harmonium
        -> Sample l (DeepHarmonium fs ms) -- ^ Initial sample
        -> Random s (Sample l (DeepHarmonium fs ms)) -- ^ Partial Gibbs resample

    -- | Generates an element of the sample spaec of a deep harmonium based by
    -- starting from a sample point from the bottom layer, and doing a naive
    -- upward sampling. This does not generate a true sample from the deep
    initialPass
        :: KnownNat l
        => Natural # DeepHarmonium fs ms -- Deep harmonium
        -> Sample l (Head ms) -- ^ Bottom layer sample
        -> Random s (Sample l (DeepHarmonium fs ms)) -- ^ Initial deep harmonium sample

-- | Harmonium transpotion. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
class Manifold (DeepHarmonium fs ms) => TransposeHarmonium fs ms where
    transposeHarmonium :: Primal c => c # DeepHarmonium fs ms -> c # DeepHarmonium (Reverse3 fs) (Reverse ms)


-- | A single pass of Gibbs sampling. Infinite, recursive application of this
-- function yields a sample from the given 'DeepHarmonium'.
gibbsPass :: ( KnownNat k, Manifold (DeepHarmonium fs (n : ms))
             , Gibbs (f : fs) (m : n : ms), Map Mean Natural f m n
             , Generative Natural m, ExponentialFamily n )
  => Natural # DeepHarmonium (f : fs) (m : n : ms) -- ^ Deep Hamonium
  -> B.Vector k (HList (z : SamplePoint n : SamplePoints ms)) -- ^ Initial Sample
  -> Random s (B.Vector k (HList (SamplePoint m : SamplePoint n : SamplePoints ms))) -- ^ Gibbs resample
{-# INLINE gibbsPass #-}
gibbsPass dhrm zyxs = do
    let yxs = snd $ hUnzip zyxs
        ys = fst $ hUnzip yxs
        f = fst $ splitBottomHarmonium dhrm
    zs <- samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip zs yxs


--- Rectification ---


-- | A rectified distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleRectified fs ms where
    -- | A true sample from a rectified harmonium.
    sampleRectified
        :: KnownNat l
        => Natural # Sum (Tail ms) -- ^ Rectification parameters
        -> Natural # DeepHarmonium fs ms -- ^ Deep harmonium
        -> Random s (Sample l (DeepHarmonium fs ms)) -- ^ Deep harmonium sample

-- | Marginalize the bottom layer out of a deep harmonium.
marginalizeRectifiedHarmonium
    :: ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
       , ExponentialFamily m, Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
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
mixtureDensity
    :: (KnownNat k, 1 <= k, Num e, Enum e, Legendre Natural z, Transition Source Natural z, AbsolutelyContinuous Natural z)
    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
    -> SamplePoint z -- ^ Observation
    -> Double -- ^ Probablity density of the observation
{-# INLINE mixtureDensity #-}
mixtureDensity hrm x =
    let (affzx,nx) = splitBottomHarmonium hrm
        nz = fst $ splitAffine affzx
        wghts = coordinates . toMean $ snd (categoricalLikelihoodRectificationParameters affzx) <+> fromOneHarmonium nx
        dxs0 = mapReplicated (`density` x) $ affzx >$>* B.enumFromN 0
        dx1 = density nz x * (1 - S.sum wghts)
     in dx1 + S.sum (S.zipWith (*) wghts dxs0)

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
        rprms = snd $ categoricalLikelihoodRectificationParameters affzx
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
        rprms = snd $ categoricalLikelihoodRectificationParameters affzx
        nx0 = fromOneHarmonium nx <+> rprms
        (nz,nzx) = splitAffine affzx
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (<+> nz) nzs
        nz0 = S.singleton nz
     in (nzs0' S.++ nz0,nx0)

-- | Computes the rectification parameters of a likelihood defined by a categorical latent variable.
categoricalLikelihoodRectificationParameters
    :: (KnownNat k, 1 <= k, Enum e, Legendre Natural z)
    => Mean ~> Natural # z <* Categorical e k -- ^ Categorical likelihood
    -> (Double, Natural # Categorical e k) -- ^ Rectification parameters
{-# INLINE categoricalLikelihoodRectificationParameters #-}
categoricalLikelihoodRectificationParameters aff =
    let (nz,nzx) = splitAffine aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz <+> Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleMixtureModel
    :: ( KnownNat k, Enum e, KnownNat n, 1 <= n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
      => Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
      -> Random s (Sample k (Harmonium Tensor o (Categorical e n))) -- ^ Sample
{-# INLINE sampleMixtureModel #-}
sampleMixtureModel hrm = do
    let rx = snd . categoricalLikelihoodRectificationParameters . fst $ splitBottomHarmonium hrm
    sampleRectified (toSingletonSum rx) hrm


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

mixtureModelExpectations
    :: ( Enum e, KnownNat k, 1 <= k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k)
    -> Mean # Harmonium Tensor o (Categorical e k)
{-# INLINE mixtureModelExpectations #-}
mixtureModelExpectations hrm =
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
          xs <- samplePoint . mapReplicatedPoint (<+> xp) $ zs *<$< f
          return . hZip zs . hZip xs $ B.replicate Null

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
          ys <- samplePoint $ g >$>* xs <+> zs *<$< f
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              yp = fst . splitAffine . fst $ splitBottomHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ zs *<$< f
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
    {-# INLINE sampleRectified #-}
    sampleRectified _ = sample

instance ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
         , ExponentialFamily n, SampleRectified fs (n : ms), Generative Natural m
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
  => SampleRectified (f : fs) (m : n : ms) where
    {-# INLINE sampleRectified #-}
    sampleRectified rprms dhrm = do
        let (pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleRectified rprms' $ biasBottom rprm dhrm'
        zs <- samplePoint $ pf >$>* ys
        return . hZip zs $ hZip ys xs

instance ( Enum e, KnownNat n, 1 <= n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
  => Generative Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE samplePoint #-}
      samplePoint hrm = do
          (smp :: Sample 1 (Harmonium Tensor o (Categorical e n))) <- sampleMixtureModel hrm
          return $ B.head smp

instance ( Enum e, KnownNat n, 1 <= n, Legendre Natural o, ExponentialFamily o
  , Manifold (Harmonium Tensor o (Categorical e n)))
  => Legendre Natural (Harmonium Tensor o (Categorical e n)) where
      {-# INLINE potential #-}
      potential hrm =
          let (lkl,nx0) = splitBottomHarmonium hrm
              nx = fromOneHarmonium nx0
           in log $ sum [ exp (sufficientStatistic i <.> nx) + potential (lkl >.>* i)
                        | i <- B.toList $ pointSampleSpace nx ]
      potentialDifferential = breakPoint . mixtureModelExpectations

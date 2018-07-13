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
    , Gibbs (upwardPass, initialPass)
    , gibbsPass
    , buildCategoricalHarmonium
    , mixtureDensity
    -- ** Transposition
    , TransposeHarmonium (transposeHarmonium)
     -- * Rectification
    , marginalizeRectifiedHarmonium
    , SampleRectified (sampleRectified)
      -- ** Categorical Harmoniums
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

-- | Adds a layer to the top of a deep harmonium.
joinBottomHarmonium
    :: (Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => Dual c ~> c # Affine f m n
    -> c # DeepHarmonium fs (n : ms)
    -> c # DeepHarmonium (f : fs) (m : n : ms)
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold m, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms)
    -> (Dual c ~> c # Affine f m n, c # DeepHarmonium fs (n : ms))
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Translate the bias of the top layer by the given 'Point'.
biasBottom
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms))
      , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
    => c # m
    -> c # DeepHarmonium fs (m : ms)
    -> c # DeepHarmonium fs (m : ms)
{-# INLINE biasBottom #-}
biasBottom pm' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double
        (pmcs,css') = S.splitAt $ coordinates dhrm
        pm = pm' <+> Point pmcs
     in Point $ coordinates pm S.++ css'

-- | Translate the bias of the top layer by the given 'Point'.
getBottomBias
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms))
      , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
    => c # DeepHarmonium fs (m : ms)
    -> c # m
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let (pmcs,_ :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double)
          = S.splitAt $ coordinates dhrm
       in Point pmcs


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
    upwardPass :: KnownNat l
           => Natural # DeepHarmonium fs ms
           -> Sample l (DeepHarmonium fs ms)
           -> Random s (Sample l (DeepHarmonium fs ms))
    initialPass :: KnownNat l
           => Natural # DeepHarmonium fs ms
           -> Sample l (Head ms)
           -> Random s (Sample l (DeepHarmonium fs ms))

-- | Harmonium transpotion. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
class Manifold (DeepHarmonium fs ms) => TransposeHarmonium fs ms where
    transposeHarmonium :: Primal c => c # DeepHarmonium fs ms -> c # DeepHarmonium (Reverse3 fs) (Reverse ms)


-- | A single pass of Gibbs sampling. Infinite, recursive application of this function yields a sample from the given 'DeepHarmonium'.
gibbsPass :: ( KnownNat k, Manifold (DeepHarmonium fs (n : ms))
             , Gibbs (f : fs) (m : n : ms), Map Mean Natural f m n
             , Generative Natural m, ExponentialFamily n )
  => Natural # DeepHarmonium (f : fs) (m : n : ms)
  -> B.Vector k (HList (x : SamplePoint n : SamplePoints ms))
  -> Random s (B.Vector k (HList (SamplePoint m : SamplePoint n : SamplePoints ms)))
{-# INLINE gibbsPass #-}
gibbsPass dhrm xyzs = do
    let yzs = snd $ hUnzip xyzs
        ys = fst $ hUnzip yzs
        f = fst $ splitBottomHarmonium dhrm
    xs <- samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip xs yzs


--- Rectification ---


-- | A rectified distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleRectified fs ms where
    sampleRectified
        :: KnownNat l
        => Natural # Sum (Tail ms)
        -> Natural # DeepHarmonium fs ms
        -> Random s (Sample l (DeepHarmonium fs ms))

marginalizeRectifiedHarmonium
    :: ( Manifold (DeepHarmonium fs (n : ms)), Map Mean Natural f m n, Manifold (Sum ms)
       , ExponentialFamily m, Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # Sum (n : ms)
      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> (Natural # Sum ms, Natural # DeepHarmonium fs (n : ms))
{-# INLINE marginalizeRectifiedHarmonium #-}
marginalizeRectifiedHarmonium rprms dhrm =
    let dhrm' = snd $ splitBottomHarmonium dhrm
        (rprm,rprms') = splitSum rprms
     in (rprms', biasBottom rprm dhrm')



-- Categorical Harmoniums --

mixtureDensity
    :: (KnownNat k, 1 <= k, Num e, Enum e, Legendre Natural z, Transition Source Natural z, AbsolutelyContinuous Natural z)
    => Natural # Harmonium Tensor z (Categorical e k)
    -> SamplePoint z
    -> Double
{-# INLINE mixtureDensity #-}
mixtureDensity hrm x =
    let (affzx,nx) = splitBottomHarmonium hrm
        nz = fst $ splitAffine affzx
        wghts = coordinates . toMean $ snd (categoricalLikelihoodRectificationParameters affzx) <+> fromOneHarmonium nx
        dxs0 = mapReplicated (`density` x) $ affzx >$>* B.enumFromN 0
        dx1 = density nz x * (1 - S.sum wghts)
     in dx1 + S.sum (S.zipWith (*) wghts dxs0)

-- | A convenience function for building a mixture model.
buildCategoricalHarmonium
    :: forall c d k e z
    . ( KnownNat k, 1 <= k, Enum e, Legendre Natural z
      , Transition c Natural z, Transition d Natural (Categorical e k) )
    => c # z -- -- ^ Component Bias
    -> S.Vector k (c # z) -- ^ Mixture components
    -> d # Categorical e k
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
        affzx = joinAffine nz nzx
        nx' = snd $ categoricalLikelihoodRectificationParameters affzx
        nx = toOneHarmonium $ toNatural mx <-> nx'
     in joinBottomHarmonium affzx nx

-- | Computes the rectification parameters of a harmonium with a categorical latent variable.
categoricalLikelihoodRectificationParameters
    :: (KnownNat k, 1 <= k, Enum e, Legendre Natural z)
    => Mean ~> Natural # z <* Categorical e k
    -> (Double, Natural # Categorical e k)
{-# INLINE categoricalLikelihoodRectificationParameters #-}
categoricalLikelihoodRectificationParameters aff =
    let (nz,nzx) = splitAffine aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz <+> Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleCategoricalHarmonium
    :: ( KnownNat k, Enum e, KnownNat n, 1 <= n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
      => Point Natural (Harmonium Tensor o (Categorical e n))
      -> Random s (Sample k (Harmonium Tensor o (Categorical e n)))
{-# INLINE sampleCategoricalHarmonium #-}
sampleCategoricalHarmonium hrm = do
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

instance (Bilinear Tensor m n, Bilinear Tensor n m, TransposeHarmonium fs (n : ms))
  => TransposeHarmonium (Tensor : fs) (m : n : ms) where
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
          (smp :: Sample 1 (Harmonium Tensor o (Categorical e n))) <- sampleCategoricalHarmonium hrm
          return $ B.head smp

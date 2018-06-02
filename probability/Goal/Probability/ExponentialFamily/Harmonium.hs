{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , DeepHarmonium
    -- ** Inference
    , (<|<)
    , (<|<*)
    -- ** Sampling
    , Gibbs (upwardPass, initialPass)
    , gibbsPass
    -- ** Transposition
    , TransposeHarmonium (transposeHarmonium)
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    -- ** Construction
    , biasBottom
    , splitBottomHarmonium
    , joinBottomHarmonium
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily

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
fromOneHarmonium :: c # DeepHarmonium '[] '[m] -> c # m
{-# INLINE fromOneHarmonium #-}
fromOneHarmonium = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toOneHarmonium :: c # m -> c # DeepHarmonium '[] '[m]
{-# INLINE toOneHarmonium #-}
toOneHarmonium = breakPoint

-- | Adds a layer to the top of a deep harmonium.
joinBottomHarmonium
    :: (Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # m
    -> Dual c ~> c # f m n
    -> c # DeepHarmonium fs (n : ms)
    -> c # DeepHarmonium (f : fs) (m : n : ms)
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pm pf dhrm =
    Point $ coordinates pm S.++ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold m, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms)
    -> (c # m, Dual c ~> c # f m n, c # DeepHarmonium fs (n : ms))
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (lcs,css') = S.splitAt $ coordinates dhrm
        (mtxcs,dcs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point dcs)

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

-- | The given deep harmonium conditioned on its bottom layer.
(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms))
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Mean # m
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (_,f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (p <.< f) dhrm'

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
        yps = sufficientStatistic ys
        (xp,f,_) = splitBottomHarmonium dhrm
    xs <- samplePoint . mapReplicatedPoint (<+> xp) $ f >$> yps
    upwardPass dhrm $ hZip xs yzs



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
           in joinBottomHarmonium pm (pm >.< pn) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n )
  => Gibbs '[f] '[m,n] where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zxs = do
          let zs = fst $ hUnzip zxs
              zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              xp = fromOneHarmonium dhrm'
          xs <- samplePoint . mapReplicatedPoint (<+> xp) $ zps <$< f
          return . hZip zs . hZip xs $ B.replicate Null
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              xp = fromOneHarmonium dhrm'
          xs <- samplePoint . mapReplicatedPoint (<+> xp) $ zps <$< f
          return . hZip zs . hZip xs $ B.replicate Null

instance ( Bilinear f m n, Map Mean Natural g n o, Manifold (DeepHarmonium fs (o : ms))
         , ExponentialFamily m, ExponentialFamily o, Generative Natural n, Gibbs (g : fs) (n : o : ms) )
  => Gibbs (f : g : fs) (m : n : o : ms) where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zyxs = do
          let (zs,yxs) = hUnzip zyxs
              (xs,xs') = hUnzip . snd $ hUnzip yxs
              xps = sufficientStatistic xs
              zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              (yp,g,_) = splitBottomHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ g >$> xps <+> zps <$< f
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              (yp,_,_) = splitBottomHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ zps <$< f
          yxs' <- initialPass dhrm' ys
          return $ hZip zs yxs'

instance Manifold m => TransposeHarmonium '[] '[m] where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium = id

instance (Bilinear Tensor m n, Bilinear Tensor n m, TransposeHarmonium fs (n : ms))
  => TransposeHarmonium (Tensor : fs) (m : n : ms) where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium dhrm =
        let (pm,pmtx,dhrm') = splitBottomHarmonium dhrm
            dhrm'' = transposeHarmonium dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

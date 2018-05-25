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


--- Types ---


-- | A hierarchical generative model defined by exponential families. Note that
-- the first elements of ms is the bottom layer of the hierachy, and each
-- subsequent element is the next layer up.
data DeepHarmonium (fs :: [* -> * -> *]) (ms :: [*])

-- | A trivial 1-layer harmonium.
type OneHarmonium m = DeepHarmonium '[] '[m]

-- | A 2-layer harmonium.
type Harmonium f n m = DeepHarmonium '[f] [n,m]


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
    :: (Manifold (f m n), Manifold (DeepHarmonium fs (m : ms)))
    => c # n
    -> Dual c ~> c # f m n
    -> c # DeepHarmonium fs (m : ms)
    -> c # DeepHarmonium (f : fs) (n : m : ms)
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pm pf dhrm =
    Point $ coordinates pm S.++ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold n, Manifold (f m n), Manifold (DeepHarmonium fs (m : ms)))
    => c # DeepHarmonium (f : fs) (n : m : ms)
    -> (c # n, Dual c ~> c # f m n, c # DeepHarmonium fs (m : ms))
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
(<|<) :: ( Map Mean Natural f m n, Manifold (DeepHarmonium fs (m : ms))
         , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
      => Natural # DeepHarmonium (f : fs) (n : m : ms)
      -> Mean # n
      -> Natural # DeepHarmonium fs (m : ms)
{-# INLINE (<|<) #-}
(<|<) dhrm p =
    let (_,f,dhrm') = splitBottomHarmonium dhrm
     in biasBottom (f >.> p) dhrm'

-- | The given deep harmonium conditioned on a sample from its bottom layer.
-- This can be interpreted as the posterior of the model given an observation of
-- the bottom layer.
(<|<*) :: ( Map Mean Natural f m n, Manifold (DeepHarmonium fs (m : ms)), ExponentialFamily n
         , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
      => Natural # DeepHarmonium (f : fs) (n : m : ms)
      -> SamplePoint n
      -> Natural # DeepHarmonium fs (m : ms)
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


-- | A single pass of Gibbs sampling. Infinite, recursive application of this function yields a sample from the given 'DeepHarmonium'.
gibbsPass :: ( KnownNat k, Manifold (DeepHarmonium fs (m : ms))
             , Gibbs (f : fs) (n : m : ms), Map Mean Natural f m n
             , Generative Natural n, ExponentialFamily m, Bilinear f m n )
  => Natural # DeepHarmonium (f : fs) (n : m : ms)
  -> B.Vector k (HList (x : SamplePoint m : SamplePoints ms))
  -> Random s (B.Vector k (HList (SamplePoint n : SamplePoint m : SamplePoints ms)))
{-# INLINE gibbsPass #-}
gibbsPass dhrm xyzs = do
    let yzs = snd $ hUnzip xyzs
        ys = fst $ hUnzip yzs
        yps = sufficientStatistic ys
        (xp,f,_) = splitBottomHarmonium dhrm
    xs <- samplePoint . mapReplicatedPoint (<+> xp) $ yps <$< f
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

instance (Bilinear f m n, Manifold (DeepHarmonium fs (m : ms)))
  => Manifold (DeepHarmonium (f : fs) (n : m : ms)) where
      type Dimension (DeepHarmonium (f : fs) (n : m : ms))
        = Dimension n + Dimension (f m n) + Dimension (DeepHarmonium fs (m : ms))

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
         , Bilinear f m n, ExponentialFamily (DeepHarmonium fs (m : ms)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (n : m : ms)) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xn :+: xm :+: xs) =
          let mdhrm = sufficientStatistic $ xm :+: xs
              pn = sufficientStatistic xn
              pm = sufficientStatistic xm
           in joinBottomHarmonium pn (pm >.< pn) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Map Mean Natural f m n, ExponentialFamily n, Generative Natural m )
  => Gibbs '[f] '[n,m] where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zxs = do
          let zs = fst $ hUnzip zxs
              zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              xp = fromOneHarmonium dhrm'
          xs <- samplePoint . mapReplicatedPoint (<+> xp) $ f >$> zps
          return . hZip zs . hZip xs $ B.replicate Null
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              xp = fromOneHarmonium dhrm'
          xs <- samplePoint . mapReplicatedPoint (<+> xp) $ f >$> zps
          return . hZip zs . hZip xs $ B.replicate Null

instance ( Map Mean Natural f n o, Bilinear g m n, Manifold (DeepHarmonium  fs ( m : ms))
         , ExponentialFamily m, ExponentialFamily o, Generative Natural n, Gibbs (g : fs) (n : m : ms) )
  => Gibbs (f : g : fs) (o : n : m : ms) where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zyxs = do
          let (zs,yxs) = hUnzip zyxs
              (xs,xs') = hUnzip . snd $ hUnzip yxs
              xps = sufficientStatistic xs
              zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              (yp,g,_) = splitBottomHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ f >$> zps <+> xps <$< g
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let zps = sufficientStatistic zs
              (_,f,dhrm') = splitBottomHarmonium dhrm
              (yp,_,_) = splitBottomHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ f >$> zps
          yxs' <- initialPass dhrm' ys
          return $ hZip zs yxs'

{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , type (<*>)
    , DeepHarmonium
    -- ** Manipulation
    , (>|>)
    -- , (<|<)
    -- ** Sampling
    , gibbsPass
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
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

-- | An infix synonym for a 'Harmonium'.
type (m <*> n) = DeepHarmonium '[Tensor] [m,n]


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
joinHeadHarmonium
    :: (Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # m
    -> Dual c ~> c # f m n
    -> c # DeepHarmonium fs (n : ms)
    -> c # DeepHarmonium (f : fs) (m : n : ms)
{-# INLINE joinHeadHarmonium #-}
joinHeadHarmonium pl pmtx dhrm =
    Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitHeadHarmonium
    :: (Manifold m, Manifold (f m n), Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms)
    -> (c # m, Dual c ~> c # f m n, c # DeepHarmonium fs (n : ms))
{-# INLINE splitHeadHarmonium #-}
splitHeadHarmonium dhrm =
    let (lcs,css') = S.splitAt $ coordinates dhrm
        (mtxcs,dcs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point dcs)

-- | Adds a layer to the top of a deep harmonium.
-- | Splits the bottom layer off of a transposed harmonium.
splitLastHarmonium
    :: (Map Mean Natural (Last3 fs) (Last (Init ms)) (Last ms), Manifold (DeepHarmonium (Init3 fs) (Init ms)))
    => c # DeepHarmonium fs ms
    -> ( c # DeepHarmonium (Init3 fs) (Init ms)
       , Dual c ~> c # (Last3 fs) (Last (Init ms)) (Last ms)
       , c # Last ms )
{-# INLINE splitLastHarmonium #-}
splitLastHarmonium dhrm =
    let (dcs,css') = S.splitAt . I.Vector . S.fromSized $ coordinates dhrm
        (mtxcs,lcs) = S.splitAt css'
     in (Point dcs, Point mtxcs, Point lcs)


-- | Translate the bias of the top layer by the given 'Point'.
biasTop
    :: forall fs m ms c
    . ( Manifold m, Manifold (DeepHarmonium fs (m : ms))
      , Dimension m <= Dimension (DeepHarmonium fs (m : ms)) )
    => c # m
    -> c # DeepHarmonium fs (m : ms)
    -> c # DeepHarmonium fs (m : ms)
{-# INLINE biasTop #-}
biasTop pm' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double
        (pmcs,css') = S.splitAt $ coordinates dhrm
        pm = pm' <+> Point pmcs
     in Point $ coordinates pm S.++ css'

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs ms c . ( Manifold (DeepHarmonium fs ms), Manifold (Last ms)
                        , Dimension (Last ms) <= Dimension (DeepHarmonium fs ms) )
    => c # Last ms
    -> c # DeepHarmonium fs ms
    -> c # DeepHarmonium fs ms
{-# INLINE biasBottom #-}
biasBottom pn' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs ms) - Dimension (Last ms)) Double
        (css',pncs) = S.splitAt $ coordinates dhrm
        pn = pn' <+> Point pncs
     in Point $ css' S.++ coordinates pn

-- | The given deep harmonium conditioned on its top layer.
(>|>) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms))
         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
      => Mean # m
      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (>|>) #-}
(>|>) p dhrm =
    let (_,f,dhrm') = splitHeadHarmonium dhrm
     in biasTop (p <.< f) dhrm'

---- | The given deep harmonium conditioned on its bottom layer. This usually
---- corresponds to the posterior of the deep harmonium.
--(<|<) :: ( Bilinear f m n, Manifold (DeepHarmonium fs (n : ms))
--         , Dimension n <= Dimension (DeepHarmonium fs (n : ms)) )
--      => Mean # m
--      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
--      -> Natural # DeepHarmonium fs (n : ms)
--{-# INLINE (>|>) #-}
--(<|<) p dhrm =
--    let (_,f,dhrm') = splitHeadHarmonium dhrm
--     in biasTop (p <.< f) dhrm'


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
    downwardPass :: KnownNat l
           => Natural # DeepHarmonium fs ms
           -> Sample l (DeepHarmonium fs ms)
           -> Random s (Sample l (DeepHarmonium fs ms))

gibbsPass :: ( KnownNat k, Manifold (DeepHarmonium fs (n : ms))
             , Gibbs (f : fs) (m : n : ms), Map Mean Natural f m n
             , Generative Natural m, ExponentialFamily n, Bilinear f m n )
  => Natural # DeepHarmonium (f : fs) (m : n : ms)
  -> B.Vector k (HList (x : SamplePoint n : SamplePoints ms))
  -> Random s (B.Vector k (HList (SamplePoint m : SamplePoint n : SamplePoints ms)))
gibbsPass dhrm xyzs = do
    let yzs = snd $ hUnzip xyzs
        ys = fst $ hUnzip yzs
        yps = sufficientStatistic ys
        (xp,f,_) = splitHeadHarmonium dhrm
    xs <- samplePoint . mapReplicatedPoint (<+> xp) $ f >$> yps
    downwardPass dhrm $ hZip xs yzs



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

instance (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
  => Manifold (DeepHarmonium (f : fs) (m : n : ms)) where
      type Dimension (DeepHarmonium (f : fs) (m : n : ms))
        = Dimension m + Dimension (f m n) + Dimension (DeepHarmonium fs (n : ms))

type family SamplePoints (ms :: [*]) where
    SamplePoints '[] = '[]
    SamplePoints (m : ms) = SamplePoint m : SamplePoints ms

instance Manifold (DeepHarmonium fs ms) => Statistical (DeepHarmonium fs ms) where
    type SamplePoint (DeepHarmonium fs ms) = HList (SamplePoints ms)

instance Generative c m => Generative c (OneHarmonium m) where
    {-# INLINE samplePoint #-}
    samplePoint = fmap (:+: Null) . samplePoint . fromOneHarmonium

instance ExponentialFamily m => ExponentialFamily (OneHarmonium m) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xl :+: Null) =
          toOneHarmonium $ sufficientStatistic xl
      {-# INLINE baseMeasure #-}
      baseMeasure = harmoniumBaseMeasure Proxy

instance ( ExponentialFamily m, ExponentialFamily n
         , Bilinear f m n, ExponentialFamily (DeepHarmonium fs (n : ms)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (m : n : ms)) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              mm = sufficientStatistic xm
              mh = sufficientStatistic xn
           in joinHeadHarmonium mm (mm >.< mh) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n )
  => Gibbs '[f] '[m,n] where
      {-# INLINE downwardPass #-}
      downwardPass dhrm xzs = do
          let xs = fst $ hUnzip xzs
              xps = sufficientStatistic xs
              (_,f,dhrm') = splitHeadHarmonium dhrm
              zp = fromOneHarmonium dhrm'
          zs <- samplePoint . mapReplicatedPoint (<+> zp) $ xps <$< f
          return . hZip xs . hZip zs $ B.replicate Null

instance ( Bilinear f m n, Map Mean Natural g n o, Bilinear g n o, Manifold (DeepHarmonium  fs ( o : ms))
         , ExponentialFamily m, ExponentialFamily o, Generative Natural n, Gibbs (g : fs) (n : o : ms) )
  => Gibbs (f : g : fs) (m : n : o : ms) where
      {-# INLINE downwardPass #-}
      downwardPass dhrm xyzs = do
          let (xs,yzs) = hUnzip xyzs
              (zs,zs') = hUnzip . snd $ hUnzip yzs
              xps = sufficientStatistic xs
              zps = sufficientStatistic zs
              (_,f,dhrm') = splitHeadHarmonium dhrm
              (yp,g,_) = splitHeadHarmonium dhrm'
          ys <- samplePoint . mapReplicatedPoint (<+> yp) $ g >$> zps <+> xps <$< f
          yzs' <- downwardPass dhrm' . hZip ys $ hZip zs zs'
          return $ hZip xs yzs'

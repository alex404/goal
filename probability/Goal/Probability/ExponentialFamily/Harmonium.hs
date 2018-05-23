{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , type (<*>)
    , DeepHarmonium
    , Hierarchical
    -- ** Manipulation
    , (>|>)
    , (<|<)
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


--- Types ---


-- | A hierarchical generative model defined by exponential families.
data DeepHarmonium (fs :: [* -> * -> *]) (ms :: [*])

-- | A hierarchical generative model defined by exponential families.
data HarmoniumTranspose (fs :: [* -> * -> *]) (ns :: [*])

-- | A trivial 1-layer harmonium.
type OneHarmonium m = DeepHarmonium '[] '[m]

-- | A 2-layer harmonium.
type Harmonium f m n = DeepHarmonium '[f] [m,n]

-- | An infix synonym for a 'Harmonium'.
type (m <*> n) = DeepHarmonium '[Tensor] [m,n]

-- | 'Hierarchical' is a collection of constraints for making arbitrarily-deep
-- harmoniums behave.
type Hierarchical fs ms =
    ( Manifold (DeepHarmonium fs ms)
    , Dimension (Head ms) <= Dimension (DeepHarmonium fs ms)
    , Dimension (Last ms) <= Dimension (DeepHarmonium fs ms)
    , Dimension (DeepHarmonium fs ms) ~ Dimension (HarmoniumTranspose (Reverse3 fs) (Reverse ms)) )

type HierarchicalTranspose fs ns =
    ( Manifold (HarmoniumTranspose fs ns)
    , Dimension (Head ns) <= Dimension (HarmoniumTranspose fs ns) )


--- Functions ---


-- | Harmonium transposition. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
transposeHarmonium
    :: Dimension (DeepHarmonium fs ms) ~ Dimension (HarmoniumTranspose (Reverse3 fs) (Reverse ms))
    => c # DeepHarmonium fs ms
    -> c # HarmoniumTranspose (Reverse3 fs) (Reverse ms)
{-# INLINE transposeHarmonium #-}
transposeHarmonium = breakPoint

-- | Harmonium transposition. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
detransposeHarmonium
    :: Dimension (DeepHarmonium fs ms) ~ Dimension (HarmoniumTranspose (Reverse3 fs) (Reverse ms))
    => c # HarmoniumTranspose (Reverse3 fs) (Reverse ms)
    -> c # DeepHarmonium fs ms
{-# INLINE detransposeHarmonium #-}
detransposeHarmonium = breakPoint

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
    :: (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
    => c # m
    -> Dual c ~> c # f m n
    -> c # DeepHarmonium fs (n : ms)
    -> c # DeepHarmonium (f : fs) (m : n : ms)
{-# INLINE joinHeadHarmonium #-}
joinHeadHarmonium pl pmtx dhrm =
    Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitHeadHarmonium
    :: (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms)
    -> (c # m, Dual c ~> c # f m n, c # DeepHarmonium fs (n : ms))
{-# INLINE splitHeadHarmonium #-}
splitHeadHarmonium dhrm =
    let (lcs,css') = S.splitAt $ coordinates dhrm
        (mtxcs,dcs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point dcs)

-- | Adds a layer to the top of a deep harmonium.
joinLastHarmonium
    :: (Bilinear f m n, Manifold (HarmoniumTranspose fs (m : ns)))
    => c # HarmoniumTranspose fs (m : ns)
    -> Dual c ~> c # f m n
    -> c # n
    -> c # HarmoniumTranspose (f : fs) (n : m : ns)
{-# INLINE joinLastHarmonium #-}
joinLastHarmonium dhrm pmtx po =
    Point $ coordinates dhrm S.++ coordinates pmtx S.++ coordinates po

-- | Splits the bottom layer off of a transposed harmonium.
splitLastHarmonium
    :: (Bilinear f m n, Manifold (HarmoniumTranspose fs (m : ns)))
    => c # HarmoniumTranspose (f : fs) (n : m : ns)
    -> (c # HarmoniumTranspose fs (m : ns), Dual c ~> c # f m n, c # n)
{-# INLINE splitLastHarmonium #-}
splitLastHarmonium dhrm =
    let (dcs,css') = S.splitAt $ coordinates dhrm
        (mtxcs,lcs) = S.splitAt css'
     in (Point dcs, Point mtxcs, Point lcs)

-- | Translate the bias of the top layer by the given 'Point'.
biasTop
    :: forall fs m ms c . ( Hierarchical fs (m : ms), Manifold m )
    => c # m
    -> c # DeepHarmonium fs (m : ms)
    -> c # DeepHarmonium fs (m : ms)
{-# INLINE biasTop #-}
biasTop pm' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double
        (pmcs,css') = S.splitAt $ coordinates dhrm
        pm = pm' <+> Point pmcs
     in Point $ coordinates pm S.++ css'

-- |  Retrieve the bias of the top layer of the given deep harmonium.
getTopBias
    :: forall fs m ms c . ( Hierarchical fs (m : ms), Manifold m )
    => c # DeepHarmonium fs (m : ms)
    -> c # m
{-# INLINE getTopBias #-}
getTopBias dhrm =
    let (pmcs, _ :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double)
            = S.splitAt $ coordinates dhrm
     in Point pmcs

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs n ns c . ( HierarchicalTranspose fs (n : ns), Manifold n )
    => c # n
    -> c # HarmoniumTranspose fs (n : ns)
    -> c # HarmoniumTranspose fs (n : ns)
{-# INLINE biasBottom #-}
biasBottom pn' dhrm =
    let css' :: S.Vector (Dimension (HarmoniumTranspose fs (n : ns)) - Dimension n) Double
        (css',pncs) = S.splitAt $ coordinates dhrm
        pn = pn' <+> Point pncs
     in Point $ css' S.++ coordinates pn

-- |  Retrieve the bias of the top layer of the given deep harmonium.
getBottomBias
    :: forall fs n ns c . ( HierarchicalTranspose fs (n : ns), Manifold n )
    => c # HarmoniumTranspose fs (n : ns)
    -> c # n
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let (_ :: S.Vector (Dimension (HarmoniumTranspose fs (n : ns)) - Dimension n) Double, pncs)
            = S.splitAt $ coordinates dhrm
     in Point pncs

-- | The given deep harmonium conditioned on its top layer.
(>|>) :: ( Bilinear f m n, Hierarchical fs (n : ms) )
      => Mean # m
      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (>|>) #-}
(>|>) p dhrm =
    let (_,f,dhrm') = splitHeadHarmonium dhrm
     in biasTop (p <.< f) dhrm'

-- | The given deep harmonium conditioned on its bottom layer.
posterior :: ( Bilinear f m n, Map Mean Natural f m n, HierarchicalTranspose fs (m : ns) )
      => Mean # n
      -> Natural # HarmoniumTranspose (f : fs) (n : m : ns)
      -> Natural # HarmoniumTranspose fs (m : ns)
{-# INLINE posterior #-}
posterior p dhrm =
    let (dhrm',f,_) = splitLastHarmonium dhrm
     in biasBottom (f >.> p) dhrm'

(<|<) :: ( Hierarchical fsg msn
         , Hierarchical fs ms
         , Reverse3 fsg ~ (f : Reverse3 fs)
         , Reverse msn ~ (n : m : ns)
         , Reverse ms ~ (m : ns)
         , Bilinear f m n
         , Map Mean Natural f m n
         , Manifold (HarmoniumTranspose (Reverse3 fs) (m : ns)) )
  => Natural # DeepHarmonium fsg msn
  -> Mean # n
  -> Natural # DeepHarmonium fs ms
{-# INLINE (<|<) #-}
(<|<) dhrm p = detransposeHarmonium . posterior p $ transposeHarmonium dhrm


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled sequentially down/up the hierarchy.
-- Gibbs sampling is then the infinite, recursive application of these upward
-- and downward passes.

class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
    downwardPass :: KnownNat l
           => Natural # DeepHarmonium fs ms
           -> Sample l (DeepHarmonium fs ms)
           -> Random s (Sample l (DeepHarmonium fs ms))
--    downwardPass :: KnownNat l
--           => Natural # DeepHarmonium fs ms
--           -> Sample l (Head ms)
--           -> Random s (Sample l (DeepHarmonium fs ms))


---- | Sample
--upwardPass :: ( Hierarchical fs ms, KnownNat l )
--        => Natural # DeepHarmonium fs ms
--        -> Sample l (Last ms)
--        -> Random s (Sample l (DeepHarmonium fs ms))
--{-# INLINE upwardPass #-}
--upwardPass dhrm z = do
--    smps0 <- downwardPass (harmoniumTranspose dhrm) z
--    return $ hReverse <$> smps0


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


instance Manifold m => Manifold (OneHarmonium m) where
    type Dimension (OneHarmonium m) = Dimension m

instance Manifold n => Manifold (HarmoniumTranspose '[] '[n]) where
    type Dimension (HarmoniumTranspose '[] '[n]) = Dimension n

instance (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
  => Manifold (DeepHarmonium (f : fs) (m : n : ms)) where
      type Dimension (DeepHarmonium (f : fs) (m : n : ms))
        = Dimension m + Dimension (f m n) + Dimension (DeepHarmonium fs (n : ms))

instance (Bilinear f m n, Manifold (HarmoniumTranspose fs (m : ns)))
  => Manifold (HarmoniumTranspose (f : fs) (n : m : ns)) where
      type Dimension (HarmoniumTranspose (f : fs) (n : m : ns))
        = Dimension n + Dimension (f m n) + Dimension (HarmoniumTranspose fs (m : ns))

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
              zp = getTopBias dhrm'
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

--instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n
--         , Hierarchical fs (n : ms), Gibbs fs (n : ms)
--         , Hierarchical (f : fs) (m : n : ms) ) => Gibbs (f : fs) (m : n : ms) where
--      {-# INLINE downwardPass #-}
--      downwardPass dhrm xs = do
--          let ps = sufficientStatistic xs
--              (_,f,dhrm') = splitHeadHarmonium dhrm
--          ys <- samplePoint . mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
--          smps <- downwardPass dhrm' ys
--          return $ B.zipWith (:+:) xs smps
--instance Gibbs '[] '[m] where
--    {-# INLINE downwardPass #-}
--    downwardPass _ = return . fmap (:+: Null)
--
--instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n
--         , Hierarchical fs (n : ms), Gibbs fs (n : ms)
--         , Hierarchical (f : fs) (m : n : ms) ) => Gibbs (f : fs) (m : n : ms) where
--      {-# INLINE downwardPass #-}
--      downwardPass dhrm xs = do
--          let ps = sufficientStatistic xs
--              (_,f,dhrm') = splitHeadHarmonium dhrm
--          ys <- samplePoint . mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
--          smps <- downwardPass dhrm' ys
--          return $ B.zipWith (:+:) xs smps

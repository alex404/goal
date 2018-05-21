{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , type (<*>)
    , DeepHarmonium
    , Hierarchical
    -- ** Construction
    , fromOneHarmonium
    , toOneHarmonium
    , splitHeadHarmonium
    , joinHeadHarmonium
    -- ** Manipulation
    , biasTop
    , getTopBias
    , (>|>)
    , (<|<)
    , HarmoniumTranspose (harmoniumTranspose)
    -- ** Sampling
    , Gibbs (downwardPass)
    , upwardPass
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


-- | A hierarchical generative model defined by exponential families.
data DeepHarmonium (fs :: [* -> * -> *]) (ms :: [*])

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
    , HarmoniumTranspose fs ms
    , Reversing (SamplePoints (Reverse ms))
    , Reverse (SamplePoints (Reverse ms)) ~ SamplePoints ms
    , Gibbs fs ms
    , Gibbs (Reverse3 fs) (Reverse ms) )


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

-- | The given deep harmonium conditioned on its top layer.
(>|>) :: ( Bilinear f m n, Hierarchical fs (n : ms) )
      => Mean # m
      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Natural # DeepHarmonium fs (n : ms)
{-# INLINE (>|>) #-}
(>|>) p dhrm =
    let (_,f,dhrm') = splitHeadHarmonium dhrm
     in biasTop (p <.< f) dhrm'

-- | The given deep harmonium conditioned on its bottom layer. This will
-- typically correspond to the posterior in the context of Bayesian inference.
(<|<) :: ( Hierarchical fs' (n : ms')
         , Reverse ms ~ (m : n : ms')
         , Reverse3 fs ~ (f : fs')
         , HarmoniumTranspose fs ms
         , HarmoniumTranspose fs' (n : ms')
         , Bilinear f m n )
  => Natural # DeepHarmonium fs ms
  -> Mean # m
  -> Natural # DeepHarmonium (Reverse3 fs') (ReverseAcc ms' '[n])
{-# INLINE (<|<) #-}
(<|<) dhrm p = harmoniumTranspose $ p >|> harmoniumTranspose dhrm


--- Classes ---


-- | Harmonium transpotion. Each defining layers are reversed, and the defining
-- bilinear functions are transposed.
class Manifold (DeepHarmonium fs ms) => HarmoniumTranspose fs ms where
    harmoniumTranspose :: Primal c => c # DeepHarmonium fs ms -> c # DeepHarmonium (Reverse3 fs) (Reverse ms)

-- | 'Gibbs' deep harmoniums can be sampled sequentially down/up the hierarchy.
-- Gibbs sampling is then the infinite, recursive application of these upward
-- and downward passes.
--
-- Ooops! This is fucked!!
class Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
    downwardPass :: KnownNat l
           => Natural # DeepHarmonium fs ms
           -> Sample l (Head ms)
           -> Random s (Sample l (DeepHarmonium fs ms))

-- | Sample
upwardPass :: ( Hierarchical fs ms, KnownNat l )
        => Natural # DeepHarmonium fs ms
        -> Sample l (Last ms)
        -> Random s (Sample l (DeepHarmonium fs ms))
{-# INLINE upwardPass #-}
upwardPass dhrm z = do
    smps0 <- downwardPass (harmoniumTranspose dhrm) z
    return $ hReverse <$> smps0


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

instance Manifold m => HarmoniumTranspose '[] '[m] where
    {-# INLINE harmoniumTranspose #-}
    harmoniumTranspose = id

instance (Bilinear f m n, Bilinear f n m, HarmoniumTranspose fs (n : ms))
  => HarmoniumTranspose (f : fs) (m : n : ms) where
    {-# INLINE harmoniumTranspose #-}
    harmoniumTranspose dhrm =
        let (pm,pmtx,dhrm') = splitHeadHarmonium dhrm
            dhrm'' = harmoniumTranspose dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

instance Gibbs '[] '[m] where
    {-# INLINE downwardPass #-}
    downwardPass _ = return . fmap (:+: Null)

instance ( Bilinear f m n, ExponentialFamily m, Generative Natural n
         , Hierarchical fs (n : ms), Gibbs fs (n : ms)
         , Hierarchical (f : fs) (m : n : ms) ) => Gibbs (f : fs) (m : n : ms) where
      {-# INLINE downwardPass #-}
      downwardPass dhrm xs = do
          let ps = sufficientStatistic xs
              (_,f,dhrm') = splitHeadHarmonium dhrm
          ys <- samplePoint . mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
          smps <- downwardPass dhrm' ys
          return $ B.zipWith (:+:) xs smps

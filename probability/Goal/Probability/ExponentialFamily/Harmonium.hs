{-# LANGUAGE UndecidableInstances #-}
-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium where
    {-
    ( -- * Harmoniums
      Harmonium
    , type (<*>)
    -- ** Structure Manipulation
    , splitHarmonium
    , joinHarmonium
    , affineToHarmonium
    , harmoniumTranspose
    -- ** Conditional Distributions
    , conditionalLatentDistribution
    , conditionalObservableDistribution
     -- * Rectification
    , categoricalHarmoniumRectificationParameters
    , rectifiedBayesRule
    ) where
-}


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


data DeepHarmonium (fs :: [* -> * -> *]) (ms :: [*])

type OneHarmonium m = DeepHarmonium '[] '[m]

type Harmonium f m n = DeepHarmonium '[f] [m,n]

type (m <*> n) = DeepHarmonium '[Tensor] [m,n]

--type Hierarchical fs ms = ( Manifold (DeepHarmonium fs ms), Dimension (Head ms) <= Dimension (DeepHarmonium fs ms) )

type Hierarchical fs ms =
    ( Manifold (DeepHarmonium fs ms)
    , Manifold (DeepHarmonium (Init3 fs) (Init ms))
    , Manifold (DeepHarmonium (Tail3 fs) (Tail ms))
    , Dimension (Head ms) <= Dimension (DeepHarmonium fs ms)
    , Dimension (Last ms) <= Dimension (DeepHarmonium fs ms)
    , HarmoniumTranspose fs ms
    )


--- Functions ---


fromOneHarmonium :: c # DeepHarmonium '[] '[m] -> c # m
fromOneHarmonium = Point . coordinates

toOneHarmonium :: c # m -> c # DeepHarmonium '[] '[m]
toOneHarmonium = Point . coordinates


-- | Assembles a 'Harmonium' out of the component parameters.
joinHeadHarmonium
    :: (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
    => c # m
    -> Dual c ~> c # f m n
    -> c # DeepHarmonium fs (n : ms)
    -> c # DeepHarmonium (f : fs) (m : n : ms)
{-# INLINE joinHeadHarmonium #-}
joinHeadHarmonium pl pmtx dhrm =
    Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates dhrm

splitHeadHarmonium
    :: (Bilinear f m n, Manifold (DeepHarmonium fs (n : ms)))
    => c # DeepHarmonium (f : fs) (m : n : ms)
    -> (c # m, Dual c ~> c # f m n, c # DeepHarmonium fs (n : ms))
{-# INLINE splitHeadHarmonium #-}
splitHeadHarmonium dhrm =
    let (lcs,css') = S.splitAt $ coordinates dhrm
        (mtxcs,dcs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point dcs)

biasTop
    :: forall fs m ms c
    . ( Hierarchical fs (m : ms), Manifold m )
    => c # m
    -> c # DeepHarmonium fs (m : ms)
    -> c # DeepHarmonium fs (m : ms)
biasTop pm' dhrm =
    let css' :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double
        (pmcs,css') = S.splitAt $ coordinates dhrm
        pm = pm' <+> Point pmcs
     in Point $ coordinates pm S.++ css'

getTopBias
    :: forall fs m ms c
    . ( Hierarchical fs (m : ms), Manifold m )
    => c # DeepHarmonium fs (m : ms)
    -> c # m
getTopBias dhrm =
    let (pmcs, _ :: S.Vector (Dimension (DeepHarmonium fs (m : ms)) - Dimension m) Double)
            = S.splitAt $ coordinates dhrm
     in Point pmcs

(>|>) :: ( Bilinear f m n, Hierarchical fs (n : ms) )
      => Mean # m
      -> Natural # DeepHarmonium (f : fs) (m : n : ms)
      -> Natural # DeepHarmonium fs (n : ms)
(>|>) p dhrm =
    let (_,f,dhrm') = splitHeadHarmonium dhrm
     in biasTop (p <.< f) dhrm'

(<|<) :: ( Hierarchical fs' (n : ms')
         , Reverse ms ~ (m : n : ms')
         , Reverse3 fs ~ (f : fs')
         , HarmoniumTranspose fs ms
         , HarmoniumTranspose fs' (n : ms')
         , Bilinear f m n )
  => Natural # DeepHarmonium fs ms
  -> Mean # m
  -> Natural # DeepHarmonium (Reverse3 fs') (ReverseAcc ms' '[n])
(<|<) dhrm p = harmoniumTranspose $ p >|> harmoniumTranspose dhrm

--(<|<*) dhrm x = dhrm <|< sufficientStatistic x
--
--(*>|>) x dhrm = sufficientStatistic x >|> dhrm

--- Classes ---


class Manifold (DeepHarmonium fs ms) => HarmoniumTranspose fs ms where
    harmoniumTranspose :: Primal c => c # DeepHarmonium fs ms -> c # DeepHarmonium (Reverse3 fs) (Reverse ms)

class Hierarchical fs ms => Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
    (>|>*) :: (KnownNat k, KnownNat l, 1 <= l)
           => Mean # Replicated k (Head ms)
           -> Natural # DeepHarmonium fs ms
           -> Random s (Sample l (Replicated k (DeepHarmonium (Tail3 fs) (Tail ms))))

(*<|<) :: ( (1 <= l), KnownNat l, KnownNat k
          , HarmoniumTranspose fs ms
          , Gibbs (Reverse3 fs) (Reverse ms)
          , Reversing (SamplePoints (Tail (Reverse ms))) )
          => Natural # DeepHarmonium fs ms
          -> Mean # Replicated k (Last ms)
          -> Random s (B.Vector l (B.Vector k (HList (Reverse (SamplePoints (Tail (Reverse ms)))))))
(*<|<) dhrm p = do
    smps0 <- p >|>* harmoniumTranspose dhrm
    return $ fmap hReverse <$> smps0

--(*<|<*) dhrm x = dhrm *<|< sufficientStatistic x

--(*>|>*) x dhrm = sufficientStatistic x >|>* dhrm


--foo :: Natural # DeepHarmonium [Tensor,Tensor] [Bernoulli,Normal,Poisson]
--foo = zero
--
--bar :: Mean # Poisson -> Natural # DeepHarmonium '[Tensor] '[Bernoulli, Normal]
--bar q = foo <|< q
--
--baz :: Mean # Replicated 4 Poisson -> Random s (Sample 3 (Replicated 4 (DeepHarmonium '[Tensor] '[Bernoulli, Normal])))
--baz q = foo *<|< q


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

instance ExponentialFamily m => ExponentialFamily (OneHarmonium m) where
      sufficientStatistic (xl :+: Null) =
          toOneHarmonium $ sufficientStatistic xl
      baseMeasure = harmoniumBaseMeasure Proxy

instance ( ExponentialFamily m, ExponentialFamily n
         , Bilinear f m n, ExponentialFamily (DeepHarmonium fs (n : ms)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (m : n : ms)) where
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              mm = sufficientStatistic xm
              mh = sufficientStatistic xn
           in joinHeadHarmonium mm (mm >.< mh) mdhrm
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance Manifold m => HarmoniumTranspose '[] '[m] where
    harmoniumTranspose = id

instance (Bilinear f m n, Bilinear f n m, HarmoniumTranspose fs (n : ms))
  => HarmoniumTranspose (f : fs) (m : n : ms) where
    harmoniumTranspose dhrm =
        let (pm,pmtx,dhrm') = splitHeadHarmonium dhrm
            dhrm'' = harmoniumTranspose dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

instance ( Bilinear f m n, Generative Natural n, Hierarchical '[f] '[m,n] ) => Gibbs '[f] '[m,n] where
      (>|>*) ps hrm = do
          let (_,f,pz) = splitHeadHarmonium hrm
          smp <- sample $ mapReplicatedPoint (<+> fromOneHarmonium pz) $ ps <$< f
          return $ fmap (:+: Null) <$> smp

instance ( Bilinear f m n, ExponentialFamily n, Generative Natural n, Gibbs (g : fs) (n : o : ms)
         , Hierarchical (f : g : fs) (m : n : o : ms) ) => Gibbs (f : g : fs) (m : n : o : ms) where
      (>|>*) ps dhrm = do
          let (_,f,dhrm') = splitHeadHarmonium dhrm
          smp <- sample $ mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
          smps <- sufficientStatisticT smp >|>* dhrm'
          return $ B.zipWith (B.zipWith (:+:)) smp smps

--instance (Bilinear f m n, Bilinear f n m) => HarmoniumTranspose '[f] '[m,n] where
--    harmoniumTranspose hrm =
--        let (pm,pmtx,pn) = splitHeadHarmonium hrm
--         in joinHeadHarmonium (fromOneHarmonium pn) (transpose pmtx) (toOneHarmonium pm)
--
--splitLastHarmonium
--    :: (Manifold (DeepHarmonium (Init3 fs) (Init ms)), Bilinear (Last3 fs) (Last (Init ms)) (Last ms))
--    => c # DeepHarmonium fs ms
--    -> ( c # DeepHarmonium (Init3 fs) (Init ms)
--       , Dual c ~> c # Last3 fs (Last (Init ms)) (Last ms)
--       , c # Last ms )
--{-# INLINE splitLastHarmonium #-}
--splitLastHarmonium dhrm =
--    let v0 = I.Vector . S.fromSized $ coordinates dhrm
--        (dcs,css') = S.splitAt v0
--        (mtxcs,ocs) = S.splitAt css'
--     in (Point dcs, Point mtxcs, Point ocs)
--
--joinLastHarmonium
--    :: (Manifold (DeepHarmonium (Init3 fs) (Init ms)), Bilinear (Last3 fs) (Last (Init ms)) (Last ms))
--    => c # DeepHarmonium (Init3 fs) (Init ms)
--    -> Dual c ~> c # Last3 fs (Last (Init ms)) (Last ms)
--    -> c # Last ms
--    -> c # DeepHarmonium fs ms
--{-# INLINE joinLastHarmonium #-}
--joinLastHarmonium (Point dhrm) (Point lmtx) (Point lbs) =
--    Point . I.Vector . S.fromSized $ dhrm S.++ lmtx S.++ lbs
--
--biasBottom
--    :: forall fs ms c
--    . ( Hierarchical fs ms, Manifold (Last ms) )
--    => c # Last ms
--    -> c # DeepHarmonium fs ms
--    -> c # DeepHarmonium fs ms
--biasBottom pm' dhrm =
--    let css' :: S.Vector (Dimension (DeepHarmonium fs ms) - Dimension (Last ms)) Double
--        (css',pmcs) = S.splitAt $ coordinates dhrm
--        pm = pm' <+> Point pmcs
--     in Point $ css' S.++ coordinates pm
--
--getBottomBias
--    :: forall fs ms c
--    . ( Manifold (DeepHarmonium fs ms), Manifold (Last ms)
--      , Dimension (Last ms) <= Dimension (DeepHarmonium fs ms) )
--    => c # DeepHarmonium fs ms
--    -> c # Last ms
--getBottomBias dhrm =
--    let (_ :: S.Vector (Dimension (DeepHarmonium fs ms) - Dimension (Last ms)) Double,pmcs)
--            = S.splitAt $ coordinates dhrm
--     in Point pmcs
--
--
--instance ( fs ~ '[f1], ms ~ '[m1,m2], Bilinear f m (Head (m : ms)), Bilinear (Last3 (f : fs)) (Last (Init (m : ms))) (Last (m : ms))
--         , Map Mean Natural (Last3 (f : fs)) (Last (Init (m : ms))) (Last (m : ms))
--         , Gibbs fs ms, Manifold (DeepHarmonium (f : fs) (m : ms)) )
--  => Gibbs (f : fs) (m : ms) where
--
--foo :: forall k l s f1 f2 m1 m2 m3 . (KnownNat k, KnownNat l, Gibbs '[f2] [m2,m3], 1 <= l, Bilinear f1 m1 m2, Generative Natural m2, ExponentialFamily m2)
--       => Mean # Replicated k m1
--       -> Natural # DeepHarmonium [f1,f2] [m1,m2,m3]
---- -> Random s (Sample l (Replicated k m2), Sample l (Replicated k (DeepHarmonium '[] '[m3])))
--       -> Random s (Sample l (Replicated k (DeepHarmonium '[f2] [m2,m3])))
--foo ps dhrm = do
--    let (_,f,dhrm') = splitHeadHarmonium dhrm
--    (smp :: Sample l (Replicated k m2)) <- sample $ mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
--    (smps :: Sample l (Replicated k (DeepHarmonium '[] '[m3]))) <- sufficientStatisticT smp >|>* dhrm'
--    --return (smp,smps)
--    return $ B.zipWith (B.zipWith (:+:)) smp smps
--
--
--(<|<) :: ( Map Mean Natural (Last3 fs) (Last (Init ms)) (Last ms)
--         , Bilinear (Last3 fs) (Last (Init ms)) (Last ms)
--         , Hierarchical (Init3 fs) (Init ms) )
--      => Natural # DeepHarmonium fs ms
--      -> Mean # Last ms
--      -> Natural # DeepHarmonium (Init3 fs) (Init ms)
--(<|<) dhrm q =
--    let (dhrm',g,_) = splitLastHarmonium dhrm
--     in biasBottom (g >.> q) dhrm'
--
--(*>|>*) :: (KnownNat k, KnownNat l)
--       => Sample k (Head ms)
--       -> Natural # DeepHarmonium fs ms
--       -> Random s (Sample l (Replicated k (DeepHarmonium (Tail3 fs) (Tail ms))))
--
--(*<|<*) :: (KnownNat k, KnownNat l)
--       => Natural # DeepHarmonium fs ms
--       -> Sample k (Last ms)
--       -> Random s (Sample l (Replicated k (DeepHarmonium (Init3 fs) (Init ms))))
--
--
--class Hierarchical fs ms => Gibbs (fs :: [* -> * -> *]) (ms :: [*]) where
--
--    (>|>*) :: (KnownNat k, KnownNat l, 1 <= l)
--           => Mean # Replicated k (Head ms)
--           -> Natural # DeepHarmonium fs ms
--           -> Random s (Sample l (Replicated k (DeepHarmonium (Tail3 fs) (Tail ms))))
--
--    (*<|<) :: (KnownNat k, KnownNat l, 1 <= l)
--           => Natural # DeepHarmonium fs ms
--           -> Mean # Replicated k (Last ms)
--           -> Random s (Sample l (Replicated k (DeepHarmonium (Init3 fs) (Init ms))))
--

--instance ( Bilinear f m n, Map Mean Natural f m n, ExponentialFamily n, ExponentialFamily m
--         , Generative Natural n, Generative Natural m
--         , Hierarchical '[f] '[m,n] )
--  => Gibbs '[f] '[m,n] where
--
--      (>|>*) ps hrm = do
--          let (_,f,pz) = splitHeadHarmonium hrm
--          smp <- sample $ mapReplicatedPoint (<+> fromOneHarmonium pz) $ ps <$< f
--          return $ fmap (:+: Null) <$> smp
--      (*<|<) hrm qs = do
--          let (px,f,_) = splitLastHarmonium hrm
--          smp <- sample $ mapReplicatedPoint (<+> fromOneHarmonium px) $ f >$> qs
--          return $ fmap (:+: Null) <$> smp
--
--instance ( Bilinear f m n
--         , Bilinear (Last3 (g : fs)) (Last (Init (n : o : ms))) (Last ms)
--         , Map Mean Natural (Last3 (g : fs)) (Last (Init (m : n : o : ms))) (Last (m : n : o : ms))
--         , Gibbs (g : fs) (n : o : ms)
--         , Generative Natural n
--         , ExponentialFamily n
--         , Generative Natural (Last (m : Init (n : o : ms)))
--         , Map Mean Natural (Last3 (f : (g : fs))) (Last (m : Init (n : o : ms))) (Last (n : o : ms))
--         , Gibbs (Init3 (f : (g : fs))) (m : Init (n : o : ms))
--         , ExponentialFamily (Last (m : Init (n : o : ms)))
--         , Appending (SamplePoints (Init (m : Init (n : o : ms)))) (SamplePoint (Last (m : Init (n : o : ms))))
--         , Hierarchical (f : (g : fs)) (m : n : o : ms)
--         , Bilinear (Last3 (f : (g : fs))) (Last (m : Init (n : o : ms))) (Last (n : o : ms))
--         , (Dimension (Last (m : Init (n : o : ms))) <= Dimension (DeepHarmonium (Init3 (f : (g : fs))) (m : Init (n : o : ms))))
--         , Append (SamplePoints (Init (m : Init (n : o : ms)))) (SamplePoint (Last (m : Init (n : o : ms))))
--             ~ (SamplePoint m : SamplePoints (Init (n : o : ms))) )
--  => Gibbs (f : (g : fs)) (m : n : o : ms) where
--
--      (>|>*) ps dhrm = do
--          let (_,f,dhrm') = splitHeadHarmonium dhrm
--          smp <- sample $ mapReplicatedPoint (<+> getTopBias dhrm') $ ps <$< f
--          smps <- sufficientStatisticT smp >|>* dhrm'
--          return $ B.zipWith (B.zipWith (:+:)) smp smps
--      (*<|<) dhrm qs = do
--          let (dhrm',f,_) = splitLastHarmonium dhrm
--          smp <- sample $ mapReplicatedPoint (<+> getBottomBias dhrm') $ f >$> qs
--          smps <- dhrm' *<|< sufficientStatisticT smp
--          return $ B.zipWith (B.zipWith append) smps smp
--
--
--foo :: Natural # DeepHarmonium [Tensor,Tensor] [Bernoulli,Bernoulli,Bernoulli]
--foo = zero
--
--bar :: (KnownNat k, KnownNat l, 1 <= l) => Mean # Replicated k Bernoulli -> Random s (Sample l (Replicated k (DeepHarmonium '[Tensor] '[Bernoulli,Bernoulli])))
--bar q = foo *<|< q
--
--
---- | An exponential family defined using a product of two other exponential
---- families. The first argument represents the so-called observable variables, and
---- the second argument the so-called latent variables.
----type Harmonium f = DeepHarmonium '[f]
----
----type (<*>) m n = Harmonium (Product m n)
----
--
------ | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
------ | Assembles a 'Harmonium' out of the component parameters.
----joinHeadHarmonium
----    :: (Map f, Domain f ~ Codomain g)
----    => Point c (Codomain f)
----    -> Point (Function (Dual c) c) f -- ^ The component parameters
----    -> Point c (DeepHarmonium (g : fs))
----    -> Point c (DeepHarmonium (f : g : fs)) -- ^ The 'Harmonium'
----{-# INLINE joinHeadHarmonium #-}
----joinHeadHarmonium pl po (Point ps) =
----     Point $ coordinates pl S.++ coordinates po S.++ ps
----
----deepHarmoniumBaseMeasure3
----    :: (ExponentialFamily (Codomain f), ExponentialFamily (Codomain g), ExponentialFamily (Domain g))
----    => Proxy (Codomain f)
----    -> Proxy (Codomain g)
----    -> Proxy (Domain g)
----    -> Proxy (DeepHarmonium [f,g])
----    -> SamplePoint (DeepHarmonium [f,g])
----    -> Double
----{-# INLINE deepHarmoniumBaseMeasure3 #-}
----deepHarmoniumBaseMeasure3 prxyf prxyg1 prxyg2 _ (xs :+: ys1 :+: ys2 :+: Null) =
----     baseMeasure prxyf xs * baseMeasure prxyg1 ys1 * baseMeasure prxyg2 ys2
----
-------- | Returns the conditional distribution of the latent variables given the sufficient statistics of
-------- the observable state.
------conditionalLatentDistributions
------    :: (Apply Mean Natural f, ExponentialFamily (Domain f), KnownNat k)
------    => Point Natural (Harmonium f)
------    -> Sample k (Domain f)
------    -> Point Natural (Replicated k (Codomain f))
------{-# INLINE conditionalLatentDistributions #-}
------conditionalLatentDistributions hrm xs =
------    let (pl,f,_) = splitHarmonium hrm
------     in  mapReplicatedPoint (<+> pl) $ f >$>* xs
------
-------- | Returns the conditional distribution of the observable variables given the sufficient
-------- statistics of the latent state.
------conditionalObservableDistributions
------    :: (Bilinear f, ExponentialFamily (Codomain f), KnownNat k)
------    => Point Natural (Harmonium f)
------    -> Sample k (Codomain f)
------    -> Point Natural (Replicated k (Domain f))
------{-# INLINE conditionalObservableDistributions #-}
------conditionalObservableDistributions hrm xs =
------    let (_,f,po) = splitHarmonium hrm
------     in mapReplicatedPoint (<+> po) $ xs *<$< f
------
------conditionalLatentDistribution
------    :: (Bilinear Mean Natural f, ExponentialFamily (Domain f))
------    => Point Natural (Harmonium f)
------    -> SamplePoint (Domain f)
------    -> Point Natural (Codomain f)
------{-# INLINE conditionalLatentDistribution #-}
------conditionalLatentDistribution hrm x =
------    let (pl,_,f) = splitHarmonium hrm
------     in pl <+> (f >.>* x)
------
------
------conditionalObservableDistribution
------    :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f))
------    => Point Natural (Harmonium f)
------    -> SamplePoint (Codomain f)
------    -> Point Natural (Domain f)
------{-# INLINE conditionalObservableDistribution #-}
------conditionalObservableDistribution hrm x =
------    let (_,lb,f) = splitHarmonium hrm
------     in lb <+> (x *<.< f)
----
----{-
----conditionalObservableDistribution
----    :: (Manifold l, Manifold o)
----    => Point Natural (Harmonium l o)
----    -> Point (Mean ~> Natural) (Affine o l)
----conditionalObservableDistribution
----  = conditionalLatentDistribution . harmoniumTranspose
-----}
----
----{-
------ | Assembles a 'Harmonium' out of an 'Affine' transformation and latent biases.
----affineToHarmonium
----    :: (Primal c, Manifold l, Manifold o)
----    => Point c l -- ^ Latent biases
----    -> Point (Function (Dual c) c) (Affine o l) -- ^ Emission distribution
----    -> Point c (Harmonium l o) -- ^ The resulting 'Harmonium'
----affineToHarmonium pl paff =
----    let (po,pmtx) = splitAffine paff
----     in joinHarmonium pl po $ transpose pmtx
-----}
----
----{-
----
-----}
----
----
--- Instances ---


-- Harmoniums --

----instance (Map f, Map g, Manifold (DeepHarmonium (g : fs)), Domain f ~ Codomain g) => Manifold (DeepHarmonium (f : g : fs)) where
----    type Dimension (DeepHarmonium (f : g : fs)) =
----        Dimension (Codomain f) + Dimension f + Dimension (DeepHarmonium (g : fs))
----
----instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), ExponentialFamily (Domain g)
----         , Codomain g ~ Domain f, Bilinear f, Bilinear g)
----  => ExponentialFamily (DeepHarmonium '[f,g]) where
----      sufficientStatistic (xl :+: xm :+: xo :+: Null) =
----          let slcs = sufficientStatistic xl
----              smcs = sufficientStatistic xm
----              socs = sufficientStatistic xo
----           in joinHeadHarmonium slcs (slcs >.< smcs) $ joinHarmonium smcs (smcs >.< socs) socs
----      baseMeasure = deepHarmoniumBaseMeasure3 Proxy Proxy Proxy
----
------instance
------    ( Apply Mean Natural f, ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Bilinear f
------    , Generative Natural (Domain f), Generative Natural (Codomain f), Domain f ~ Codomain g)
------    => Hierarchical (f : g : fs) where
------      type TopLayer (f : g : fs) = Codomain f
------      type BottomLayer (f : g : fs) = Domain (Last (g : fs))
------      type TopSection (f : g : fs) = DeepHarmonium (Init (f : g : fs))
------      type BottomSection (f : g : fs) = DeepHarmonium (g : fs)
------      (>|>) ps dhrm =
------          let (_,f,Point cs) = splitHeadHarmonium dhrm
------              (qcs',dhrmcs) = S.splitAt cs
------              qcs'' = mapReplicated coordinates . mapReplicatedPoint (<+> Point qcs') $ ps <$< f
------           in joinReplicated $ S.map (Point . (S.++ dhrmcs)) qcs''
-------      (<|<) hrm qs =
-------        let (pl,f,_) = splitHarmonium hrm
-------         in mapReplicatedPoint (<+> pl) $ f >$> qs
-------      (>|>*) ps hrm = sample $ ps >|> hrm
-------      (*<|<) hrm qs = sample $ hrm <|< qs
----
----
----
----
------instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), ExponentialFamily (Domain g)
------         , Codomain g ~ Domain f, Bilinear f, Bilinear g)
------  => ExponentialFamily (DeepHarmonium '[f,g]) where
------      sufficientStatistic (xl :+: xm :+: xo :+: Null) =
------          let slcs = sufficientStatistic xl
------              smcs = sufficientStatistic xm
------              socs = sufficientStatistic xo
------           in joinHeadHarmonium slcs (slcs >.< smcs) $ joinHarmonium smcs (smcs >.< socs) socs
------      baseMeasure = deepHarmoniumBaseMeasure3 Proxy Proxy Proxy
----
----
------instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Codomain g ~ Domain f
------         , ExponentialFamily (DeepHarmonium (g : fs) ))
------  => ExponentialFamily (DeepHarmonium (f : g : fs)) where
------    sufficientStatistic (xl :+: xo :+: xs) =
------        let sls = sufficientStatistic xl
------            sos = sufficientStatistic xo
------            slos = slcs >.< socs
------            sdhrm = sufficientStatistic (xo :+: xs)
------         in joinDeepHarmonium sls slos sdhrm
------    baseMeasure = harmoniumBaseMeasure Proxy Proxy
----
----
------- Graveyard ---
----
----
----
----{-
----instance (Enum e, 1 <= n, KnownNat n, Generative Natural o, ClosedFormExponentialFamily o)
----    => Generative Natural (Harmonium (Categorical e n) o) where
----        generate hrm = S.head <$> sampleCategoricalHarmonium0 hrm
----
-----}
---- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
--splitHarmonium
--    :: Bilinear f m n
--    => c # Harmonium f m n -- ^ The 'Harmonium'
--    -> (c # m, Function (Dual c) c # f m n, c # n) -- ^ The component parameters
--{-# INLINE splitHarmonium #-}
--splitHarmonium plo =
--    let (lcs,css') = S.splitAt $ coordinates plo
--        (mtxcs,ocs) = S.splitAt css'
--     in (Point lcs, Point mtxcs, Point ocs)
--
---- | Assembles a 'Harmonium' out of the component parameters.
--joinHarmonium
--    :: Bilinear f m n
--    => c # m
--    -> Dual c ~> c # f m n -- ^ The component parameters
--    -> c # n
--    -> c # Harmonium f m n -- ^ The 'Harmonium'
--{-# INLINE joinHarmonium #-}
--joinHarmonium pl pmtx po =
--     Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates po
--


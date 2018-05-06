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

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Internal as I


--- Types ---


type Harmonium f m n = DeepHarmonium '[f] '[] m n

data DeepHarmonium (fs :: [* -> * -> *]) (hs :: [*]) m n

type family Append as b where
    Append '[] b = '[b]
    Append (a ': as) b = a ': Append as b

type family Append3 (as :: [* -> * -> *]) (b :: * -> * -> *) where
    Append3 '[] b = '[b]
    Append3 (a ': as) b = a ': Append3 as b

type family Init3 (as :: [* -> * -> *]) where
    Init3 '[a] = '[]
    Init3 (a ': as) = a ': Init3 as

type family Last3 (as :: [* -> * -> *]) where
    Last3 '[a] = a
    Last3 (a ': as) = Last3 as

type family Head3 (as :: [* -> * -> *]) where
    Head3 (a ': as) = a




--- De/construction ---


-- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
splitHarmonium
    :: Bilinear f m n
    => c # Harmonium f m n -- ^ The 'Harmonium'
    -> (c # m, Function (Dual c) c # f m n, c # n) -- ^ The component parameters
{-# INLINE splitHarmonium #-}
splitHarmonium plo =
    let (lcs,css') = S.splitAt $ coordinates plo
        (mtxcs,ocs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point ocs)

-- | Assembles a 'Harmonium' out of the component parameters.
joinHarmonium
    :: Bilinear f m n
    => c # m
    -> Dual c ~> c # f m n -- ^ The component parameters
    -> c # n
    -> c # Harmonium f m n -- ^ The 'Harmonium'
{-# INLINE joinHarmonium #-}
joinHarmonium pl pmtx po =
     Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates po

-- | Assembles a 'Harmonium' out of the component parameters.
joinHeadHarmonium
    :: (Bilinear f m h, Manifold (DeepHarmonium fs hs h n))
    => c # m
    -> Dual c ~> c # f m h -- ^ The component parameters
    -> c # DeepHarmonium fs hs h n
    -> c # DeepHarmonium (f : fs) (h : hs) m n -- ^ The 'Harmonium'
{-# INLINE joinHeadHarmonium #-}
joinHeadHarmonium pl pmtx (Point ps) =
     Point $ coordinates pl S.++ coordinates pmtx S.++ ps

splitLastHarmonium
    :: (Manifold (DeepHarmonium (Init3 fs) (Init hs) m (Last hs)), Bilinear (Last3 fs) (Last hs) n)
    => c # DeepHarmonium fs hs m n
    -> ( c # DeepHarmonium (Init3 fs) (Init hs) m (Last hs)
       , Dual c ~> c # Last3 fs (Last hs) n
       , c # n )
{-# INLINE splitLastHarmonium #-}
splitLastHarmonium plo =
    let v0 = I.Vector . S.fromSized $ coordinates plo
        (dcs,css') = S.splitAt v0
        (mtxcs,ocs) = S.splitAt css'
     in (Point dcs, Point mtxcs, Point ocs)

joinLastHarmonium
    :: (Manifold (DeepHarmonium (Init3 fs) (Init hs) m (Last hs)), Bilinear (Last3 fs) (Last hs) n)
    => c # DeepHarmonium (Init3 fs) (Init hs) m (Last hs)
    -> Dual c ~> c # Last3 fs (Last hs) n
    -> c # n
    -> c # DeepHarmonium fs hs m n
{-# INLINE joinLastHarmonium #-}
joinLastHarmonium (Point dhrm) (Point lmtx) (Point lbs) =
    Point . I.Vector . S.fromSized $ dhrm S.++ lmtx S.++ lbs

--joinLastHarmonium
--    :: (Bilinear (Last3 [f1,f2,f3]) (Last [m1,m2]) n, Manifold (DeepHarmonium (Init3 [f1,f2,f3]) (Init [m1,m2]) m (Last [m1,m2])))
--    => c # DeepHarmonium (Init3 [f1,f2,f3]) (Init [m1,m2]) m (Last [m1,m2])
--    -> Dual c ~> c # (Last3 [f1,f2,f3]) (Last [m1,m2]) n
--    -> c # n
--    -> c # DeepHarmonium [f1,f2,f3] [m1,m2] m n
--{-# INLINE joinLastHarmonium #-}
--joinLastHarmonium (Point dhrm) (Point lmtx) (Point lbs) =
--    Point $ dhrm S.++ lmtx S.++ lbs

-- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
splitHeadHarmonium
    :: (Bilinear f m h, Manifold (DeepHarmonium fs hs h n))
    => c # DeepHarmonium (f : fs) (h : hs) m n -- ^ The 'Harmonium'
    -> (c # m, Dual c ~> c # f m h, c # DeepHarmonium fs hs h n)
{-# INLINE splitHeadHarmonium #-}
splitHeadHarmonium plo =
    let (lcs,css') = S.splitAt $ coordinates plo
        (mtxcs,dcs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point dcs)

class Manifold (DeepHarmonium fs hs m n) => Hierarchical (fs :: [* -> * -> *]) (hs :: [*]) m n where

    type TopSection fs hs m n :: *
    type BottomSection fs hs m n :: *

    splitHeadHarmonium'
        :: c # DeepHarmonium fs hs m n
        -> (c # m, Dual c ~> c # (Head3 fs) m (Head (Append hs n)), BottomSection fs hs m n)
    joinHeadHarmonium'
        :: c # m
        -> Dual c ~> c # (Head3 fs) m (Head (Append hs n))
        -> BottomSection fs hs m n
        -> c # DeepHarmonium fs hs m n
    splitLastHarmonium'
        :: c # DeepHarmonium fs hs m n
        -> ( c # DeepHarmonium (Init3 fs) (Init hs) m (Last hs)
       , Dual c ~> c # Last3 fs (Last hs) n
       , c # n )

    (>|>*)
        :: (KnownNat k, KnownNat l)
        => Mean # Replicated k m
        -> Natural # DeepHarmonium fs hs m n
        -> Random s (Sample l (Replicated k (BottomSection fs hs m n)))
    (*<|<)
        :: (KnownNat k, KnownNat l)
        => Natural # DeepHarmonium fs hs m n
        -> Mean # Replicated k n
        -> Random s (Sample l (Replicated k (TopSection fs hs m n)))



--    biasTop :: c # m -> c # DeepHarmonium fs hs m n -> c # DeepHarmonium fs hs m n
--    biasBottom :: c # n -> c # DeepHarmonium fs hs m n -> c # DeepHarmonium fs hs m n

--    (>|>)
--        :: Mean # m
--        -> Natural # DeepHarmonium fs hs m n
--        -> Natural # BottomSection fs hs m n
--    (<|<)
--        :: Natural # DeepHarmonium fs hs m n
--        -> Mean # n
--        -> Natural # TopSection fs hs m n


--- Internal Functions ---


harmoniumBaseMeasure
    :: (ExponentialFamily m, ExponentialFamily n)
    => Proxy m
    -> Proxy n
    -> Proxy (Harmonium f m n)
    -> SamplePoint (Harmonium f m n)
    -> Double
{-# INLINE harmoniumBaseMeasure #-}
harmoniumBaseMeasure prxyl prxyo _ (ls :+: os :+: Null) =
     baseMeasure prxyl ls * baseMeasure prxyo os

deepHarmoniumBaseMeasure
    :: (ExponentialFamily m, ExponentialFamily (DeepHarmonium fs hs h n))
    => Proxy m
    -> Proxy (DeepHarmonium fs hs h n)
    -> Proxy (DeepHarmonium (f : fs) (h : hs) m n)
    -> SamplePoint (DeepHarmonium (f : fs) (h : hs) m n)
    -> Double
{-# INLINE deepHarmoniumBaseMeasure #-}
deepHarmoniumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs


----- Instances ---


instance Bilinear f m n => Manifold (Harmonium f m n) where
    type Dimension (Harmonium f m n) = Dimension m + Dimension (f m n) + Dimension n

instance (Bilinear f m h, Manifold (DeepHarmonium fs hs h n))
  => Manifold (DeepHarmonium (f : fs) (h : hs) m n) where
    type Dimension (DeepHarmonium (f : fs) (h : hs) m n)
      = Dimension m + Dimension (f m h) + Dimension (DeepHarmonium fs hs h n)

type family SamplePoints (ms :: [*]) where
    SamplePoints '[] = '[]
    SamplePoints (m : ms) = SamplePoint m : SamplePoints ms

instance Manifold (DeepHarmonium fs hs m n) => Statistical (DeepHarmonium fs hs m n) where
    type SamplePoint (DeepHarmonium fs hs m n) = HList (SamplePoints (Append (m : hs) n))

--instance (Bilinear f m h, Statistical (DeepHarmonium fs hs h n))
--  => Statistical (DeepHarmonium (f : fs) (h : hs) m n) where
--    type SamplePoint (DeepHarmonium (f : fs) (h : hs) m n)
--      = (SamplePoint m, SamplePoint (DeepHarmonium fs hs h n))
--
--instance Bilinear f m n => Statistical (Harmonium f m n) where
--    type SamplePoint (Harmonium f m n) = (SamplePoint m, SamplePoint n)

instance (ExponentialFamily m, ExponentialFamily n, Bilinear f m n)
  => ExponentialFamily (Harmonium f m n) where
      sufficientStatistic (xl :+: xo :+: Null) =
        let slcs = sufficientStatistic xl
            socs = sufficientStatistic xo
         in joinHarmonium slcs (slcs >.< socs) socs
      baseMeasure = harmoniumBaseMeasure Proxy Proxy

instance ( ExponentialFamily m, ExponentialFamily h, Bilinear f m h
         , ExponentialFamily (DeepHarmonium fs hs h n) )
  => ExponentialFamily (DeepHarmonium (f : fs) (h : hs) m n) where
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              mm = sufficientStatistic xm
              mh = sufficientStatistic xn
           in joinHeadHarmonium mm (mm >.< mh) mdhrm
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f m n, Map Mean Natural f m n, ExponentialFamily n, ExponentialFamily m
         , Generative Natural n, Generative Natural m )
  => Hierarchical '[f] '[] m n where
      type TopSection '[f] '[] m n = m
      type BottomSection '[f] '[] m n = n
--      biasTop pm' hrm =
--          let (pm,mtx,pn) = splitHarmonium hrm
--           in joinHarmonium (pm <+> pm') mtx pn
--      biasBottom pn' hrm =
--          let (pm,mtx,pn) = splitHarmonium hrm
--           in joinHarmonium pm mtx (pn <+> pn')
--      (>|>) p hrm =
--          let (_,f,po) = splitHarmonium hrm
--           in po <+> (p <.< f)
--      (<|<) hrm q =
--          let (pl,f,_) = splitHarmonium hrm
--           in pl <+> (f >.> q)
--      (>|>*) ps hrm =
--          let (_,f,po) = splitHarmonium hrm
--           in sample $ mapReplicatedPoint (<+> po) $ ps <$< f
--      (*<|<) hrm qs =
--          let (pl,f,_) = splitHarmonium hrm
--           in sample $ mapReplicatedPoint (<+> pl) $ f >$> qs

instance ( Bilinear f m h, Bilinear (Last3 (f : fs)) (Last (h : hs)) n
         , Map Mean Natural (Last3 (f : fs)) (Last (h : hs)) n
         , Hierarchical fs hs h n, Hierarchical (Init3 (f : fs)) (Init (h : hs)) m (Last (h : hs)) )
  => Hierarchical (f : fs) (h : hs) m n where

      type TopSection (f : fs) (h : hs) m n
        = DeepHarmonium (Init3 (f : fs)) (Init (h : hs)) m (Last (h : hs))
      type BottomSection (f : fs) (h : hs) m n = DeepHarmonium fs hs h n

--      biasTop pm' dhrm =
--          let (pm,mtx,dhrm') = splitHeadHarmonium dhrm
--           in joinHeadHarmonium (pm <+> pm') mtx dhrm'
--      biasBottom pn' dhrm =
--          let (dhrm',mtx,pn) = splitLastHarmonium dhrm
--           in joinLastHarmonium dhrm' mtx (pn <+> pn')
--
--      (>|>) p dhrm =
--          let (_,f,dhrm') = splitHeadHarmonium dhrm
--           in biasTop (p <.< f) dhrm'
--      (<|<) dhrm q =
--          let (dhrm',g,_) = splitLastHarmonium dhrm
--           in biasBottom (g >.> q) dhrm'
--
--      (>|>*) ps dhrm =
--          let (_,f,dhrm') = splitHeadHarmonium hrm
--           in sample $ mapReplicatedPoint (<+> po) $ ps <$< f

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
------ | Makes the observable variables the latent variables and vice versa by transposing the component
------ 'Tensor' and swapping the biases.
----harmoniumTranspose
----    :: (Manifold l, Manifold o, Primal c)
----    => Point c (Harmonium l o)
----    -> Point c (Harmonium o l)
----harmoniumTranspose plo =
----    let (pl,po,pmtx) = splitHarmonium plo
----     in joinHarmonium po pl (transpose pmtx)
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

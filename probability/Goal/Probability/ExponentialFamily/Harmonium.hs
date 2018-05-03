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


--- Types ---


data Harmonium f

data HarmoniumLayer f g



--- De/construction ---


-- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
splitHarmonium
    :: Map f
    => c # Harmonium f -- ^ The 'Harmonium'
    -> (c # Codomain f, Function (Dual c) c # f, c # Domain f) -- ^ The component parameters
{-# INLINE splitHarmonium #-}
splitHarmonium plo =
    let (lcs,css') = S.splitAt $ coordinates plo
        (mtxcs,ocs) = S.splitAt css'
     in (Point lcs, Point mtxcs, Point ocs)

-- | Assembles a 'Harmonium' out of the component parameters.
joinHarmonium
    :: Map f
    => Point c (Codomain f)
    -> Point (Function (Dual c) c) f -- ^ The component parameters
    -> Point c (Domain f)
    -> Point c (Harmonium f) -- ^ The 'Harmonium'
{-# INLINE joinHarmonium #-}
joinHarmonium pl pmtx po =
     Point $ coordinates pl S.++ coordinates pmtx S.++ coordinates po


--- Internal Functions ---


harmoniumBaseMeasure
    :: (ExponentialFamily (Codomain f), ExponentialFamily (Domain f))
    => Proxy (Codomain f)
    -> Proxy (Domain f)
    -> Proxy (Harmonium f)
    -> SamplePoint (Harmonium f)
    -> Double
{-# INLINE harmoniumBaseMeasure #-}
harmoniumBaseMeasure prxyl prxyo _ (ls,os) =
     baseMeasure prxyl ls * baseMeasure prxyo os


--- Instances ---


instance (Map f) => Manifold (Harmonium f) where
    type Dimension (Harmonium f) = Dimension (Codomain f) + Dimension (Domain f) + Dimension f

instance Map f => Statistical (Harmonium f) where
    type SamplePoint (Harmonium f) = (SamplePoint (Codomain f), SamplePoint (Domain f))

instance (ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Bilinear f)
  => ExponentialFamily (Harmonium f) where
    sufficientStatistic (xl,xo) =
        let slcs = sufficientStatistic xl
            socs = sufficientStatistic xo
         in joinHarmonium slcs (slcs >.< socs) socs
    baseMeasure = harmoniumBaseMeasure Proxy Proxy

-- | An exponential family defined using a product of two other exponential
-- families. The first argument represents the so-called observable variables, and
-- the second argument the so-called latent variables.
--type Harmonium f = DeepHarmonium '[f]
--
--type (<*>) m n = Harmonium (Product m n)
--

--class Hierarchical fs where
--    type TopLayer fs :: *
--    type BottomLayer fs :: *
--    type TopSection fs :: *
--    type BottomSection fs :: *
--    (>|>)
--        :: KnownNat k
--        => Mean # Replicated k (TopLayer fs)
--        -> Natural # DeepHarmonium fs
--        -> Natural # Replicated k (BottomSection fs)
--    (<|<)
--        :: KnownNat k
--        => Natural # DeepHarmonium fs
--        -> Mean # Replicated k (BottomLayer fs)
--        -> Natural # Replicated k (TopSection fs)
--    (>|>*)
--        :: (KnownNat k, KnownNat l)
--        => Mean # Replicated k (TopLayer fs)
--        -> Natural # DeepHarmonium fs
--        -> Random s (Sample l (Replicated k (BottomSection fs)))
--    (*<|<)
--        :: (KnownNat k, KnownNat l)
--        => Natural # DeepHarmonium fs
--        -> Mean # Replicated k (BottomLayer fs)
--        -> Random s (Sample l (Replicated k (TopSection fs)))
--
---- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
--splitHeadHarmonium
--    :: (Map f, Manifold (DeepHarmonium (g : fs)))
--    => c # DeepHarmonium (f : g : fs) -- ^ The 'Harmonium'
--    -> (c # Codomain f, Function (Dual c) c # f, c # DeepHarmonium (g : fs)) -- ^ The component parameters
--{-# INLINE splitHeadHarmonium #-}
--splitHeadHarmonium plo =
--    let (lcs,css') = S.splitAt $ coordinates plo
--        (mtxcs,dcs) = S.splitAt css'
--     in (Point lcs, Point mtxcs, Point dcs)
--
---- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
--splitLastHarmonium
--    :: (Manifold (DeepHarmonium (Init [f,g,h])), Map (Last [f,g,h]), Codomain h ~ Domain g)
--    => c # DeepHarmonium [f,g,h] -- ^ The 'Harmonium'
--    -> ( c # DeepHarmonium (Init [f,g,h])
--       , Function (Dual c) c # Last [f,g,h]
--       , c # Domain (Last [f,g,h]) ) -- ^ The component parameters
--{-# INLINE splitLastHarmonium #-}
--splitLastHarmonium plo =
--    let v0 = coordinates plo
--        (dcs,css') = S.splitAt v0
--        (mtxcs,ocs) = S.splitAt css'
--     in (Point dcs, Point mtxcs, Point ocs)
--
---- | Assembles a 'Harmonium' out of the component parameters.
--joinHeadHarmonium
--    :: (Map f, Domain f ~ Codomain g)
--    => Point c (Codomain f)
--    -> Point (Function (Dual c) c) f -- ^ The component parameters
--    -> Point c (DeepHarmonium (g : fs))
--    -> Point c (DeepHarmonium (f : g : fs)) -- ^ The 'Harmonium'
--{-# INLINE joinHeadHarmonium #-}
--joinHeadHarmonium pl po (Point ps) =
--     Point $ coordinates pl S.++ coordinates po S.++ ps
--
---- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
--joinLastHarmonium
--    :: (Manifold (DeepHarmonium [f,g,h,s]), Map (Last [f,g,h,s]), Codomain (Last [f,g,h,s]) ~ Domain (Last (Init [f,g,h,s])))
--    => c # DeepHarmonium (Init [f,g,h,s])
--    -> Function (Dual c) c # Last [f,g,h,s]
--    -> c # Domain (Last [f,g,h,s])
--    -> c # DeepHarmonium [f,g,h,s]
--{-# INLINE joinLastHarmonium #-}
--joinLastHarmonium (Point dhrm) (Point lmtx) (Point lbs) =
--    Point $ dhrm S.++ lmtx S.++ lbs
--
--deepHarmoniumBaseMeasure3
--    :: (ExponentialFamily (Codomain f), ExponentialFamily (Codomain g), ExponentialFamily (Domain g))
--    => Proxy (Codomain f)
--    -> Proxy (Codomain g)
--    -> Proxy (Domain g)
--    -> Proxy (DeepHarmonium [f,g])
--    -> SamplePoint (DeepHarmonium [f,g])
--    -> Double
--{-# INLINE deepHarmoniumBaseMeasure3 #-}
--deepHarmoniumBaseMeasure3 prxyf prxyg1 prxyg2 _ (xs :+: ys1 :+: ys2 :+: Null) =
--     baseMeasure prxyf xs * baseMeasure prxyg1 ys1 * baseMeasure prxyg2 ys2
--
------ | Returns the conditional distribution of the latent variables given the sufficient statistics of
------ the observable state.
----conditionalLatentDistributions
----    :: (Apply Mean Natural f, ExponentialFamily (Domain f), KnownNat k)
----    => Point Natural (Harmonium f)
----    -> Sample k (Domain f)
----    -> Point Natural (Replicated k (Codomain f))
----{-# INLINE conditionalLatentDistributions #-}
----conditionalLatentDistributions hrm xs =
----    let (pl,f,_) = splitHarmonium hrm
----     in  mapReplicatedPoint (<+> pl) $ f >$>* xs
----
------ | Returns the conditional distribution of the observable variables given the sufficient
------ statistics of the latent state.
----conditionalObservableDistributions
----    :: (Bilinear f, ExponentialFamily (Codomain f), KnownNat k)
----    => Point Natural (Harmonium f)
----    -> Sample k (Codomain f)
----    -> Point Natural (Replicated k (Domain f))
----{-# INLINE conditionalObservableDistributions #-}
----conditionalObservableDistributions hrm xs =
----    let (_,f,po) = splitHarmonium hrm
----     in mapReplicatedPoint (<+> po) $ xs *<$< f
----
----conditionalLatentDistribution
----    :: (Bilinear Mean Natural f, ExponentialFamily (Domain f))
----    => Point Natural (Harmonium f)
----    -> SamplePoint (Domain f)
----    -> Point Natural (Codomain f)
----{-# INLINE conditionalLatentDistribution #-}
----conditionalLatentDistribution hrm x =
----    let (pl,_,f) = splitHarmonium hrm
----     in pl <+> (f >.>* x)
----
----
----conditionalObservableDistribution
----    :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f))
----    => Point Natural (Harmonium f)
----    -> SamplePoint (Codomain f)
----    -> Point Natural (Domain f)
----{-# INLINE conditionalObservableDistribution #-}
----conditionalObservableDistribution hrm x =
----    let (_,lb,f) = splitHarmonium hrm
----     in lb <+> (x *<.< f)
--
--{-
--conditionalObservableDistribution
--    :: (Manifold l, Manifold o)
--    => Point Natural (Harmonium l o)
--    -> Point (Mean ~> Natural) (Affine o l)
--conditionalObservableDistribution
--  = conditionalLatentDistribution . harmoniumTranspose
---}
--
--{-
---- | Assembles a 'Harmonium' out of an 'Affine' transformation and latent biases.
--affineToHarmonium
--    :: (Primal c, Manifold l, Manifold o)
--    => Point c l -- ^ Latent biases
--    -> Point (Function (Dual c) c) (Affine o l) -- ^ Emission distribution
--    -> Point c (Harmonium l o) -- ^ The resulting 'Harmonium'
--affineToHarmonium pl paff =
--    let (po,pmtx) = splitAffine paff
--     in joinHarmonium pl po $ transpose pmtx
---}
--
--{-
---- | Makes the observable variables the latent variables and vice versa by transposing the component
---- 'Tensor' and swapping the biases.
--harmoniumTranspose
--    :: (Manifold l, Manifold o, Primal c)
--    => Point c (Harmonium l o)
--    -> Point c (Harmonium o l)
--harmoniumTranspose plo =
--    let (pl,po,pmtx) = splitHarmonium plo
--     in joinHarmonium po pl (transpose pmtx)
--
---}
--
--
----- Instances ---
--
--
---- Harmoniums --
--
--instance Manifold (DeepHarmonium '[]) where
--    type Dimension (DeepHarmonium '[]) = 0
--
--instance Map f => Manifold (DeepHarmonium '[f]) where
--    type Dimension (DeepHarmonium '[f]) = Dimension f + Dimension (Codomain f) + Dimension (Domain f)
--
--instance (Map f, Map g, Manifold (DeepHarmonium (g : fs)), Domain f ~ Codomain g) => Manifold (DeepHarmonium (f : g : fs)) where
--    type Dimension (DeepHarmonium (f : g : fs)) =
--        Dimension (Codomain f) + Dimension f + Dimension (DeepHarmonium (g : fs))
--
--instance Manifold (DeepHarmonium fs) => Statistical (DeepHarmonium fs) where
--    type SamplePoint (DeepHarmonium fs) = HList (SamplePoints fs)
--
--type family SamplePoints (ms :: [*]) where
--    SamplePoints '[f] = '[SamplePoint (Codomain f), SamplePoint (Domain f)]
--    SamplePoints (f : fs) = SamplePoint (Codomain f) : SamplePoints fs
--
--instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), ExponentialFamily (Domain g)
--         , Codomain g ~ Domain f, Bilinear f, Bilinear g)
--  => ExponentialFamily (DeepHarmonium '[f,g]) where
--      sufficientStatistic (xl :+: xm :+: xo :+: Null) =
--          let slcs = sufficientStatistic xl
--              smcs = sufficientStatistic xm
--              socs = sufficientStatistic xo
--           in joinHeadHarmonium slcs (slcs >.< smcs) $ joinHarmonium smcs (smcs >.< socs) socs
--      baseMeasure = deepHarmoniumBaseMeasure3 Proxy Proxy Proxy
--
--instance ( Apply Mean Natural f, ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Bilinear f
--         , Generative Natural (Domain f), Generative Natural (Codomain f) )
--  => Hierarchical ('[f]) where
--      type TopLayer ('[f]) = Codomain f
--      type BottomLayer ('[f]) = Domain f
--      type TopSection ('[f]) = Codomain f
--      type BottomSection ('[f]) = Domain f
--      (>|>) ps hrm =
--          let (_,f,po) = splitHarmonium hrm
--           in mapReplicatedPoint (<+> po) $ ps <$< f
--      (<|<) hrm qs =
--          let (pl,f,_) = splitHarmonium hrm
--           in mapReplicatedPoint (<+> pl) $ f >$> qs
--      (>|>*) ps hrm = sample $ ps >|> hrm
--      (*<|<) hrm qs = sample $ hrm <|< qs
--
----instance
----    ( Apply Mean Natural f, ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Bilinear f
----    , Generative Natural (Domain f), Generative Natural (Codomain f), Domain f ~ Codomain g)
----    => Hierarchical (f : g : fs) where
----      type TopLayer (f : g : fs) = Codomain f
----      type BottomLayer (f : g : fs) = Domain (Last (g : fs))
----      type TopSection (f : g : fs) = DeepHarmonium (Init (f : g : fs))
----      type BottomSection (f : g : fs) = DeepHarmonium (g : fs)
----      (>|>) ps dhrm =
----          let (_,f,Point cs) = splitHeadHarmonium dhrm
----              (qcs',dhrmcs) = S.splitAt cs
----              qcs'' = mapReplicated coordinates . mapReplicatedPoint (<+> Point qcs') $ ps <$< f
----           in joinReplicated $ S.map (Point . (S.++ dhrmcs)) qcs''
-----      (<|<) hrm qs =
-----        let (pl,f,_) = splitHarmonium hrm
-----         in mapReplicatedPoint (<+> pl) $ f >$> qs
-----      (>|>*) ps hrm = sample $ ps >|> hrm
-----      (*<|<) hrm qs = sample $ hrm <|< qs
--
--
--
--
----instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), ExponentialFamily (Domain g)
----         , Codomain g ~ Domain f, Bilinear f, Bilinear g)
----  => ExponentialFamily (DeepHarmonium '[f,g]) where
----      sufficientStatistic (xl :+: xm :+: xo :+: Null) =
----          let slcs = sufficientStatistic xl
----              smcs = sufficientStatistic xm
----              socs = sufficientStatistic xo
----           in joinHeadHarmonium slcs (slcs >.< smcs) $ joinHarmonium smcs (smcs >.< socs) socs
----      baseMeasure = deepHarmoniumBaseMeasure3 Proxy Proxy Proxy
--
--
----instance ( ExponentialFamily (Domain f), ExponentialFamily (Codomain f), Codomain g ~ Domain f
----         , ExponentialFamily (DeepHarmonium (g : fs) ))
----  => ExponentialFamily (DeepHarmonium (f : g : fs)) where
----    sufficientStatistic (xl :+: xo :+: xs) =
----        let sls = sufficientStatistic xl
----            sos = sufficientStatistic xo
----            slos = slcs >.< socs
----            sdhrm = sufficientStatistic (xo :+: xs)
----         in joinDeepHarmonium sls slos sdhrm
----    baseMeasure = harmoniumBaseMeasure Proxy Proxy
--
--
----- Graveyard ---
--
--
--
--{-
--instance (Enum e, 1 <= n, KnownNat n, Generative Natural o, ClosedFormExponentialFamily o)
--    => Generative Natural (Harmonium (Categorical e n) o) where
--        generate hrm = S.head <$> sampleCategoricalHarmonium0 hrm
--
---}

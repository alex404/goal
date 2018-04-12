{-# LANGUAGE UndecidableInstances,ExplicitNamespaces #-}
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

-- | An exponential family defined using a product of two other exponential
-- families. The first argument represents the so-called observable variables, and
-- the second argument the so-called latent variables.
data Harmonium f

type H f = Harmonium f

type (<*>) m n = Harmonium (Product m n)

-- | Splits a 'Harmonium' into its components parts of a pair of biases and a 'Tensor'.
splitHarmonium
    :: Map f
    => Point c (Harmonium f) -- ^ The 'Harmonium'
    -> (Point c (Codomain f), Point c (Domain f), Point (Function (Dual c) c) f) -- ^ The component parameters
{-# INLINE splitHarmonium #-}
splitHarmonium plo =
    let (lcs,css') = S.splitAt $ coordinates plo
        (ocs,mtxcs) = S.splitAt css'
     in (Point lcs, Point ocs, Point mtxcs)

-- | Assembles a 'Harmonium' out of the component parameters.
joinHarmonium
    :: Map f
    => Point c (Codomain f)
    -> Point c (Domain f)
    -> Point (Function (Dual c) c) f -- ^ The component parameters
    -> Point c (Harmonium f) -- ^ The 'Harmonium'
{-# INLINE joinHarmonium #-}
joinHarmonium pl po pmtx =
     Point $ coordinates pl S.++ coordinates po S.++ coordinates pmtx

harmoniumBaseMeasure0
    :: (ExponentialFamily (Codomain f), ExponentialFamily (Domain f))
    => Proxy (Codomain f) -> Proxy (Domain f) -> Proxy (Harmonium f) -> SamplePoint (Harmonium f) -> Double
{-# INLINE harmoniumBaseMeasure0 #-}
harmoniumBaseMeasure0 prxyl prxyo _ (ls,os) =
     baseMeasure prxyl ls * baseMeasure prxyo os

-- | Returns the conditional distribution of the latent variables given the sufficient statistics of
-- the observable state.
conditionalLatentDistributions
    :: (Bilinear Mean Natural f, ExponentialFamily (Domain f), KnownNat k)
    => Point Natural (Harmonium f)
    -> Sample k (Domain f)
    -> Point Natural (Replicated k (Codomain f))
{-# INLINE conditionalLatentDistributions #-}
conditionalLatentDistributions hrm xs =
    let (pl,_,f) = splitHarmonium hrm
     in  mapReplicatedPoint (<+> pl) $ f >$>* xs

conditionalLatentDistribution
    :: (Bilinear Mean Natural f, ExponentialFamily (Domain f))
    => Point Natural (Harmonium f)
    -> SamplePoint (Domain f)
    -> Point Natural (Codomain f)
{-# INLINE conditionalLatentDistribution #-}
conditionalLatentDistribution hrm x =
    let (pl,_,f) = splitHarmonium hrm
     in pl <+> (f >.>* x)

-- | Returns the conditional distribution of the observable variables given the sufficient
-- statistics of the latent state.
conditionalObservableDistributions
    :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f), KnownNat k)
    => Point Natural (Harmonium f)
    -> Sample k (Codomain f)
    -> Point Natural (Replicated k (Domain f))
{-# INLINE conditionalObservableDistributions #-}
conditionalObservableDistributions hrm xs =
    let (_,lb,f) = splitHarmonium hrm
     in mapReplicatedPoint (<+> lb) $ xs *<$< f

conditionalObservableDistribution
    :: (Bilinear Mean Natural f, ExponentialFamily (Codomain f))
    => Point Natural (Harmonium f)
    -> SamplePoint (Codomain f)
    -> Point Natural (Domain f)
{-# INLINE conditionalObservableDistribution #-}
conditionalObservableDistribution hrm x =
    let (_,lb,f) = splitHarmonium hrm
     in lb <+> (x *<.< f)

{-
conditionalObservableDistribution
    :: (Manifold l, Manifold o)
    => Point Natural (Harmonium l o)
    -> Point (Mean ~> Natural) (Affine o l)
conditionalObservableDistribution
  = conditionalLatentDistribution . harmoniumTranspose
-}

{-
-- | Assembles a 'Harmonium' out of an 'Affine' transformation and latent biases.
affineToHarmonium
    :: (Primal c, Manifold l, Manifold o)
    => Point c l -- ^ Latent biases
    -> Point (Function (Dual c) c) (Affine o l) -- ^ Emission distribution
    -> Point c (Harmonium l o) -- ^ The resulting 'Harmonium'
affineToHarmonium pl paff =
    let (po,pmtx) = splitAffine paff
     in joinHarmonium pl po $ transpose pmtx
-}

{-
-- | Makes the observable variables the latent variables and vice versa by transposing the component
-- 'Tensor' and swapping the biases.
harmoniumTranspose
    :: (Manifold l, Manifold o, Primal c)
    => Point c (Harmonium l o)
    -> Point c (Harmonium o l)
harmoniumTranspose plo =
    let (pl,po,pmtx) = splitHarmonium plo
     in joinHarmonium po pl (transpose pmtx)


--- Rectification ---


sampleRectifiedHarmonium
    :: (ExponentialFamily l, Generative Natural o, Generative Natural l, KnownNat k)
    => Point Natural l x -- ^ Rectification Parameters
    -> Point Natural (Harmonium l o) x -- ^ Rectified Harmonium Parameters
    -> Random s (S.Vector k (Sample l, Sample o)) -- ^ Samples
sampleRectifiedHarmonium rx hrm = do
    xs <- replicateMV $ generate rx
    zs <- mapM generate $ conditionalObservableDistribution hrm >$>* xs
    return $ zipV xs zs

categoricalHarmoniumRectificationParameters
    :: (Enum e, 1 <= n, KnownNat n, ClosedFormExponentialFamily o)
    => Point Natural (Harmonium (Categorical e n) o) x
    -> (x, Point Natural (Categorical e n) x)
categoricalHarmoniumRectificationParameters hrm =
    let (nx,nz,nxz) = splitHarmonium hrm
        rho0 = potential nz
        generator nxzi = subtract rho0 . potential $ nz <+> Point nxzi
        rprms = Point $ generator <$> toRows (toMatrix nxz)
     in (rho0,nx <+> rprms)

-- | Bayes' rule given a rectified harmonium generative model.
rectifiedBayesRule
    :: (ExponentialFamily l, ExponentialFamily o)
    => Point (Function Mean Natural) (Affine o l) x -- ^ Likelihood
    -> Point Natural l x
    -> Sample o -- ^ Observation
    -> Point Natural l x -- ^ Prior
    -> Point Natural l x -- ^ Posterior
rectifiedBayesRule lklhd rprms z p0 =
    let mtx = transpose . snd $ splitAffine lklhd
     in mtx >.>* z <+> p0 <+> rprms

sampleCategoricalHarmonium0
    :: (Enum e, 1 <= n, KnownNat n, Generative Natural o, ClosedFormExponentialFamily o)
    => Point Natural (Harmonium (Categorical e n) o) x
    -> Random s (S.Vector 1 (e, Sample o))
sampleCategoricalHarmonium0 hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    sampleRectifiedHarmonium rx hrm


-}


--- Instances ---


-- Harmoniums --

instance (Map f) => Manifold (Harmonium f) where
    type Dimension (Harmonium f) = Dimension (Codomain f) + Dimension (Domain f) + Dimension f

instance Map f => Statistical (Harmonium f) where
    type SamplePoint (Harmonium f) = (SamplePoint (Codomain f), SamplePoint (Domain f))

instance (ExponentialFamily l, ExponentialFamily o) => ExponentialFamily (l <*> o) where
    sufficientStatistic (xl,xo) =
        let slcs = sufficientStatistic xl
            socs = sufficientStatistic xo
         in joinHarmonium slcs socs (slcs >.< socs)
    baseMeasure = harmoniumBaseMeasure0 Proxy Proxy




--- Graveyard ---



{-
instance (Enum e, 1 <= n, KnownNat n, Generative Natural o, ClosedFormExponentialFamily o)
    => Generative Natural (Harmonium (Categorical e n) o) where
        generate hrm = S.head <$> sampleCategoricalHarmonium0 hrm

-}
-- Datatype manipulation --

{-
{-
-- | The gradient for training a rectified harmonium.
rectificationDifferentials
    :: (SourceGenerative Natural l, ExponentialFamily l)
    => [[Int]] -- ^ Stochastic Gradient Sample
    -> Natural :#: Harmonium l (Replicated Poisson)
    -> RandST r (Differentials :#: Tangent Natural (Harmonium l (Replicated Poisson)))
rectificationDifferentials ns p = do
    xs <- mapM standardGenerate $ conditionalLatentDistribution p >$>* ns
    let pos = conditionalObservableDistribution p >$>* xs
        Harmonium l o = manifold p
        mls = sufficientStatistic l <$> xs
        mos = sufficientStatistic o <$> ns
        c = sum . listCoordinates $ averagePoint mos
        cs' = subtract c . potential <$> pos
        diprms = [ ml >.< mo | (ml,mo) <- zip mls mos ]
        dps = [ (c' .>) . fromCoordinates (Tangent p)
            $ coordinates ml C.++ C.replicate (dimension o) 0 C.++ coordinates diprm
            | (c',ml,diprm) <- zip3 cs' mls diprms ]
    return $ averagePoint dps

     -}

harmoniumGradientCalculator
    :: [Mixture :#: x]
    -> [Mixture :#: x]
    -> [Mixture :#: z]
    -> [Mixture :#: z]
    -> Natural :#: Harmonium x z
    -> Differentials :#: Tangent Natural (Harmonium x z)
harmoniumGradientCalculator mxs mxs' mzs mzs' hrm =
    let dxs = zipWith (<->) mxs mxs'
        dzs = zipWith (<->) mzs mzs'
        dis = [ (mx >.< mz) <-> (mx' >.< mz') | (mx,mx',mz,mz') <- zip4 mxs mxs' mzs mzs' ]
     in averagePoint [ fromCoordinates (Tangent hrm)
        $ coordinates dx C.++ coordinates dz C.++ coordinates di | (dx,dz,di) <- zip3 dxs dzs dis ]

-}

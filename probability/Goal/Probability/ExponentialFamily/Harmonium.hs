{-# LANGUAGE TypeApplications,UndecidableInstances #-}
-- | Exponential Family Harmoniums and Conjugation.
module Goal.Probability.ExponentialFamily.Harmonium
    ( -- * Harmoniums
      OneHarmonium
    , Harmonium
    , Mixture
    , DeepHarmonium
    , unnormalizedHarmoniumObservableDensity
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    , fromNullMixture
    , toNullMixture
    -- ** Construction
    , biasBottom
    , getBottomBias
    , splitBottomHarmonium
    , joinBottomHarmonium
    -- ** Transposition
    , TransposeHarmonium (transposeHarmonium)
    -- ** Sampling
    , Gibbs (upwardPass,initialPass)
    , gibbsPass
    , harmoniumEmpiricalExpectations
    , harmoniumEmpiricalExpectations0
    -- * Conjugated Harmoniums
    , conjugatedHarmoniumDensity
    , logConjugatedHarmoniumDensity
    , marginalizeConjugatedHarmonium
    , SampleConjugated (sampleConjugatedHarmonium)
    -- * Mixture Models
    , joinMixture
    , splitMixture
    , mixtureDensity
    , logMixtureDensity
    , mixtureParameters
    , mixtureExpectations
    , mixtureRelativeEntropyUpperBound
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium.Conjugation
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Internal as I
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H


--- Types ---


-- | A hierarchical generative model defined by exponential families. Note that
-- the first elements of ms is the bottom layer of the hierachy, and each
-- subsequent element is the next layer up.
data DeepHarmonium (fs :: [Type -> Type -> Type]) (xs :: [Type])

-- | A trivial 1-layer harmonium.
type OneHarmonium x = DeepHarmonium '[] '[x]

-- | A 2-layer harmonium.
type Harmonium f z x = DeepHarmonium '[f] [z,x]

type Mixture z k = Harmonium Tensor z (Categorical Int k)


--- Functions ---


-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromNullMixture :: c # Mixture z 0 -> c # z
{-# INLINE fromNullMixture #-}
fromNullMixture = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toNullMixture :: c # z -> c # Mixture z 0
{-# INLINE toNullMixture #-}
toNullMixture = breakPoint

-- | Converts a 'OneHarmonium' into a standard exponential family distribution.
fromOneHarmonium :: c # OneHarmonium x -> c # x
{-# INLINE fromOneHarmonium #-}
fromOneHarmonium = breakPoint

-- | Converts an exponential family distribution into a 'OneHarmonium'.
toOneHarmonium :: c # x -> c # OneHarmonium x
{-# INLINE toOneHarmonium #-}
toOneHarmonium = breakPoint

-- | Adds a layer defined by an affine function to the bottom of a deep harmonium.
joinBottomHarmonium
    :: Dual c #> c # Affine f z y -- ^ Affine function
    -> c # DeepHarmonium fs (y : xs) -- ^ Upper part of the deep harmonium
    -> c # DeepHarmonium (f : fs) (z : y : xs) -- ^ Combined deep harmonium
{-# INLINE joinBottomHarmonium #-}
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the top layer off of a harmonium.
splitBottomHarmonium
    :: (Manifold z, Manifold (f z y))
    => c # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep Harmonium
    -> (Dual c #> c # Affine f z y, c # DeepHarmonium fs (y : xs)) -- ^ Affine function and upper part
{-# INLINE splitBottomHarmonium #-}
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fs z xs c . Manifold z
    => c # z -- ^ Bias step
    -> c # DeepHarmonium fs (z : xs) -- ^ Deep Harmonium
    -> c # DeepHarmonium fs (z : xs) -- ^ Biased deep harmonium
{-# INLINE biasBottom #-}
biasBottom pm0 dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
        (pmcs,css') = VS.splitAt dz . S.fromSized $ coordinates dhrm
        pmcs' = H.add pmcs $ S.fromSized (coordinates pm0)
     in Point . I.Vector $ pmcs' VS.++ css'

-- | Get the bias of the bottom layer of the given 'DeepHarmonium'.
getBottomBias
    :: forall fs z xs c . Manifold z
    => c # DeepHarmonium fs (z : xs) -- ^ Deep harmonium
    -> c # z -- ^ Bottom layer bias
{-# INLINE getBottomBias #-}
getBottomBias dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
     in Point . I.Vector . VS.take dz . S.fromSized $ coordinates dhrm


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs (fs :: [Type -> Type -> Type]) (xs :: [Type]) ss where

    -- | Given a 'DeepHarmonium' and an element of its sample space, partially
    -- updates the sample by resampling from the bottom to the top layer, but
    -- without updating the bottom layer itself.
    upwardPass
        :: Natural # DeepHarmonium fs xs -- ^ Deep harmonium
        -> [HList ss] -- ^ Initial sample
        -> Random r [HList ss] -- ^ Partial Gibbs resample

    -- | Generates an element of the sample space of a deep harmonium based by
    -- starting from a sample point from the bottom layer, and doing a naive
    -- upward sampling. This does not generate a true sample from the deep
    -- harmonium.
    initialPass
        :: Natural # DeepHarmonium fs xs -- Deep harmonium
        -> [Head ss] -- ^ Bottom layer sample
        -> Random r [HList ss] -- ^ Initial deep harmonium sample

-- | Harmonium transposition. Each defining layers are reversed, and the
-- defining bilinear functions are transposed.
class Manifold (DeepHarmonium fs xs) => TransposeHarmonium fs xs where
    transposeHarmonium :: Primal c => c # DeepHarmonium fs xs -> c # DeepHarmonium (Reverse fs) (Reverse xs)


-- | A single pass of Gibbs sampling. Infinite iteration of this function yields
-- a sample from the given 'DeepHarmonium'.
gibbsPass :: ( Manifold (DeepHarmonium fs (y : xs)), Gibbs (f : fs) (z : y : xs) ss
             , Map Mean Natural f z y, Generative Natural z, ExponentialFamily y )
  => Natural # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep Hamonium
  -> [HList ss] -- ^ Initial Sample
  -> Random r [HList ss] -- ^ Gibbs resample
{-# INLINE gibbsPass #-}
gibbsPass dhrm zyxs = do
    let yxs = snd $ hUnzip zyxs
        ys = fst $ hUnzip yxs
        f = fst $ splitBottomHarmonium dhrm
    zs <- mapM samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip zs yxs


--- Conjugation ---


-- | A conjugated distribution has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleConjugated fs xs ss where
    -- | A true sample from a conjugated harmonium.
    sampleConjugatedHarmonium
        :: Int -- ^ Sample Size
        -> Natural # Sum (Tail xs) -- ^ Conjugation parameters
        -> Natural # DeepHarmonium fs xs -- ^ Deep harmonium
        -> Random r [HList ss] -- ^ Deep harmonium sample

-- | Marginalize the bottom layer out of a deep harmonium.
marginalizeConjugatedHarmonium
    :: ( Manifold (DeepHarmonium fs (y : xs)), Map Mean Natural f z y, Manifold (Sum xs) )
      => Natural # Sum (y : xs) -- ^ Conjugation Parameters
      -> Natural # DeepHarmonium (f : fs) (z : y : xs) -- ^ Deep harmonium
      -> (Natural # Sum xs, Natural # DeepHarmonium fs (y : xs)) -- ^ Marginalized deep harmonium
{-# INLINE marginalizeConjugatedHarmonium #-}
marginalizeConjugatedHarmonium rprms dhrm =
    let dhrm' = snd $ splitBottomHarmonium dhrm
        (rprm,rprms') = splitSum rprms
     in (rprms', biasBottom rprm dhrm')

-- Mixture Models --

-- | The stochastic cross entropy differential of a mixture model.
mixtureRelativeEntropyUpperBound
    :: forall z n . ( ClosedFormExponentialFamily z, KnownNat n )
      => Natural # Mixture z n -- ^ Categorical harmonium
      -> Natural # Mixture z n -- ^ Categorical harmonium
      -> Double -- ^ Upper bound
{-# INLINE mixtureRelativeEntropyUpperBound #-}
mixtureRelativeEntropyUpperBound phrm qhrm =
    let pzc = fst $ splitBottomHarmonium phrm
        npc = snd $ splitMixture phrm
        spc = toSource npc
        qzc = fst $ splitBottomHarmonium qhrm
        qc = snd $ splitMixture qhrm
        wghts = (1 - S.sum (coordinates spc)) : listCoordinates spc
        smps = sampleSpace $ Proxy @ (Categorical Int n)
        pzs = pzc >$>* smps
        qzs = qzc >$>* smps
        dvg0 = weightedAverage (zip wghts . zipWith divergence pzs $ dualTransition <$> qzs)
     in divergence npc (transition qc) + dvg0


-- | A convenience function for building a categorical harmonium/mixture model.
joinMixture
    :: forall k z . ( KnownNat k, Legendre Natural z )
    => S.Vector (k+1) (Natural # z) -- ^ Mixture components
    -> Natural # Categorical Int k -- ^ Weights
    -> Natural # Mixture z k -- ^ Mixture Model
{-# INLINE joinMixture #-}
joinMixture nzs0 nx0 =
    let nz0 :: S.Vector 1 (Natural # z)
        (nz0,nzs0') = S.splitAt nzs0
        nz = S.head nz0
        nzs = S.map (<-> nz) nzs0'
        nzx = fromMatrix . S.fromColumns $ S.map coordinates nzs
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx = toOneHarmonium $ nx0 <-> rprms
     in joinBottomHarmonium affzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitMixture
    :: forall k z . ( KnownNat k, Legendre Natural z )
    => Natural # Mixture z k -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # z), Natural # Categorical Int k) -- ^ (components, weights)
{-# INLINE splitMixture #-}
splitMixture hrm =
    let (affzx,nx) = splitBottomHarmonium hrm
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx0 = fromOneHarmonium nx <+> rprms
        (nz,nzx) = splitAffine affzx
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (<+> nz) nzs
     in (S.cons nz nzs0',nx0)

-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleMixture
    :: ( Enum e, KnownNat n, Legendre Natural o
       , Generative Natural o, Manifold (Harmonium Tensor o (Categorical e n) ) )
       => Int
       -> Natural # Harmonium Tensor o (Categorical e n) -- ^ Categorical harmonium
       -> Random r (Sample (Harmonium Tensor o (Categorical e n))) -- ^ Sample
{-# INLINE sampleMixture #-}
sampleMixture k hrm = do
    let rx = snd . mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
    sampleConjugatedHarmonium k (toSingletonSum rx) hrm


--- Internal Functions ---


harmoniumBaseMeasure
    :: ExponentialFamily x
    => Proxy x
    -> Proxy (OneHarmonium x)
    -> SamplePoint (OneHarmonium x)
    -> Double
{-# INLINE harmoniumBaseMeasure #-}
harmoniumBaseMeasure prxyl _ (x :+: Null) =
     baseMeasure prxyl x

deepHarmoniumBaseMeasure
    :: (ExponentialFamily x, ExponentialFamily (DeepHarmonium fs xs))
    => Proxy x
    -> Proxy (DeepHarmonium fs xs)
    -> Proxy (DeepHarmonium (f : fs) (x : xs))
    -> SamplePoint (DeepHarmonium (f : fs) (x : xs))
    -> Double
{-# INLINE deepHarmoniumBaseMeasure #-}
deepHarmoniumBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     baseMeasure prxym xm * baseMeasure prxydhrm xs

mixtureExpectations
    :: ( KnownNat k, Legendre Natural z, ExponentialFamily z )
    => Natural # Mixture z k
    -> Mean # Mixture z k
{-# INLINE mixtureExpectations #-}
mixtureExpectations hrm =
    let (nzs,nx) = splitMixture hrm
        mx = toMean nx
        mzs = S.map dualTransition nzs
        wghts = categoricalWeights mx
        wmzs = S.zipWith (.>) wghts mzs
        mz = S.foldr1 (<+>) wmzs
        twmzs = S.tail wmzs
        mzx = transpose . fromRows $ twmzs
     in joinBottomHarmonium (joinAffine mz mzx) $ toOneHarmonium mx

mixtureParameters
    :: ( KnownNat k, Legendre Mean z, Legendre Natural z, ExponentialFamily z )
    => Mean # Mixture z k
    -> Natural # Mixture z k
{-# INLINE mixtureParameters #-}
mixtureParameters hrm =
    let (maff,mx0) = splitBottomHarmonium hrm
        (mz,mzx) = splitAffine maff
        mx = fromOneHarmonium mx0
        twmzs = toRows $ transpose mzx
        wmzs = S.cons (mz <-> S.foldr (<+>) zero twmzs) twmzs
        wghts = categoricalWeights mx
        mzs = S.zipWith (/>) wghts wmzs
     in joinMixture (S.map dualTransition mzs) (dualTransition mx)

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
conjugatedHarmoniumDensity
    :: forall f z x . ( Bilinear f z x, ExponentialFamily (Harmonium f z x), Map Mean Natural f z x
       , Legendre Natural z, Legendre Natural x, ExponentialFamily z, ExponentialFamily x )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f z x
      -> SamplePoint z
      -> Double
{-# INLINE conjugatedHarmoniumDensity #-}
conjugatedHarmoniumDensity (rho0,rprms) hrm ox =
    let (f,nl0) = splitBottomHarmonium hrm
        (no,nlo) = splitAffine f
        nl = fromOneHarmonium nl0
     in (* baseMeasure (Proxy @ z) ox) . exp $ sum
            [ sufficientStatistic ox <.> no
            , potential (nl <+> ox *<.< nlo)
            , negate $ potential (nl <+> rprms) + rho0 ]

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedHarmoniumDensity
    :: forall f z x . ( Bilinear f z x, ExponentialFamily (Harmonium f z x), Map Mean Natural f z x
       , Legendre Natural x, Legendre Natural z, ExponentialFamily x, ExponentialFamily z )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f z x
      -> SamplePoint z
      -> Double
{-# INLINE logConjugatedHarmoniumDensity #-}
logConjugatedHarmoniumDensity (rho0,rprms) hrm z =
    let (f,nx0) = splitBottomHarmonium hrm
        (nz,nzx) = splitAffine f
        nx = fromOneHarmonium nx0
     in log (baseMeasure (Proxy @ z) z) + sum
            [ sufficientStatistic z <.> nz
            , potential (nx <+> z *<.< nzx)
            , negate $ potential (nx <+> rprms) + rho0 ]


-- Why the heck does this produce different answers than the conjugated version in particular contexts!?!
---- | Computes the negative log-likelihood of a sample point of a mixture model.
--mixtureDensity
--    :: ( KnownNat n, AbsolutelyContinuous Natural z, Legendre Natural z )
--    => Natural # Mixture z n -- ^ Categorical Harmonium
--    -> SamplePoint z -- ^ Observation
--    -> Double -- ^ Negative log likelihood
--{-# INLINE mixtureDensity #-}
--mixtureDensity hrm x =
--    let (ncmpnts,nwghts) = splitMixture hrm
--        dnss = (`density` x) <$> S.toList ncmpnts
--        wghts = S.toList $ categoricalWeights nwghts
--     in weightedAverage $ zip wghts dnss

-- | Computes the negative log-likelihood of a sample point of a mixture model.
mixtureDensity
    :: ( Enum e, KnownNat k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE mixtureDensity #-}
mixtureDensity hrm =
    let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
     in conjugatedHarmoniumDensity rh0rx hrm

-- | Computes the negative log-likelihood of a sample point of a mixture model.
logMixtureDensity
    :: ( Enum e, KnownNat k, Legendre Natural o, ExponentialFamily o )
    => Natural # Harmonium Tensor o (Categorical e k) -- ^ Categorical Harmonium
    -> SamplePoint o -- ^ Observation
    -> Double -- ^ Negative log likelihood
{-# INLINE logMixtureDensity #-}
logMixtureDensity hrm =
    let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
     in logConjugatedHarmoniumDensity rh0rx hrm

-- Misc --

unnormalizedHarmoniumObservableDensity
    :: forall f z x . (ExponentialFamily z, Legendre Natural x, Bilinear f z x)
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableDensity hrm z =
    let (affzx,nx0) = splitBottomHarmonium hrm
        (nz,nzx) = splitAffine affzx
        nx = fromOneHarmonium nx0
        mz = sufficientStatistic z
     in baseMeasure (Proxy @ z) z * exp (nz <.> mz + potential (nx <+> mz <.< nzx))

-- | This might be inefficient due to the use of average point instead of a
-- slightly more complicated foldr.
harmoniumEmpiricalExpectations
    :: ( ExponentialFamily z, Map Mean Natural f x z
       , Bilinear f z x, Map Mean Natural f z x, Legendre Natural x)
    => Sample z -- ^ Model Samples
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Mean # Harmonium f z x -- ^ Harmonium expected sufficient statistics
{-# INLINE harmoniumEmpiricalExpectations #-}
harmoniumEmpiricalExpectations zs =
    harmoniumEmpiricalExpectations0 (sufficientStatistic <$> zs)

-- | This might be inefficient due to the use of average point instead of a
-- slightly more complicated foldr.
harmoniumEmpiricalExpectations0
    :: ( ExponentialFamily z, Map Mean Natural f x z
       , Bilinear f z x, Map Mean Natural f z x, Legendre Natural x)
    => [Mean # z] -- ^ Model Samples
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Mean # Harmonium f z x -- ^ Harmonium expected sufficient statistics
{-# INLINE harmoniumEmpiricalExpectations0 #-}
harmoniumEmpiricalExpectations0 mzs hrm =
    let aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        mxs = dualTransition <$> aff >$> mzs
        mzx = averagePoint $ zipWith (>.<) mzs mxs
        maff = joinAffine (averagePoint mzs) mzx
     in joinBottomHarmonium maff . toOneHarmonium $ averagePoint mxs


--- Instances ---


instance Manifold m => Manifold (DeepHarmonium fs '[m]) where
    type Dimension (DeepHarmonium fs '[m]) = Dimension m

instance (Manifold z, Manifold y, Manifold (f z y), Manifold (DeepHarmonium fs (y : xs)))
  => Manifold (DeepHarmonium (f : fs) (z : y : xs)) where
      type Dimension (DeepHarmonium (f : fs) (z : y : xs))
        = Dimension z + Dimension (f z y) + Dimension (DeepHarmonium fs (y : xs))

instance Generative c x => Generative c (OneHarmonium x) where
    {-# INLINE samplePoint #-}
    samplePoint = fmap (:+: Null) . samplePoint . fromOneHarmonium

instance ExponentialFamily x => ExponentialFamily (OneHarmonium x) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (x :+: Null) =
          toOneHarmonium $ sufficientStatistic x
      {-# INLINE baseMeasure #-}
      baseMeasure = harmoniumBaseMeasure Proxy

instance ( ExponentialFamily z, ExponentialFamily y, Bilinear f z y
         , ExponentialFamily (DeepHarmonium fs (y : xs)) )
  => ExponentialFamily (DeepHarmonium (f : fs) (z : y : xs)) where
      {-# INLINE sufficientStatistic #-}
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              pm = sufficientStatistic xm
              pn = sufficientStatistic xn
           in joinBottomHarmonium (joinAffine pm $ pm >.< pn) mdhrm
      {-# INLINE baseMeasure #-}
      baseMeasure = deepHarmoniumBaseMeasure Proxy Proxy

instance ( Bilinear f z x, ExponentialFamily z, Generative Natural x )
  => Gibbs '[f] '[z,x] where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zxs = initialPass dhrm . fst $ hUnzip zxs
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              xp = fromOneHarmonium dhrm'
          xs <- mapM (samplePoint . (<+>) xp) $ zs *<$< f
          return . hZip zs . hZip xs $ repeat Null

instance ( Bilinear f z y, Map Mean Natural g y x, Manifold (DeepHarmonium fs (x : xs))
         , ExponentialFamily z, ExponentialFamily x, Generative Natural y, Gibbs (g : fs) (y : x : xs) )
  => Gibbs (f : g : fs) (z : y : x : xs) where
      {-# INLINE upwardPass #-}
      upwardPass dhrm zyxs = do
          let (zs,yxs) = hUnzip zyxs
              (xs,xs') = hUnzip . snd $ hUnzip yxs
              (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              (g,_) = splitBottomHarmonium dhrm'
          ys <- mapM samplePoint $ zipWith (<+>) (g >$>* xs) (zs *<$< f)
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      {-# INLINE initialPass #-}
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              yp = fst . splitAffine . fst $ splitBottomHarmonium dhrm'
          ys <- mapM (samplePoint . (<+> yp)) $ zs *<$< f
          yxs' <- initialPass dhrm' ys
          return $ hZip zs yxs'

instance Manifold x => TransposeHarmonium '[] '[x] where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium = id

instance (Bilinear f z y, Bilinear f z y, TransposeHarmonium fs (y : xs))
  => TransposeHarmonium (f : fs) (z : y : xs) where
    {-# INLINE transposeHarmonium #-}
    transposeHarmonium dhrm =
        let (aff,dhrm') = splitBottomHarmonium dhrm
            (pm,pmtx) = splitAffine aff
            dhrm'' = transposeHarmonium dhrm'
         in Point . I.Vector . S.fromSized $ coordinates dhrm'' S.++ coordinates (transpose pmtx) S.++ coordinates pm

instance Generative Natural x => SampleConjugated '[] '[x] where
    {-# INLINE sampleConjugatedHarmonium #-}
    sampleConjugatedHarmonium k _ = sample k

instance ( Manifold (DeepHarmonium fs (y : xs)), Map Mean Natural f z y, Manifold (Sum xs)
         , ExponentialFamily y, SampleConjugated fs (y : xs), Generative Natural z )
  => SampleConjugated (f : fs) (z : y : xs) where
    {-# INLINE sampleConjugatedHarmonium #-}
    sampleConjugatedHarmonium k rprms dhrm = do
        let (pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleConjugatedHarmonium k rprms' $ biasBottom rprm dhrm'
        zs <- mapM samplePoint $ pf >$>* ys
        return . hZip zs $ hZip ys xs

instance ( KnownNat n, Legendre Natural z
       , Generative Natural z, Manifold (Mixture z n ) )
  => Generative Natural (Mixture z n) where
      {-# INLINE samplePoint #-}
      samplePoint hrm = head <$> sampleMixture 1 hrm

instance (KnownNat n, Legendre Natural z, ExponentialFamily z, Manifold (Mixture z n))
  => Legendre Natural (Mixture z n) where
      {-# INLINE potential #-}
      potential hrm =
          let (lkl,nx0) = splitBottomHarmonium hrm
              (rho0,rprms) = mixtureLikelihoodConjugationParameters lkl
              nx = fromOneHarmonium nx0
           in potential (nx <+> rprms) + rho0
      potentialDifferential = primalIsomorphism . mixtureExpectations

--- Graveyard

---- | The observable density of a categorical harmonium.
--mixtureDensity0
--    :: (KnownNat k, Enum e, Legendre Natural z, AbsolutelyContinuous Natural z)
--    => Natural # Harmonium Tensor z (Categorical e k) -- ^ Categorical harmonium
--    -> SamplePoint z -- ^ Observation
--    -> Double -- ^ Probablity density of the observation
--{-# INLINE mixtureDensity0 #-}
--mixtureDensity0 hrm x =
--    let (affzx,nx) = splitBottomHarmonium hrm
--        nz = fst $ splitAffine affzx
--        wghts = listCoordinates . toMean
--            $ snd (mixtureLikelihoodConjugationParameters affzx) <+> fromOneHarmonium nx
--        dxs0 = (`density` x) <$> affzx >$>* pointSampleSpace (fromOneHarmonium nx)
--        dx1 = density nz x * (1 - sum wghts)
--     in dx1 + sum (zipWith (*) wghts dxs0)
--
--

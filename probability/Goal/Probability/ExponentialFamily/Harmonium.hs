{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE TypeApplications,UndecidableInstances #-}
-- | An Exponential Family 'Harmonium' is a product exponential family with a
-- particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
-- A 'Mixture' model is a special case of harmonium, and a 'DeepHarmonium' is a
-- hierarchical extension of it.
module Goal.Probability.ExponentialFamily.Harmonium
    (
    -- * Harmoniums
      Harmonium
    -- ** Manipulation
    , splitHarmonium
    , joinHarmonium
    , transposeHarmonium
    -- ** Evaluation
    , unnormalizedHarmoniumObservableDensity
    , harmoniumExpectationStep
    -- * Mixture Models
    , Mixture
    , joinNaturalMixture
    , splitNaturalMixture
    , joinMeanMixture
    , splitMeanMixture
    , mixtureDensity
    , logMixtureDensity
    -- * Deep Harmoniums
    , OneHarmonium
    , DeepHarmonium
    -- ** Conversion
    , fromOneHarmonium
    , toOneHarmonium
    -- ** Construction
    , biasBottom
    , getBottomBias
    , splitBottomHarmonium
    , joinBottomHarmonium
    -- ** Sampling
    , initialPass
    , gibbsPass
    -- * Conjugated Harmoniums
    , marginalizeConjugatedHarmonium
    , sampleConjugatedHarmonium
    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Internal as I
import qualified Data.Vector.Storable as VS
import qualified Numeric.LinearAlgebra as H


--- Types ---


-- | A hierarchical generative model defined by exponential families. Note that
-- the first element @z@ is the bottom variable of the hierachy, and each
-- subsequent element @(f,x)@ represents the 'Bilinear' relation and subsequent
-- variable.
data DeepHarmonium z (fxs :: [(Type -> Type -> Type,Type)])

-- | A trivial 1-layer harmonium.
type OneHarmonium z = DeepHarmonium z '[]

-- | A 2-layer harmonium.
type Harmonium z f x = DeepHarmonium z '[ '(f,x)]

-- | A 'Mixture' model is simply a 'Harmonium' where the latent variable is
-- 'Categorical'.
type Mixture z k = Harmonium z Tensor (Categorical k)

-- | A type family for gathering the conjugation parameters from the type of a
-- 'DeepHarmonium'.
type family ConjugationParameters (fxs :: [(Type -> Type -> Type,Type)]) where
    ConjugationParameters '[] = '[]
    ConjugationParameters ('(f,z) ': fxs) = z ': ConjugationParameters fxs



--- Functions ---


-- | Converts a 'OneHarmonium' into a standard 'ExponentialFamily' distribution.
fromOneHarmonium :: c # OneHarmonium x -> c # x
fromOneHarmonium = breakPoint

-- | Converts an 'ExponentialFamily' distribution into a 'OneHarmonium'.
toOneHarmonium :: c # x -> c # OneHarmonium x
toOneHarmonium = breakPoint

-- | Adds a variable to the bottom of a 'DeepHarmonium' given an 'Affine'
-- function.
joinBottomHarmonium
    :: c #> Affine f z y -- ^ Affine function
    -> c # DeepHarmonium y gxs -- ^ Upper part of the deep harmonium
    -> c # DeepHarmonium z ('(f,y) : gxs) -- ^ Combined deep harmonium
joinBottomHarmonium pf dhrm =
    Point $ coordinates pf S.++ coordinates dhrm

-- | Splits the bottom 'Affine' function off of a 'DeepHarmonium'.
splitBottomHarmonium
    :: (Manifold z, Manifold (f z y))
    => c # DeepHarmonium z ('(f,y) : gxs) -- ^ Deep Harmonium
    -> (c #> Affine f z y, c # DeepHarmonium y gxs) -- ^ Affine function and upper part
splitBottomHarmonium dhrm =
    let (affcs,dcs) = S.splitAt $ coordinates dhrm
     in (Point affcs, Point dcs)

-- | Creates a 'Harmonium' from component parameters.
joinHarmonium
    :: (Manifold z, Manifold x, Manifold (f z x))
    => c # z -- ^ Visible layer biases
    -> c #> f z x -- ^ ^ Interaction parameters
    -> c # x -- ^ Hidden layer Biases
    -> c # Harmonium z f x -- ^ Harmonium
joinHarmonium nz nzx nx =
    let nx0 = toOneHarmonium nx
        aff = joinAffine nz nzx
     in joinBottomHarmonium aff nx0

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium
    :: (Manifold z, Manifold x, Manifold (f z x))
    => c # Harmonium z f x -- ^ Harmonium
    -> (c # z, c #> f z x, c # x) -- ^ Biases and interaction parameters
splitHarmonium hrm =
    let (f,nx0) = splitBottomHarmonium hrm
        (nz,nzx) = splitAffine f
        nx = fromOneHarmonium nx0
     in (nz,nzx,nx)

-- | Translate the bias of the bottom layer by the given 'Point'.
biasBottom
    :: forall fxs z c . Manifold z
    => c # z -- ^ Bias step
    -> c # DeepHarmonium z fxs -- ^ Deep Harmonium
    -> c # DeepHarmonium z fxs -- ^ Biased deep harmonium
biasBottom pm0 dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
        (pmcs,css') = VS.splitAt dz . S.fromSized $ coordinates dhrm
        pmcs' = H.add pmcs $ S.fromSized (coordinates pm0)
     in Point . I.Vector $ pmcs' VS.++ css'

-- | Get the bias of the bottom layer of the given 'DeepHarmonium'.
getBottomBias
    :: forall fxs z c . Manifold z
    => c # DeepHarmonium z fxs -- ^ Deep harmonium
    -> c # z -- ^ Bottom layer bias
getBottomBias dhrm =
    let dz = natValInt (Proxy @ (Dimension z))
     in Point . I.Vector . VS.take dz . S.fromSized $ coordinates dhrm


--- Classes ---


-- | 'Gibbs' deep harmoniums can be sampled through Gibbs sampling.
class Gibbs z (fxs :: [(Type -> Type -> Type,Type)]) where

    -- | Given a 'DeepHarmonium' and an element of its sample space, partially
    -- resamples the sample by resampling from the bottom to the top layer, but
    -- without updating the bottom layer itself.
    upwardPass
        :: Natural # DeepHarmonium z fxs -- ^ Deep harmonium
        -> Sample (DeepHarmonium z fxs) -- ^ Initial sample
        -> Random r (Sample (DeepHarmonium z fxs)) -- ^ Partial Gibbs resample

    -- | Generates an element of the sample space of a deep harmonium by
    -- starting from a sample point from the bottom layer, and doing naive
    -- upward sampling. This does not generate a true sample from the deep
    -- harmonium.
    initialPass
        :: Natural # DeepHarmonium z fxs -- Deep harmonium
        -> Sample z -- ^ Bottom layer sample
        -> Random r (Sample (DeepHarmonium z fxs)) -- ^ Initial deep harmonium sample

-- | A single pass of Gibbs sampling. Infinite iteration of this function yields
-- a sample from the given 'DeepHarmonium'.
gibbsPass :: ( Manifold (DeepHarmonium y gxs), Gibbs z ('(f,y) : gxs), Map Mean Natural f z y
             , Bilinear f z y, Generative Natural z, ExponentialFamily y )
  => Natural # DeepHarmonium z ('(f,y) : gxs) -- ^ Deep Hamonium
  -> Sample (DeepHarmonium z ('(f,y) : gxs)) -- ^ Initial Sample
  -> Random r (Sample (DeepHarmonium z ('(f,y) : gxs))) -- ^ Gibbs resample
gibbsPass dhrm zyxs = do
    let yxs = snd $ hUnzip zyxs
        ys = fst $ hUnzip yxs
        f = fst $ splitBottomHarmonium dhrm
    zs <- mapM samplePoint $ f >$>* ys
    upwardPass dhrm $ hZip zs yxs


--- Conjugation ---


-- | A conjugated harmonium has a number of computational features, one of
-- which is being able to generate samples from the model with a single downward
-- pass.
class SampleConjugated z (fxs :: [(Type -> Type -> Type,Type)]) where
    -- | A true sample from a conjugated harmonium.
    sampleConjugatedHarmonium
        :: Int -- ^ Sample Size
        -> Natural # Sum (ConjugationParameters fxs) -- ^ Conjugation parameters
        -> Natural # DeepHarmonium z fxs -- ^ Deep harmonium
        -> Random r (Sample (DeepHarmonium z fxs))  -- ^ Deep harmonium sample

-- | Marginalize the bottom layer out of a deep harmonium.
marginalizeConjugatedHarmonium
    :: ( Manifold (DeepHarmonium y gxs), Map Mean Natural f z y
       , Manifold (Sum (ConjugationParameters gxs)) )
      => Natural # Sum (ConjugationParameters ('(f,y) : gxs)) -- ^ Conjugation Parameters
      -> Natural # DeepHarmonium z ('(f,y) : gxs) -- ^ Deep harmonium
      -> ( Natural # Sum (ConjugationParameters gxs)
         , Natural # DeepHarmonium y gxs ) -- ^ Marginalized deep harmonium
marginalizeConjugatedHarmonium rprms dhrm =
    let dhrm' = snd $ splitBottomHarmonium dhrm
        (rprm,rprms') = splitSum rprms
     in (rprms', biasBottom rprm dhrm')

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedHarmoniumDensities
    :: forall f z x .  ( Bilinear f z x, Map Mean Natural f x z
                       , LegendreExponentialFamily z, LegendreExponentialFamily x )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium z f x
      -> Sample z
      -> [Double]
logConjugatedHarmoniumDensities (rho0,rprms) hrm zs =
    let (nz,nzx,nx) = splitHarmonium hrm
        cnds = potential . (nx +) <$> (zs *<$< nzx)
        dts = dotMap nz $ sufficientStatistic <$> zs
        scls = [ logBaseMeasure (Proxy @ z) z - potential (nx + rprms) - rho0 | z <- zs ]
        zipper a b c = a + b + c
     in zipWith3 zipper cnds dts scls


-- Misc --

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableDensity
    :: forall f z x . ( ExponentialFamily z, LegendreExponentialFamily x
                      , Map Mean Natural f x z, Bilinear f z x)
    => Natural # Harmonium z f x
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableDensity hrm z =
    let (affzx,nx0) = splitBottomHarmonium hrm
        (nz,nzx) = splitAffine affzx
        nx = fromOneHarmonium nx0
        mz = sufficientStatistic z
     in exp (nz <.> mz + potential (nx + mz <.< nzx) + logBaseMeasure (Proxy @ z) z)


---- Mixture Models --


-- | A convenience function for building a categorical harmonium/mixture model.
joinNaturalMixture
    :: forall k z . ( KnownNat k, LegendreExponentialFamily z )
    => S.Vector (k+1) (Natural # z) -- ^ Mixture components
    -> Natural # Categorical k -- ^ Weights
    -> Natural # Mixture z k -- ^ Mixture Model
joinNaturalMixture nzs0 nx0 =
    let nz0 :: S.Vector 1 (Natural # z)
        (nz0,nzs0') = S.splitAt nzs0
        nz = S.head nz0
        nzs = S.map (subtract nz) nzs0'
        nzx = fromMatrix . S.fromColumns $ S.map coordinates nzs
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx = toOneHarmonium $ nx0 - rprms
     in joinBottomHarmonium affzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitNaturalMixture
    :: forall k z . ( KnownNat k, LegendreExponentialFamily z )
    => Natural # Mixture z k -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # z), Natural # Categorical k) -- ^ (components, weights)
splitNaturalMixture hrm =
    let (affzx,nx) = splitBottomHarmonium hrm
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx0 = fromOneHarmonium nx + rprms
        (nz,nzx) = splitAffine affzx
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (+ nz) nzs
     in (S.cons nz nzs0',nx0)

-- | The density over the observable variables of a mixture model.
mixtureDensity
    :: (ExponentialFamily z, LegendreExponentialFamily z, KnownNat k)
    => (Natural # Mixture z k)
    -> SamplePoint z
    -> Double
mixtureDensity mxt z =
    let rho0rprms = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium mxt
     in exp . head $ logConjugatedHarmoniumDensities rho0rprms mxt [z]

-- | The log-density over the observable variables of a mixture model.
logMixtureDensity
    :: (ExponentialFamily z, LegendreExponentialFamily z, KnownNat k)
    => (Natural # Mixture z k)
    -> SamplePoint z
    -> Double
logMixtureDensity mxt z =
    let rho0rprms = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium mxt
     in head $ logConjugatedHarmoniumDensities rho0rprms mxt [z]

-- | Swap the biases and 'transpose' the interaction parameters of the given 'Harmonium'.
transposeHarmonium
    :: Bilinear f z x
    => Natural # Harmonium z f x
    -> Natural # Harmonium x f z
transposeHarmonium dhrm =
        let (affzx,nx') = splitBottomHarmonium dhrm
            nx = fromOneHarmonium nx'
            (nz,nzx) = splitAffine affzx
            affxz = joinAffine nx $ transpose nzx
         in joinBottomHarmonium affxz $ toOneHarmonium nz


-- | Computes the joint expectations of a harmonium based on a sample from the
-- observable layer.
harmoniumExpectationStep
    :: ( ExponentialFamily z, Map Mean Natural f x z
       , Bilinear f z x, LegendreExponentialFamily x )
    => Sample z -- ^ Model Samples
    -> Natural # Harmonium z f x -- ^ Harmonium
    -> Mean # Harmonium z f x -- ^ Harmonium expected sufficient statistics
harmoniumExpectationStep zs hrm =
    let mzs = sufficientStatistic <$> zs
        aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        mxs = transition <$> aff >$> mzs
        mzx = (>$<) mzs mxs
     in joinHarmonium (average mzs) mzx $ average mxs


--- Internal ---


-- | Generates a sample from a categorical harmonium, a.k.a a mixture distribution.
sampleMixture
    :: ( KnownNat n, LegendreExponentialFamily z, Generative Natural z, Manifold (Mixture z n) )
       => Int
       -> Natural # Mixture z n -- ^ Categorical harmonium
       -> Random r (Sample (Mixture z n)) -- ^ Sample
sampleMixture k hrm = do
    let rx = snd . mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
    sampleConjugatedHarmonium k (toSingletonSum rx) hrm

harmoniumLogBaseMeasure
    :: ExponentialFamily x
    => Proxy x
    -> Proxy (OneHarmonium x)
    -> HList '[SamplePoint x]
    -> Double
harmoniumLogBaseMeasure prxyl _ (x :+: Null) =
     logBaseMeasure prxyl x

deepHarmoniumLogBaseMeasure
    :: (ExponentialFamily z, ExponentialFamily (DeepHarmonium y gxs))
    => Proxy z
    -> Proxy (DeepHarmonium y gxs)
    -> Proxy (DeepHarmonium z ('(f,y) : gxs))
    -> SamplePoint (DeepHarmonium z ('(f,y) : gxs))
    -> Double
deepHarmoniumLogBaseMeasure prxym prxydhrm _ (xm :+: xs) =
     logBaseMeasure prxym xm + logBaseMeasure prxydhrm xs

joinMeanMixture
    :: (KnownNat k, Manifold z)
    => S.Vector (k+1) (Mean # z) -- ^ Mixture components
    -> Mean # Categorical k -- ^ Weights
    -> Mean # Mixture z k
joinMeanMixture mzs mx =
    let wghts = categoricalWeights mx
        wmzs = S.zipWith (.>) wghts mzs
        mz = S.foldr1 (+) wmzs
        twmzs = S.tail wmzs
        mzx = transpose . fromRows $ twmzs
     in joinBottomHarmonium (joinAffine mz mzx) $ toOneHarmonium mx

splitMeanMixture
    :: ( KnownNat k, DuallyFlatExponentialFamily z )
    => Mean # Mixture z k
    -> (S.Vector (k+1) (Mean # z), Mean # Categorical k)
splitMeanMixture hrm =
    let (maff,mx0) = splitBottomHarmonium hrm
        (mz,mzx) = splitAffine maff
        mx = fromOneHarmonium mx0
        twmzs = toRows $ transpose mzx
        wmzs = S.cons (mz - S.foldr (+) 0 twmzs) twmzs
        wghts = categoricalWeights mx
        mzs = S.zipWith (/>) wghts wmzs
     in (mzs,mx)

-- | Computes the conjugation parameters of a likelihood defined by a categorical latent variable.
mixtureLikelihoodConjugationParameters
    :: (KnownNat k, LegendreExponentialFamily z)
    => Natural #> (z <* Categorical k) -- ^ Categorical likelihood
    -> (Double, Natural # Categorical k) -- ^ Conjugation parameters
mixtureLikelihoodConjugationParameters aff =
    let (nz,nzx) = splitAffine aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz + Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)


--- Instances ---


instance Manifold x => Manifold (DeepHarmonium x '[]) where
    type Dimension (DeepHarmonium x '[]) = Dimension x

instance (Manifold z, Manifold y, Manifold (f z y), Manifold (DeepHarmonium y gxs))
  => Manifold (DeepHarmonium z ('(f,y) : gxs)) where
      type Dimension (DeepHarmonium z ('(f,y) : gxs))
        = Dimension z + Dimension (f z y) + Dimension (DeepHarmonium y gxs)

instance Manifold (DeepHarmonium z fxs) => Statistical (DeepHarmonium z fxs) where
    type SamplePoint (DeepHarmonium z fxs) = HList (SamplePoints (z : ConjugationParameters fxs))

instance Generative c x => Generative c (OneHarmonium x) where
    samplePoint = fmap (:+: Null) . samplePoint . fromOneHarmonium

instance ExponentialFamily x => ExponentialFamily (OneHarmonium x) where
      sufficientStatistic (x :+: Null) =
          toOneHarmonium $ sufficientStatistic x
      averageSufficientStatistic hxs =
          let xs = hHead <$> hxs
           in toOneHarmonium $ averageSufficientStatistic xs
      logBaseMeasure = harmoniumLogBaseMeasure Proxy

instance ( ExponentialFamily z, ExponentialFamily y, Bilinear f z y
         , ExponentialFamily (DeepHarmonium y gxs) )
  => ExponentialFamily (DeepHarmonium z ('(f,y) : gxs)) where
      sufficientStatistic (xm :+: xn :+: xs) =
          let mdhrm = sufficientStatistic $ xn :+: xs
              pm = sufficientStatistic xm
              pn = sufficientStatistic xn
           in joinBottomHarmonium (joinAffine pm $ pm >.< pn) mdhrm
      averageSufficientStatistic xmxnxs =
          let (xms,xns,xss) = unzip3 [ (xm,xn,xs) | (xm :+: xn :+: xs) <- xmxnxs ]
              mdhrm = averageSufficientStatistic $ zipWith (:+:) xns xss
              pms = sufficientStatistic <$> xms
              pns = sufficientStatistic <$> xns
           in joinBottomHarmonium (joinAffine (average pms) $ pms >$< pns) mdhrm
      logBaseMeasure = deepHarmoniumLogBaseMeasure Proxy Proxy

instance ( Map Mean Natural f x z, Bilinear f z x, ExponentialFamily z, Generative Natural x )
  => Gibbs z '[ '(f,x)] where
      upwardPass dhrm zxs = initialPass dhrm . fst $ hUnzip zxs
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              xp = fromOneHarmonium dhrm'
          xs <- mapM (samplePoint . (+) xp) $ zs *<$< f
          return . hZip zs . hZip xs $ repeat Null

instance ( Map Mean Natural f y z, Bilinear f z y, Map Mean Natural g y x, Manifold (DeepHarmonium y fxs)
         , ExponentialFamily z, Generative Natural y, ExponentialFamily x, Gibbs y ('(g,x) : fxs) )
  => Gibbs z ('(f,y) : '(g,x) : fxs) where
      upwardPass dhrm zyxs = do
          let (zs,yxs) = hUnzip zyxs
              (xs,xs') = hUnzip . snd $ hUnzip yxs
              (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              (g,_) = splitBottomHarmonium dhrm'
          ys <- mapM samplePoint $ zipWith (+) (g >$>* xs) (zs *<$< f)
          yxs' <- upwardPass dhrm' . hZip ys $ hZip xs xs'
          return $ hZip zs yxs'
      initialPass dhrm zs = do
          let (aff,dhrm') = splitBottomHarmonium dhrm
              f = snd $ splitAffine aff
              yp = fst . splitAffine . fst $ splitBottomHarmonium dhrm'
          ys <- mapM (samplePoint . (+ yp)) $ zs *<$< f
          yxs' <- initialPass dhrm' ys
          return $ hZip zs yxs'

instance Generative Natural z => SampleConjugated z '[] where
    sampleConjugatedHarmonium k _ = sample k

instance ( Manifold (DeepHarmonium y gxs), Map Mean Natural f z y
         , Manifold (Sum (ConjugationParameters gxs)), ExponentialFamily y
         , SampleConjugated y gxs, Generative Natural z )
  => SampleConjugated z ('(f,y) : gxs) where
    sampleConjugatedHarmonium k rprms dhrm = do
        let (pf,dhrm') = splitBottomHarmonium dhrm
            (rprm,rprms') = splitSum rprms
        (ys,xs) <- fmap hUnzip . sampleConjugatedHarmonium k rprms' $ biasBottom rprm dhrm'
        zs <- mapM samplePoint $ pf >$>* ys
        return . hZip zs $ hZip ys xs

instance ( KnownNat n, LegendreExponentialFamily z, Generative Natural z, Manifold (Mixture z n ) )
         => Generative Natural (Mixture z n) where
    sample = sampleMixture

instance (KnownNat k, LegendreExponentialFamily z) => Transition Natural Mean (Mixture z k) where
    transition nhrm =
        let (nzs,nx) = splitNaturalMixture nhrm
            mx = toMean nx
            mzs = S.map transition nzs
         in joinMeanMixture mzs mx

instance (KnownNat k, DuallyFlatExponentialFamily z) => Transition Mean Natural (Mixture z k) where
    transition mhrm =
        let (mzs,mx) = splitMeanMixture mhrm
            nx = transition mx
            nzs = S.map transition mzs
         in joinNaturalMixture nzs nx

instance (KnownNat n, LegendreExponentialFamily z) => Legendre (Mixture z n) where
      type PotentialCoordinates (Mixture z n) = Natural
      potential hrm =
          let (lkl,nx0) = splitBottomHarmonium hrm
              (rho0,rprms) = mixtureLikelihoodConjugationParameters lkl
              nx = fromOneHarmonium nx0
           in potential (nx + rprms) + rho0

instance (KnownNat n, DuallyFlatExponentialFamily z) => DuallyFlat (Mixture z n) where
    dualPotential mhrm =
        let nhrm = toNatural mhrm
         in mhrm <.> nhrm - potential nhrm

instance (KnownNat n, LegendreExponentialFamily z)
  => AbsolutelyContinuous Natural (Mixture z n) where
    densities = exponentialFamilyDensities

instance ( LegendreExponentialFamily z, KnownNat n, SamplePoint z ~ t )
  => LogLikelihood Natural (Mixture z n) t where
    logLikelihood xs hrm =
        let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
         in average $ logConjugatedHarmoniumDensities rh0rx hrm xs
    logLikelihoodDifferential zs hrm =
        let pxs = harmoniumExpectationStep zs hrm
            qxs = transition hrm
         in pxs - qxs

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE TypeApplications, UndecidableInstances #-}
-- | An Exponential Family 'Harmonium' is a product exponential family with a
-- particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
-- A 'Mixture' model is a special case of harmonium.
module Goal.Graphical.Models.Harmonium
    (
    -- * Harmoniums
      AffineHarmonium (AffineHarmonium)
    , Harmonium
    -- ** Constuction
    , splitHarmonium
    , joinHarmonium
    -- ** Manipulation
    , transposeHarmonium
    -- ** Evaluation
    , expectationStep
    -- ** Sampling
    , initialPass
    , gibbsPass
    -- ** Mixture Models
    , Mixture
    , AffineMixture
    , joinNaturalMixture
    , splitNaturalMixture
    , joinMeanMixture
    , splitMeanMixture
    , joinSourceMixture
    , splitSourceMixture
    -- ** Linear Gaussian Harmoniums
    , LinearGaussianHarmonium
    , FullGaussianHarmonium
    , DiagonalGaussianHarmonium
    , IsotropicGaussianHarmonium
    -- ** Conjugated Harmoniums
    , ConjugatedLikelihood (conjugationParameters)
    , harmoniumConjugationParameters
    , joinConjugatedHarmonium
    , splitConjugatedHarmonium
    -- * Hierarchical Mixtures of Gaussians
    , HierarchicalMixtureOfGaussians
    , FullHMoG
    , DiagonalHMoG
    , IsotropicHMoG
    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models

import qualified Goal.Core.Vector.Storable as S


--- Types ---


-- | A 2-layer harmonium.
newtype AffineHarmonium f x0 z0 x z = AffineHarmonium (Affine f x0 x z0, z)

type Harmonium f x z = AffineHarmonium f x z x z

type instance Observation (AffineHarmonium f x0 z0 x z) = SamplePoint x

-- | A 'Mixture' model is simply a 'AffineHarmonium' where the latent variable is
-- 'Categorical'.
type Mixture x k = Harmonium Tensor x (Categorical k)

-- | A 'Mixture' where only a subset of the component parameters are mixed.
type AffineMixture x0 x k =
    AffineHarmonium Tensor x0 (Categorical k) x (Categorical k)

-- | A `MultivariateNormal` reintrepreted as a join distribution over two component `MultivariateNormal`s.
type LinearGaussianHarmonium f n k =
    AffineHarmonium Tensor (MVNMean n) (MVNMean k) (MultivariateNormal f n) (FullNormal k)

-- | A `LinearGaussianHarmonium` with all covariances.
type FullGaussianHarmonium n k = LinearGaussianHarmonium MVNCovariance n k
-- | A `LinearGaussianHarmonium` with a diagonal covariance between oservable variables.
type IsotropicGaussianHarmonium n k = LinearGaussianHarmonium Scale n k
-- | A `LinearGaussianHarmonium` with a scale covariance between oservable variables.
type DiagonalGaussianHarmonium n k = LinearGaussianHarmonium Diagonal n k

type HierarchicalMixtureOfGaussians f n m k =
    AffineHarmonium Tensor (MVNMean n) (MVNMean m) (MultivariateNormal f n)
    (Mixture (FullNormal m) k)

type FullHMoG n m k = HierarchicalMixtureOfGaussians MVNCovariance n m k
type DiagonalHMoG n m k = HierarchicalMixtureOfGaussians Diagonal n m k
type IsotropicHMoG n m k = HierarchicalMixtureOfGaussians Scale n m k

--- Classes ---


-- | The conjugation parameters of a conjugated likelihood.
class ( ExponentialFamily x, ExponentialFamily z, Map Natural f x0 z0
      , Translation x x0 , Translation z z0 )
  => ConjugatedLikelihood f x0 z0 x z where
    conjugationParameters
        :: Natural # Affine f x0 x z0 -- ^ Categorical likelihood
        -> (Double, Natural # z) -- ^ Conjugation parameters


--- Functions ---


-- Construction --

-- | Creates a 'Harmonium' from component parameters.
joinHarmonium
    :: (Manifold x, Manifold z, Manifold (f x0 z0))
    => c # x -- ^ Visible layer biases
    -> c # f x0 z0 -- ^ ^ Interaction parameters
    -> c # z -- ^ Hidden layer Biases
    -> c # AffineHarmonium f x0 z0 x z -- ^ Harmonium
joinHarmonium nx nx0z0 = join (join nx nx0z0)

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium
    :: (Manifold x, Manifold z, Manifold (f x0 z0))
    => c # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> (c # x, c # f x0 z0, c # z) -- ^ Biases and interaction parameters
splitHarmonium hrm =
    let (fxz0,nz) = split hrm
        (nx,nx0z0) = split fxz0
     in (nx,nx0z0,nz)

-- | Build a mixture model in source coordinates.
joinSourceMixture
    :: (KnownNat k, Manifold x)
    => S.Vector (k+1) (Source # x) -- ^ Mixture components
    -> Source # Categorical k -- ^ Weights
    -> Source # Mixture x k
joinSourceMixture sxs sk =
    let (sx,sxs') = S.splitAt sxs
        aff = join (S.head sx) (fromColumns sxs')
     in join aff sk

-- | Build a mixture model in source coordinates.
splitSourceMixture
    :: (KnownNat k, Manifold x)
    => Source # Mixture x k
    -> (S.Vector (k+1) (Source # x), Source # Categorical k)
splitSourceMixture mxmdl =
    let (aff,sk) = split mxmdl
        (sx0,sxs0') = split aff
     in (S.cons sx0 $ toColumns sxs0' ,sk)

-- | Build a mixture model in mean coordinates.
joinMeanMixture
    :: (KnownNat k, Manifold x)
    => S.Vector (k+1) (Mean # x) -- ^ Mixture components
    -> Mean # Categorical k -- ^ Weights
    -> Mean # Mixture x k
joinMeanMixture mxs mk =
    let wghts = categoricalWeights mk
        wmxs = S.zipWith (.>) wghts mxs
        mx = S.foldr1 (+) wmxs
        twmxs = S.tail wmxs
        mxk = transpose . fromRows $ twmxs
     in joinHarmonium mx mxk mk

-- | Split a mixture model in mean coordinates.
splitMeanMixture
    :: ( KnownNat k, LegendreExponentialFamily x )
    => Mean # Mixture x k
    -> (S.Vector (k+1) (Mean # x), Mean # Categorical k)
splitMeanMixture hrm =
    let (mx,mxz,mk) = splitHarmonium hrm
        twmxs = toRows $ transpose mxz
        wmxs = S.cons (mx - S.foldr (+) 0 twmxs) twmxs
        wghts = categoricalWeights mk
        mxs = S.zipWith (/>) wghts wmxs
     in (mxs,mk)

-- | A convenience function for building a categorical harmonium/mixture model.
joinNaturalMixture
    :: forall k x . ( KnownNat k, LegendreExponentialFamily x )
    => S.Vector (k+1) (Natural # x) -- ^ Mixture components
    -> Natural # Categorical k -- ^ Weights
    -> Natural # Mixture x k -- ^ Mixture Model
joinNaturalMixture nxs0 nk0 =
    let nx0 :: S.Vector 1 (Natural # x)
        (nx0,nxs0') = S.splitAt nxs0
        nx = S.head nx0
        nxs = S.map (subtract nx) nxs0'
        nxk = fromColumns nxs
        affxk = join nx nxk
        rprms = snd $ conjugationParameters affxk
        nk = nk0 - rprms
     in joinHarmonium nx nxk nk

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitNaturalMixture
    :: forall k x . ( KnownNat k, LegendreExponentialFamily x )
    => Natural # Mixture x k -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # x), Natural # Categorical k) -- ^ (components, weights)
splitNaturalMixture hrm =
    let (nx,nxk,nk) = splitHarmonium hrm
        affxk = join nx nxk
        rprms = snd $ conjugationParameters affxk
        nk0 = nk + rprms
        nxs = toColumns nxk
        nxs0' = S.map (+ nx) nxs
     in (S.cons nx nxs0',nk0)


-- Manipulation --

-- | Swap the biases and 'transpose' the interaction parameters of the given 'Harmonium'.
transposeHarmonium
    :: (Bilinear c f x0 z0, Manifold x, Manifold z)
    => c # AffineHarmonium f x0 z0 x z
    -> c # AffineHarmonium f z0 x0 z x
transposeHarmonium hrm =
        let (nz,nyx,nw) = splitHarmonium hrm
         in joinHarmonium nw (transpose nyx) nz

-- Evaluation --

-- | Computes the joint expectations of a harmonium based on a sample from the
-- observable layer.
expectationStep
    :: ( ExponentialFamily x, LegendreExponentialFamily z
       , Translation x x0, Translation z z0
       , Bilinear Mean f x0 z0, Bilinear Natural f x0 z0 )
    => Sample x -- ^ Model Samples
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> Mean # AffineHarmonium f x0 z0 x z -- ^ Harmonium expected sufficient statistics
expectationStep xs hrm =
    let mxs = sufficientStatistic <$> xs
        mx0s = anchor <$> mxs
        pstr = fst . split $ transposeHarmonium hrm
        mzs = transition <$> pstr >$> mx0s
        mz0s = anchor <$> mzs
        mx0z0 = (>$<) mx0s mz0s
     in joinHarmonium (average mxs) mx0z0 $ average mzs

---- Sampling --

-- | Initialize a Gibbs chain from a set of observations.
initialPass
    :: forall f x0 z0 x z
    . ( Manifold z, ExponentialFamily x, Translation x x0, Translation z z0
      , Generative Natural z, Bilinear Natural f x0 z0)
    => Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> Sample x -- ^ Model Samples
    -> Random (Sample (x, z))
initialPass hrm xs = do
    let pstr = fst . split $ transposeHarmonium hrm
        mxs :: [Mean # x]
        mxs = sufficientStatistic <$> xs
        mx0s = anchor <$> mxs
    zs <- mapM samplePoint $ pstr >$> mx0s
    return $ zip xs zs

-- | Update a 'Sample' with Gibbs sampling.
gibbsPass
    :: forall f x0 z0 x z
    . ( ExponentialFamily z, ExponentialFamily x, Translation x x0, Translation z z0
      , Generative Natural z,Generative Natural x, Bilinear Natural f x0 z0 )
    => Natural # AffineHarmonium f x0 z0 x z -- ^ Harmonium
    -> Sample (x, z)
    -> Random (Sample (x, z))
gibbsPass hrm xzs = do
    let zs = snd <$> xzs
        mzs :: [Mean # z]
        mzs = sufficientStatistic <$> zs
        mz0s = anchor <$> mzs
        pstr = fst . split $ transposeHarmonium hrm
        lkl = fst $ split hrm
    xs' <- mapM samplePoint $ lkl >$> mz0s
    let mxs' :: [Mean # x]
        mxs' = sufficientStatistic <$> xs'
        mx0s' = anchor <$> mxs'
    zs' <- mapM samplePoint $ pstr >$> mx0s'
    return $ zip xs' zs'

-- Conjugation --

-- | The conjugation parameters of a conjugated `Harmonium`.
harmoniumConjugationParameters
    :: ConjugatedLikelihood f x0 z0 x z
    => Natural # AffineHarmonium f x0 z0 x z -- ^ Categorical likelihood
    -> (Double, Natural # z) -- ^ Conjugation parameters
harmoniumConjugationParameters hrm =
    conjugationParameters . fst $ split hrm

-- | The conjugation parameters of a conjugated `Harmonium`.
splitConjugatedHarmonium
    :: ConjugatedLikelihood f x0 z0 x z
    => Natural # AffineHarmonium f x0 z0 x z
    -> (Natural # Affine f x0 x z0, Natural # z)
splitConjugatedHarmonium hrm =
    let (lkl,nw) = split hrm
        cw = snd $ conjugationParameters lkl
     in (lkl,nw + cw)

-- | The conjugation parameters of a conjugated `Harmonium`.
joinConjugatedHarmonium
    :: ConjugatedLikelihood f x0 z0 x z
    => Natural # Affine f x0 x z0 -- ^ Conjugation parameters
    -> Natural # z
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Categorical likelihood
joinConjugatedHarmonium lkl nw =
    let cw = snd $ conjugationParameters lkl
     in join lkl $ nw - cw

-- | The conjugation parameters of a conjugated `Harmonium`.
sampleConjugated
    :: forall f x0 z0 x z
     . ( ConjugatedLikelihood f x0 z0 x z, Generative Natural x
       , Generative Natural z, Map Natural f x0 z0 )
    => Int
    -> Natural # AffineHarmonium f x0 z0 x z -- ^ Categorical likelihood
    -> Random (Sample (x,z)) -- ^ Conjugation parameters
sampleConjugated n hrm = do
    let (lkl,nz) = splitConjugatedHarmonium hrm
    zs <- sample n nz
    let mzs :: [Mean # z]
        mzs = sufficientStatistic <$> zs
    xs <- mapM samplePoint $ lkl >$+> mzs
    return $ zip xs zs

-- | The conjugation parameters of a conjugated `Harmonium`.
conjugatedPotential
    :: ( LegendreExponentialFamily z, ConjugatedLikelihood f x0 z0 x z )
    => Natural # AffineHarmonium f x0 z0 x z -- ^ Categorical likelihood
    -> Double -- ^ Conjugation parameters
conjugatedPotential hrm = do
    let (lkl,nz) = split hrm
        (rho0,rprms) = conjugationParameters lkl
     in potential (nz + rprms) + rho0


--- Internal ---


-- Conjugation --

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableLogDensity
    :: forall f x0 z0 x z
    . ( LegendreExponentialFamily z, ExponentialFamily x, Translation z z0
      , Translation x x0, Bilinear Natural f x0 z0 )
    => Natural # AffineHarmonium f x0 z0 x z
    -> Sample x
    -> [Double]
unnormalizedHarmoniumObservableLogDensity hrm xs =
    let (pstr,nx) = split $ transposeHarmonium hrm
        mxs = sufficientStatistic <$> xs
        nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
     in zipWith (+) nrgs $ logBaseMeasure (Proxy @x) <$> xs

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedDensities
    :: ( LegendreExponentialFamily z, ExponentialFamily x, Translation z z0
       , Translation x x0, Bilinear Natural f x0 z0 )
    => (Double, Natural # z) -- ^ Conjugation Parameters
    -> Natural # AffineHarmonium f x0 z0 x z
    -> Sample x
    -> [Double]
logConjugatedDensities (rho0,rprms) hrm x =
    let udns = unnormalizedHarmoniumObservableLogDensity hrm x
        nz = snd $ split hrm
     in subtract (potential (nz + rprms) + rho0) <$> udns

-- Mixtures --

mixtureLikelihoodConjugationParameters
    :: (KnownNat k, LegendreExponentialFamily x, Translation x x0)
    => Natural # Affine Tensor x0 x (Categorical k) -- ^ Categorical likelihood
    -> (Double, Natural # Categorical k) -- ^ Conjugation parameters
mixtureLikelihoodConjugationParameters aff =
    let (nx,nx0z0) = split aff
        rho0 = potential nx
        rprms = S.map (\nx0z0i -> subtract rho0 . potential $ nx >+> nx0z0i) $ toColumns nx0z0
     in (rho0, Point rprms)

affineMixtureToMixture
    :: (KnownNat k, Manifold x0, Manifold x, Translation x x0)
    => Natural # AffineMixture x0 x k
    -> Natural # Mixture x k
affineMixtureToMixture lmxmdl =
    let (flsk,nk) = split lmxmdl
        (nls,nlk) = split flsk
        nlsk = fromColumns . S.map (0 >+>) $ toColumns nlk
     in join (join nls nlsk) nk

mixtureToAffineMixture
    :: (KnownNat k, Manifold x, Manifold x0, Translation x x0)
    => Mean # Mixture x k
    -> Mean # AffineMixture x0 x k
mixtureToAffineMixture mxmdl =
    let (flsk,mk) = split mxmdl
        (mls,mlsk) = split flsk
        mlk = fromColumns . S.map anchor $ toColumns mlsk
     in join (join mls mlk) mk


-- Linear Gaussian Harmoniums --


linearGaussianHarmoniumConjugationParameters
    :: forall n k f .
        ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
        , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k) )
    => Natural # LinearModel f n k
    -> (Double, Natural # FullNormal k) -- ^ Conjugation parameters
linearGaussianHarmoniumConjugationParameters aff =
    let (thts,tht3) = split aff
        (tht1,tht2) = split thts
        (itht20,lndt,_) = inverseLogDeterminant . negate $ 2 .> tht2
        itht2 = -2 .> itht20
        tht21 = itht2 >.> tht1
        rho0 = -0.25 * (tht1 <.> tht21) -0.5 * lndt
        rho1 = -0.5 .> (transpose tht3 >.> tht21)
        rho2 = -0.25 .> changeOfBasis tht3 itht2
     in (rho0, join rho1 $ fromTensor rho2)

harmoniumLogBaseMeasure
    :: forall f y x z w . (ExponentialFamily z, ExponentialFamily w)
    => Proxy (AffineHarmonium f y x z w)
    -> SamplePoint (z,w)
    -> Double
harmoniumLogBaseMeasure _ (z,w) =
    logBaseMeasure (Proxy @z) z + logBaseMeasure (Proxy @w) w


-- HMOGs --

type HMoG2 f n m k = AffineMixture (FullNormal m) (LinearGaussianHarmonium f n m) k

hmog1to2
    :: ( KnownNat n, KnownNat m, KnownNat k
       , Manifold (f (MVNMean n) (MVNMean n)) )
    => c # HierarchicalMixtureOfGaussians f n m k
    -> c # HMoG2 f n m k
hmog1to2 hmog =
    let (lkl,mog) = split hmog
        (aff,cats) = split mog
        (mvn,tns) = split aff
     in join (join (lkl `join` mvn) tns) cats

hmog2to1
    :: ( KnownNat n, KnownNat m, KnownNat k
       , Manifold (f (MVNMean n) (MVNMean n)) )
    => c # HMoG2 f n m k
    -> c # HierarchicalMixtureOfGaussians f n m k
hmog2to1 hmog =
    let (bigaff,cats) = split hmog
        (iso,tns) = split bigaff
        (lkl,mvn) = split iso
     in join lkl (join (join mvn tns) cats)


--- Instances ---


--- Deriving ---

deriving instance (Manifold (Affine f x0 x z0), Manifold z)
  => Manifold (AffineHarmonium f x0 z0 x z)
deriving instance (Manifold (Affine f x0 x z0), Manifold z)
  => Product (AffineHarmonium f x0 z0 x z)

--- Harmonium ---

instance Manifold (AffineHarmonium f y x z w) => Statistical (AffineHarmonium f y x z w) where
    type SamplePoint (AffineHarmonium f y x z w) = SamplePoint (z,w)

type instance PotentialCoordinates (AffineHarmonium f y x z w) = Natural

instance ( ExponentialFamily x, ExponentialFamily z, Translation x x0
         , Translation z z0, Bilinear Mean f x0 z0 )
  => ExponentialFamily (AffineHarmonium f x0 z0 x z) where
      sufficientStatistic (z,w) =
          let mz = sufficientStatistic z
              mw = sufficientStatistic w
              my = anchor mz
              mx = anchor mw
           in joinHarmonium mz (my >.< mx) mw
      averageSufficientStatistic zws =
          let (zs,ws) = unzip zws
              mzs = sufficientStatistic <$> zs
              mws = sufficientStatistic <$> ws
              mys = anchor <$> mzs
              mxs = anchor <$> mws
           in joinHarmonium (average mzs) (mys >$< mxs) (average mws)
      logBaseMeasure = harmoniumLogBaseMeasure

instance ( ConjugatedLikelihood f x0 z0 x z, Generative Natural x
         , Generative Natural z, Map Natural f x0 z0 )
         => Generative Natural (AffineHarmonium f x0 z0 x z) where
    sample = sampleConjugated

instance ( Manifold (f y x), LegendreExponentialFamily w, ConjugatedLikelihood f y x z w )
  => Legendre (AffineHarmonium f y x z w) where
      potential = conjugatedPotential

instance ( Manifold (f y x), LegendreExponentialFamily w
         , Transition Mean Natural (AffineHarmonium f y x z w), ConjugatedLikelihood f y x z w )
  => DuallyFlat (AffineHarmonium f y x z w) where
    dualPotential mhrm =
        let nhrm = toNatural mhrm
         in mhrm <.> nhrm - potential nhrm

instance ( LegendreExponentialFamily z, ExponentialFamily x
         , Bilinear Mean f x0 z0, ConjugatedLikelihood f x0 z0 x z )
  => AbsolutelyContinuous Natural (AffineHarmonium f x0 z0 x z) where
    logDensities = exponentialFamilyLogDensities

instance ( LegendreExponentialFamily z, ConjugatedLikelihood f x0 z0 x z
         , Bilinear Natural f x0 z0 )
  => ObservablyContinuous Natural (AffineHarmonium f x0 z0 x z) where
    logObservableDensities hrm zs =
        let rho0rprms = harmoniumConjugationParameters hrm
         in logConjugatedDensities rho0rprms hrm zs

instance ( LegendreExponentialFamily z, ExponentialFamily x, SamplePoint x ~ t
         , ConjugatedLikelihood f x0 z0 x z, Bilinear Natural f x0 z0
         , Bilinear Mean f x0 z0, Transition Natural Mean (AffineHarmonium f x0 z0 x z) )
  => LogLikelihood Natural (AffineHarmonium f x0 z0 x z) t where
    logLikelihood xs hrm =
         average $ logObservableDensities hrm xs
    logLikelihoodDifferential zs hrm =
        let pxs = expectationStep zs hrm
            qxs = transition hrm
         in pxs - qxs

--- Mixture ---

instance ( KnownNat k, LegendreExponentialFamily x, Translation x x0)
  => ConjugatedLikelihood Tensor x0 (Categorical k) x (Categorical k) where
    conjugationParameters = mixtureLikelihoodConjugationParameters

instance ( KnownNat k, Manifold y, Manifold z
         , LegendreExponentialFamily z, Translation z y )
  => Transition Natural Mean (AffineMixture y z k) where
    transition mxmdl0 =
        let mxmdl = affineMixtureToMixture mxmdl0
            (nzs,nx) = splitNaturalMixture mxmdl
            mx = toMean nx
            mzs = S.map transition nzs
         in mixtureToAffineMixture $ joinMeanMixture mzs mx

instance (KnownNat k, DuallyFlatExponentialFamily x)
  => Transition Mean Natural (Mixture x k) where
    transition mhrm =
        let (mxs,mk) = splitMeanMixture mhrm
            nk = transition mk
            nxs = S.map transition mxs
         in joinNaturalMixture nxs nk

instance (KnownNat k, LegendreExponentialFamily x, Transition Mean Source x)
  => Transition Mean Source (Mixture x k) where
    transition mhrm =
        let (mxs,mk) = splitMeanMixture mhrm
            nk = transition mk
            nxs = S.map transition mxs
         in joinSourceMixture nxs nk

instance (KnownNat k, LegendreExponentialFamily x, Transition Natural Source x)
  => Transition Natural Source (Mixture x k) where
    transition nhrm =
        let (nxs,nk) = splitNaturalMixture nhrm
            sk = transition nk
            sxs = S.map transition nxs
         in joinSourceMixture sxs sk

instance (KnownNat k, LegendreExponentialFamily x, Transition Source Natural x)
  => Transition Source Natural (Mixture x k) where
    transition shrm =
        let (sxs,sk) = splitSourceMixture shrm
            nk = transition sk
            nxs = S.map transition sxs
         in joinNaturalMixture nxs nk


--- Linear Guassian Harmonium ---


instance ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
         , ExponentialFamily (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k) )
    => ConjugatedLikelihood Tensor (MVNMean n) (MVNMean k)
    (MultivariateNormal f n) (FullNormal k) where
        conjugationParameters = linearGaussianHarmoniumConjugationParameters

instance ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
  => Transition Natural Source (LinearGaussianHarmonium f n k) where
      transition nlgh =
          let (nfxz,nz) = split nlgh
              (nx,nvrxz) = split nfxz
              (nmux,nvrx) = split nx
              (nmuz,nvrz) = split nz
              (svrx0,svrxz0,svrz0) = blockSymmetricMatrixInversion nvrx (2 /> nvrxz) (toTensor nvrz)
              svrx = -0.5 .> svrx0
              svrxz = -0.5 .> svrxz0
              svrz = -0.5 .> svrz0
              smux = svrx >.> nmux + svrxz >.> nmuz
              smuz = svrz >.> nmuz + transpose svrxz >.> nmux
              sx = join smux $ fromTensor svrx
              sz = join smuz $ fromTensor svrz
              sfxz = join sx $ fromTensor svrxz
              slgh = join sfxz sz
           in breakChart slgh


instance (KnownNat n, KnownNat k, Bilinear Source f (MVNMean n) (MVNMean n)
         , Transition Source Mean (MultivariateNormal f n) )
  => Transition Source Mean (LinearGaussianHarmonium f n k) where
      transition slgh =
          let (sfxz,sz) = split slgh
              (sx,svrxz) = split sfxz
              smux = fst $ split sx
              smuz = fst $ split sz
              mvrxz = breakChart $ svrxz + smux >.< smuz
              mx = toMean sx
              mz = toMean sz
              mfxz = join mx mvrxz
           in join mfxz mz

instance ( KnownNat n, KnownNat k, Manifold (f (MVNMean n) (MVNMean n))
         , Transition Mean Source (MultivariateNormal f n) )
  => Transition Mean Source (LinearGaussianHarmonium f n k) where
      transition mlgh =
          let (mfxz,mz) = split mlgh
              (mx,mvrxz) = split mfxz
              mmux = fst $ split mx
              mmuz = fst $ split mz
              svrxz = breakChart $ mvrxz - mmux >.< mmuz
              sx = toSource mx
              sz = toSource mz
              sfxz = join sx svrxz
           in join sfxz sz

instance ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
         , Bilinear Source f (MVNMean n) (MVNMean n)
         , Transition Source Mean (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
         => Transition Natural Mean (LinearGaussianHarmonium f n k) where
    transition = toMean . toSource

instance ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
         , Transition Mean Source (MultivariateNormal f n)
         , Square Source f (MVNMean n), ExponentialFamily (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
  => Transition Source Natural (LinearGaussianHarmonium f n k) where
    transition slgh =
        let (sx,svrxz,sz) = splitHarmonium slgh
            (smux,svrx) = split sx
            (smuz,svrz) = split sz
            invsvrz = inverse svrz
            nvrx0 = inverse $ svrx - fromTensor (changeOfBasis (transpose svrxz) invsvrz)
            nvrx = -0.5 .> nvrx0
            nxz = dualComposition nvrx0 svrxz invsvrz
            nmux = nvrx0 >.> smux -nxz >.> smuz
            nmvn = join nmux nvrx
            nlkl = join nmvn nxz
            nz = toNatural sz
         in joinConjugatedHarmonium (breakChart nlkl) nz

instance ( KnownNat n, KnownNat k, Square Natural f (MVNMean n), Square Source f (MVNMean n)
         , Transition Mean Source (MultivariateNormal f n)
         , ExponentialFamily (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
  => Transition Mean Natural (LinearGaussianHarmonium f n k) where
    transition = toNatural . toSource


--- Hierarchical Mixture of Gaussians ---

instance ( KnownNat n, KnownNat m, KnownNat k, ExponentialFamily (MultivariateNormal f n)
         , Square Natural f (MVNMean n), Translation (Mixture (FullNormal m) k) (MVNMean m)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean m) )
         => ConjugatedLikelihood
    Tensor (MVNMean n) (MVNMean m) (MultivariateNormal f n) (Mixture (FullNormal m) k)
        where conjugationParameters lm =
                let (rho0,rprms) = conjugationParameters lm
                 in (rho0,join (join rprms 0) 0)

instance (KnownNat n, KnownNat m, Manifold (f (MVNMean n) (MVNMean n)))
  => Translation (LinearGaussianHarmonium f n m) (FullNormal m) where
      (>+>) hrm ny =
          let (nz,nyx,nw) = splitHarmonium hrm
           in joinHarmonium nz nyx (nw >+> ny)
      anchor hrm =
          let (_,_,nw) = splitHarmonium hrm
           in anchor nw


instance (KnownNat n, KnownNat k) => Translation (Mixture (FullNormal n) k) (MVNMean n) where
      (>+>) hrm ny =
          let (nz,nyx,nw) = splitHarmonium hrm
           in joinHarmonium (nz >+> ny) nyx nw
      anchor hrm =
          let (nz,_,_) = splitHarmonium hrm
           in anchor nz

instance ( KnownNat n, KnownNat m, KnownNat k, ExponentialFamily (MultivariateNormal f n)
         , Square Natural f (MVNMean n), Translation (Mixture (FullNormal m) k) (MVNMean m)
         , Bilinear Source f (MVNMean n) (MVNMean n), Transition Source Mean (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean m) )
  => Transition Natural Mean (HierarchicalMixtureOfGaussians f n m k) where
      transition = hmog2to1 . transition . hmog1to2

instance ( KnownNat n, KnownNat k, KnownNat m, Manifold (f (MVNMean n) (MVNMean n))
         , Transition Mean Source (MultivariateNormal f n) )
  => Transition Mean Source (HierarchicalMixtureOfGaussians f n m k) where
      transition mhmog =
          let (mfxy,myz) = split mhmog
              mmuy = fst . split . snd . split $ transposeHarmonium myz
              syz = toSource myz
              (mx,mvrxy) = split mfxy
              mmux = fst $ split mx
              svrxy = breakChart $ mvrxy - mmux >.< mmuy
              sx = toSource mx
              sfxy = join sx svrxy
           in join sfxy syz

instance ( KnownNat n, KnownNat k, KnownNat m, Manifold (f (MVNMean n) (MVNMean n))
         , Square Source f (MVNMean n), ExponentialFamily (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean m)
         , Transition Mean Source (MultivariateNormal f n), Square Natural f (MVNMean n) )
  => Transition Source Natural (HierarchicalMixtureOfGaussians f n m k) where
    transition shmog =
        let (sfxy,syz) = split shmog
            (sx,svrxy) = split sfxy
            (smux,svrx) = split sx
            (smuy,svry) = split . snd . split $ transposeHarmonium syz
            invsvry = inverse svry
            nvrx0 = inverse $ svrx - fromTensor (changeOfBasis (transpose svrxy) invsvry)
            nvrx = -0.5 .> nvrx0
            nxy = dualComposition nvrx0 svrxy invsvry
            nmux = nvrx0 >.> smux -nxy >.> smuy
            nmvn = join nmux nvrx
            nlkl = join nmvn nxy
            nyz = toNatural syz
         in joinConjugatedHarmonium (breakChart nlkl) nyz

instance ( KnownNat n, KnownNat k, KnownNat m, Manifold (f (MVNMean n) (MVNMean n))
         , Square Source f (MVNMean n), ExponentialFamily (MultivariateNormal f n)
         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean m)
         , Transition Mean Source (MultivariateNormal f n), Square Natural f (MVNMean n) )
  => Transition Mean Natural (HierarchicalMixtureOfGaussians f n m k) where
    transition = toNatural . toSource


--instance ( KnownNat n, KnownNat m, KnownNat k
--         , Manifold (f (MVNMean n) (MVNMean n)) )
--  => Transition Mean Natural (HierarchicalMixtureOfGaussians f n m k) where
--      transition = hmog2to1 . transition . hmog1to2

--instance (KnownNat n, KnownNat m, KnownNat k)
--  => Generative Natural (IsotropicHMOG n m k) where
--      sample n hmog = do
--          let (pca,mog) = splitConjugatedHarmonium hmog
--          yzs <- sample n mog
--          xs <- mapM samplePoint $ pca >$>* (fst <$> yzs)
--          return $ zip xs yzs
--


--- Graveyard

--instance (KnownNat n, KnownNat k) => Transition Natural Mean (IsotropicGaussianHarmonium n k) where
--      transition = linearGaussianHarmoniumToIsotropic . transition . isotropicGaussianHarmoniumToLinear
--
--instance (KnownNat n, KnownNat k) => Transition Natural Mean (DiagonalGaussianHarmonium n k) where
--      transition = linearGaussianHarmoniumToDiagonal . transition . diagonalGaussianHarmoniumToLinear
--
--
--instance (KnownNat n, KnownNat k) => Transition Mean Natural (LinearGaussianHarmonium n k) where
--      transition = naturalJointToLinearGaussianHarmonium . transition
--        . meanLinearGaussianHarmoniumToJoint
--
--sourcePCAMaximizationStep'
--    :: forall n k . (KnownNat n, KnownNat k)
--    => Mean # LinearGaussianHarmonium n k
--    -> Source # PrincipleComponentAnalysis n k
--sourcePCAMaximizationStep' hrm =
--    let (mz,mzx,mx) = splitHarmonium hrm
--        (muz,etaz) = splitMeanMultivariateNormal mz
--        (mux,etax) = splitMeanMultivariateNormal mx
--        outrs = toMatrix mzx - S.outerProduct muz mux
--        wmtx = S.matrixMatrixMultiply outrs $ S.inverse etax
--        wmtxtr = S.transpose wmtx
--        n = fromIntegral $ natVal (Proxy @n)
--        zcvr = etaz - S.outerProduct muz muz
--        vr = S.trace $ zcvr - 2*S.matrixMatrixMultiply wmtx (S.transpose outrs)
--            + S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtx etax) wmtxtr
--        iso = join (Point muz) $ singleton vr / n
--     in join iso $ fromMatrix wmtx


--linearGaussianHarmoniumConjugationParameters
--    :: (KnownNat n, KnownNat k)
--    => Natural # Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)
--    -> (Double, Natural # MultivariateNormal k) -- ^ Conjugation parameters
--linearGaussianHarmoniumConjugationParameters aff =
--    let (thts,tht30) = split aff
--        (tht1,tht2) = splitNaturalMultivariateNormal thts
--        tht3 = toMatrix tht30
--        ttht3 = S.transpose tht3
--        itht2 = S.pseudoInverse tht2
--        rho0 = -0.25 * tht1 `S.dotProduct` (itht2 `S.matrixVectorMultiply` tht1)
--            -0.5 * (log . S.determinant . negate $ 2*tht2)
--        rho1 = -0.5 * ttht3 `S.matrixVectorMultiply` (itht2 `S.matrixVectorMultiply` tht1)
--        rho2 = -0.25 * ttht3 `S.matrixMatrixMultiply` (itht2 `S.matrixMatrixMultiply` tht3)
--     in (rho0, joinNaturalMultivariateNormal rho1 rho2)


--univariateToLinearGaussianHarmonium
--    :: c # AffineHarmonium Tensor NormalMean NormalMean Normal Normal
--    -> c # LinearGaussianHarmonium 1 1
--univariateToLinearGaussianHarmonium hrm =
--    let (z,zx,x) = splitHarmonium hrm
--     in joinHarmonium (breakChart z) (breakChart zx) (breakChart x)
--
--linearGaussianHarmoniumToUnivariate
--    :: c # LinearGaussianHarmonium 1 1
--    -> c # AffineHarmonium Tensor NormalMean NormalMean Normal Normal
--linearGaussianHarmoniumToUnivariate hrm =
--    let (z,zx,x) = splitHarmonium hrm
--     in joinHarmonium (breakChart z) (breakChart zx) (breakChart x)
--
--univariateToLinearModel
--    :: Natural # Affine Tensor NormalMean Normal NormalMean
--    -> Natural # Affine Tensor (MVNMean 1) (MultivariateNormal 1) (MVNMean 1)
--univariateToLinearModel aff =
--    let (z,zx) = split aff
--     in join (breakChart z) (breakChart zx)
--
--naturalLinearGaussianHarmoniumToJoint
--    :: (KnownNat n, KnownNat k)
--    => Natural # LinearGaussianHarmonium n k
--    -> Natural # MultivariateNormal (n+k)
--naturalLinearGaussianHarmoniumToJoint hrm =
--    let (z,zx,x) = splitHarmonium hrm
--        zxmtx = toMatrix zx/2
--        mvnz = splitNaturalMultivariateNormal z
--        mvnx = splitNaturalMultivariateNormal x
--        (mu,cvr) = fromLinearGaussianHarmonium0 mvnz zxmtx mvnx
--     in joinNaturalMultivariateNormal mu cvr
--
--naturalJointToLinearGaussianHarmonium
--    :: (KnownNat n, KnownNat k)
--    => Natural # MultivariateNormal (n+k)
--    -> Natural # LinearGaussianHarmonium n k
--naturalJointToLinearGaussianHarmonium mvn =
--    let (mu,cvr) = splitNaturalMultivariateNormal mvn
--        ((muz,cvrz),zxmtx,(mux,cvrx)) = toLinearGaussianHarmonium0 mu cvr
--        zx = 2*fromMatrix zxmtx
--        z = joinNaturalMultivariateNormal muz cvrz
--        x = joinNaturalMultivariateNormal mux cvrx
--     in joinHarmonium z zx x
--
--meanLinearGaussianHarmoniumToJoint
--    :: (KnownNat n, KnownNat k)
--    => Mean # LinearGaussianHarmonium n k
--    -> Mean # MultivariateNormal (n+k)
--meanLinearGaussianHarmoniumToJoint hrm =
--    let (z,zx,x) = splitHarmonium hrm
--        zxmtx = toMatrix zx
--        mvnz = splitMeanMultivariateNormal z
--        mvnx = splitMeanMultivariateNormal x
--        (mu,cvr) = fromLinearGaussianHarmonium0 mvnz zxmtx mvnx
--     in joinMeanMultivariateNormal mu cvr
--
--meanJointToLinearGaussianHarmonium
--    :: (KnownNat n, KnownNat k)
--    => Mean # MultivariateNormal (n+k)
--    -> Mean # LinearGaussianHarmonium n k
--meanJointToLinearGaussianHarmonium mvn =
--    let (mu,cvr) = splitMeanMultivariateNormal mvn
--        ((muz,cvrz),zxmtx,(mux,cvrx)) = toLinearGaussianHarmonium0 mu cvr
--        zx = fromMatrix zxmtx
--        z = joinMeanMultivariateNormal muz cvrz
--        x = joinMeanMultivariateNormal mux cvrx
--     in joinHarmonium z zx x
--
--fromLinearGaussianHarmonium0
--    :: (KnownNat n, KnownNat k)
--    => (S.Vector n Double, S.Matrix n n Double)
--    -> S.Matrix n k Double
--    -> (S.Vector k Double, S.Matrix k k Double)
--    -> (S.Vector (n+k) Double, S.Matrix (n+k) (n+k) Double)
--fromLinearGaussianHarmonium0 (muz,cvrz) zxmtx (mux,cvrx) =
--    let mu = muz S.++ mux
--        top = S.horizontalConcat cvrz zxmtx
--        btm = S.horizontalConcat (S.transpose zxmtx) cvrx
--     in (mu, S.verticalConcat top btm)
--
--toLinearGaussianHarmonium0
--    :: (KnownNat n, KnownNat k)
--    => S.Vector (n+k) Double
--    -> S.Matrix (n+k) (n+k) Double
--    -> ( (S.Vector n Double, S.Matrix n n Double)
--       , S.Matrix n k Double
--       , (S.Vector k Double, S.Matrix k k Double) )
--toLinearGaussianHarmonium0 mu cvr =
--    let (muz,mux) = S.splitAt mu
--        (tops,btms) = S.splitAt $ S.toRows cvr
--        (cvrzs,zxmtxs) = S.splitAt . S.toColumns $ S.fromRows tops
--        cvrz = S.fromColumns cvrzs
--        zxmtx = S.fromColumns zxmtxs
--        cvrx = S.fromColumns . S.drop . S.toColumns $ S.fromRows btms
--     in ((muz,cvrz),zxmtx,(mux,cvrx))

--isotropicGaussianHarmoniumToLinear
--    :: (KnownNat n, KnownNat k)
--    => Natural # IsotropicGaussianHarmonium n k
--    -> Natural # LinearGaussianHarmonium n k
--isotropicGaussianHarmoniumToLinear isohrm =
--    let (lkl,prr) = split isohrm
--        (iso,tns) = split lkl
--        lkl' = join (isotropicNormalToFull iso) tns
--     in join lkl' prr
--
--linearGaussianHarmoniumToIsotropic
--    :: (KnownNat n, KnownNat k)
--    => Mean # LinearGaussianHarmonium n k
--    -> Mean # IsotropicGaussianHarmonium n k
--linearGaussianHarmoniumToIsotropic lnrhrm =
--    let (lkl,prr) = split lnrhrm
--        (lnr,tns) = split lkl
--        lkl' = join (fullNormalToIsotropic lnr) tns
--     in join lkl' prr
--
--diagonalGaussianHarmoniumToLinear
--    :: (KnownNat n, KnownNat k)
--    => Natural # DiagonalGaussianHarmonium n k
--    -> Natural # LinearGaussianHarmonium n k
--diagonalGaussianHarmoniumToLinear isohrm =
--    let (lkl,prr) = split isohrm
--        (iso,tns) = split lkl
--        lkl' = join (diagonalNormalToFull iso) tns
--     in join lkl' prr
--
--linearGaussianHarmoniumToDiagonal
--    :: (KnownNat n, KnownNat k)
--    => Mean # LinearGaussianHarmonium n k
--    -> Mean # DiagonalGaussianHarmonium n k
--linearGaussianHarmoniumToDiagonal lnrhrm =
--    let (lkl,prr) = split lnrhrm
--        (lnr,tns) = split lkl
--        lkl' = join (fullNormalToDiagonal lnr) tns
--     in join lkl' prr

--type IsotropicHMOG2 n m k = AffineMixture (MultivariateNormal m) (IsotropicGaussianHarmonium n m) k
--type DiagonalHMOG2 n m k = AffineMixture (MultivariateNormal m) (DiagonalGaussianHarmonium n m) k

--hmog1to2'
--    :: ( KnownNat n, KnownNat m, KnownNat k )
--    => c # DiagonalHMOG n m k
--    -> c # DiagonalHMOG2 n m k
--hmog1to2' hmog =
--    let (lkl,mog) = split hmog
--        (aff,cats) = split mog
--        (mvn,tns) = split aff
--     in join (join (lkl `join` mvn) tns) cats
--
--hmog2to1'
--    :: ( KnownNat n, KnownNat m, KnownNat k )
--    => c # DiagonalHMOG2 n m k
--    -> c # DiagonalHMOG n m k
--hmog2to1' hmog =
--    let (bigaff,cats) = split hmog
--        (iso,tns) = split bigaff
--        (lkl,mvn) = split iso
--     in join lkl (join (join mvn tns) cats)



--instance ConjugatedLikelihood Tensor NormalMean NormalMean Normal Normal where
--    conjugationParameters aff =
--        let rprms :: Natural # MultivariateNormal 1
--            (rho0,rprms) = conjugationParameters $ univariateToLinearModel aff
--         in (rho0,breakChart rprms)


--instance ( KnownNat k, LegendreExponentialFamily z
--         , Generative Natural z, Manifold (Mixture z k) )
--         => Generative Natural (Mixture z k) where
--    sample = sampleConjugated

--instance (KnownNat k, LegendreExponentialFamily z)
--  => Transition Natural Mean (Mixture z k) where
--    transition nhrm =
--        let (nzs,nx) = splitNaturalMixture nhrm
--            mx = toMean nx
--            mzs = S.map transition nzs
--         in joinMeanMixture mzs mx

--instance ( KnownNat k, Manifold y, Manifold z, LegendreExponentialFamily z
--         , Generative Natural z, Translation z y )
--  => Generative Natural (AffineMixture y z k) where
--      sample n = sampleConjugated n . affineMixtureToMixture

--instance Transition Natural Mean
--  (AffineHarmonium Tensor NormalMean NormalMean Normal Normal) where
--      transition = linearGaussianHarmoniumToUnivariate . transition . univariateToLinearGaussianHarmonium
--
--instance Transition Mean Natural
--  (AffineHarmonium Tensor NormalMean NormalMean Normal Normal) where
--      transition =  linearGaussianHarmoniumToUnivariate . transition . univariateToLinearGaussianHarmonium

--instance Generative Natural (LinearGaussianHarmonium f n k) where
--    sample n lgh = do
--        let (aff,prr) = splitConjugatedHarmonium lgh
--        zs <- sample n prr
--        let (mus,sgms) = unzip $ split . toSource <$> aff >$>* zs
--            nrm0 :: Source # MultivariateNormal f n
--            nrm0 = join 0 $ head sgms
--        xs0 <- sample n nrm0
--        let xs = zipWith (+) xs0 $ coordinates <$> mus
--        return $ zip xs zs

--instance (KnownNat n, KnownNat k) => Generative Natural (IsotropicGaussianHarmonium n k) where
--    sample n lgh = do
--        let (aff,prr) = splitConjugatedHarmonium lgh
--        zs <- sample n prr
--        xs <- mapM samplePoint $ aff >$>* zs
--        return $ zip xs zs

--instance (KnownNat n, KnownNat k) => Transition Mean Natural (IsotropicGaussianHarmonium n k) where
--      transition igh =
--          let prr = snd $ split igh
--              pca = transition $ sourcePCAMaximizationStep igh
--           in joinConjugatedHarmonium pca $ transition prr


--type instance PotentialCoordinates (Mixture z k) = Natural
--
--instance (KnownNat k, LegendreExponentialFamily z) => Legendre (Mixture z k) where
--      potential = conjugatedPotential

--instance (KnownNat n, KnownNat k) => Translation (Mixture (MultivariateNormal n) k) (MVNMean n) where
--      (>+>) hrm ny =
--          let (nz,nyx,nw) = splitHarmonium hrm
--           in joinHarmonium (nz >+> ny) nyx nw
--      anchor hrm =
--          let (nz,_,_) = splitHarmonium hrm
--           in anchor nz

--instance (KnownNat n, KnownNat m)
--  => Translation (DiagonalGaussianHarmonium n m) (MultivariateNormal m) where
--      (>+>) hrm ny =
--          let (nz,nyx,nw) = splitHarmonium hrm
--           in joinHarmonium nz nyx (nw >+> ny)
--      anchor hrm =
--          let (_,_,nw) = splitHarmonium hrm
--           in anchor nw
--
--instance (KnownNat n, KnownNat m, KnownNat k) => ConjugatedLikelihood
--    Tensor (MVNMean n) (MVNMean m) (DiagonalNormal n) (Mixture (MultivariateNormal m) k)
--        where conjugationParameters pca =
--                let (rho0,rprms) = conjugationParameters pca
--                 in (rho0,join (join rprms 0) 0)
--
--instance (KnownNat n, KnownNat m)
--  => Translation (IsotropicGaussianHarmonium n m) (MultivariateNormal m) where
--      (>+>) hrm ny =
--          let (nz,nyx,nw) = splitHarmonium hrm
--           in joinHarmonium nz nyx (nw >+> ny)
--      anchor hrm =
--          let (_,_,nw) = splitHarmonium hrm
--           in anchor nw
--
--instance (KnownNat n, KnownNat m, KnownNat k)
--  => Transition Natural Mean (IsotropicHMOG n m k) where
--      transition = hmog2to1 . transition . hmog1to2
--
--instance (KnownNat n, KnownNat m, KnownNat k)
--  => Generative Natural (IsotropicHMOG n m k) where
--      sample n hmog = do
--          let (pca,mog) = splitConjugatedHarmonium hmog
--          yzs <- sample n mog
--          xs <- mapM samplePoint $ pca >$>* (fst <$> yzs)
--          return $ zip xs yzs
--
--instance (KnownNat n, KnownNat m, KnownNat k)
--  => Transition Natural Mean (DiagonalHMOG n m k) where
--      transition = hmog2to1' . transition . hmog1to2'
--
--instance (KnownNat n, KnownNat m, KnownNat k)
--  => Generative Natural (DiagonalHMOG n m k) where
--      sample n hmog = do
--          let (pca,mog) = splitConjugatedHarmonium hmog
--          yzs <- sample n mog
--          xs <- mapM samplePoint $ pca >$>* (fst <$> yzs)
--          return $ zip xs yzs


--instance ( KnownNat k, Square Source f (MVNMean n), Bilinear Natural f (MVNMean n) (MVNMean n)
--         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
--         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
--         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
--  => Transition Source Natural (LinearGaussianHarmonium f n k) where
--      transition slgh =
--          let (sfxz,sz) = split slgh
--              (sx,svrxz) = split sfxz
--              (smux,svrx) = split sx
--              (smuz,svrz) = split sz
--              (nvrx0,nvrxz0,nvrz0) = blockSymmetricMatrixInversion svrx svrxz (toTensor svrz)
--              nvrx1 :: Natural # Tensor (MVNMean n) (MVNMean n)
--              nvrx1 = breakChart $ nvrx0
--              nvrxz1 :: Natural # Tensor (MVNMean n) (MVNMean k)
--              nvrxz1 = breakChart $ nvrxz0
--              nvrz1 :: Natural # Tensor (MVNMean k) (MVNMean k)
--              nvrz1 = breakChart $ nvrz0
--              nvrx2 = fromTensor nvrx1
--              nvrz2 = fromTensor nvrz1
--              nmux = nvrx2 >.> breakChart smux + nvrxz1 >.> breakChart smuz
--              nmuz = breakChart $ nvrz2 >.> breakChart smuz + transpose nvrxz1 >.> breakChart smux
--              nx = join nmux (0.5 .> nvrx2)
--              nz = join nmuz (0.5 .> nvrz2)
--              nfxz = join nx (-nvrxz1)
--           in join nfxz nz

--instance ( KnownNat n, KnownNat k)
--  => Transition Source Natural (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)) where
--    transition sfa =
--        let (snrms,smtx) = split sfa
--            (smu,ssg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates
--                $ splitReplicated snrms
--            invsg = recip ssg
--            nmu = invsg * smu
--            nsg = -0.5 * invsg
--            nmtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) $ toMatrix smtx
--            nnrms = joinReplicated $ S.zipWith (curry fromTuple) nmu nsg
--         in join nnrms $ fromMatrix nmtx


--instance ( KnownNat k, Square Source f (MVNMean n), Bilinear Natural f (MVNMean n) (MVNMean n)
--         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
--         , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
--         , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
--  => Transition Source Natural (LinearGaussianHarmonium f n k) where
--      transition slgh =
--          let (sfxz,sz) = split slgh
--              (sx,svrxz) = split sfxz
--              (smux,svrx) = split sx
--              (smuz,svrz) = split sz
--              (nvrx0,nvrxz0,nvrz0) = blockSymmetricMatrixInversion svrx svrxz (toTensor svrz)
--              nmux :: Natural # MVNMean n
--              nmux = breakChart $ nvrx0 >.> smux + nvrxz0 >.> smuz
--              nmuz :: Natural # MVNMean k
--              nmuz = breakChart $ nvrz0 >.> smuz + transpose nvrxz0 >.> smux
--              nvrx :: Natural # Tensor (MVNMean n) (MVNMean n)
--              nvrx = breakChart $ -0.5 .> nvrx0
--              nvrxz :: Natural # Tensor (MVNMean n) (MVNMean k)
--              nvrxz = breakChart $ -nvrxz0
--              nvrz :: Natural # Tensor (MVNMean k) (MVNMean k)
--              nvrz = breakChart $ -0.5 .> nvrz0
--              nx = join nmux (fromTensor nvrx)
--              nz = join nmuz (fromTensor nvrz)
--              nfxz = join nx nvrxz
--           in join nfxz nz



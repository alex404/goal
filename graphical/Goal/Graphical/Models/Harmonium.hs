{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    TypeApplications,
    UndecidableInstances,
    NoStarIsType,
    GeneralizedNewtypeDeriving,
    StandaloneDeriving,
    ScopedTypeVariables,
    ExplicitNamespaces,
    TypeOperators,
    KindSignatures,
    DataKinds,
    RankNTypes,
    TypeFamilies,
    FlexibleContexts,
    MultiParamTypeClasses,
    ConstraintKinds,
    FlexibleInstances
#-}
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
    , IsotropicGaussianHarmonium
    -- ** Conjugated Harmoniums
    , ConjugatedLikelihood (conjugationParameters)
    , joinConjugatedHarmonium
    , splitConjugatedHarmonium
    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models

import qualified Goal.Core.Vector.Storable as S


--- Types ---


-- | A 2-layer harmonium.
newtype AffineHarmonium f y x z w = AffineHarmonium (Affine f y z x, w)

deriving instance (Manifold z, Manifold (f y x), Manifold w)
  => Manifold (AffineHarmonium f y x z w)
deriving instance (Manifold z, Manifold (f y x), Manifold w)
  => Product (AffineHarmonium f y x z w)

type Harmonium f z w = AffineHarmonium f z w z w

type instance Observation (AffineHarmonium f y x z w) = SamplePoint z

-- | A 'Mixture' model is simply a 'AffineHarmonium' where the latent variable is
-- 'Categorical'.
type Mixture z k = Harmonium Tensor z (Categorical k)

-- | A 'Mixture' where only a subset of the component parameters are mixed.
type AffineMixture y z k =
    AffineHarmonium Tensor y (Categorical k) z (Categorical k)

type LinearGaussianHarmonium n k =
    AffineHarmonium Tensor (MVNMean n) (MVNMean k) (MultivariateNormal n) (MultivariateNormal k)

type IsotropicGaussianHarmonium n k =
    AffineHarmonium Tensor (MVNMean n) (MVNMean k) (IsotropicNormal n) (MultivariateNormal k)


--- Classes ---


-- | The conjugation parameters of a conjugated likelihood.
class ( ExponentialFamily z, ExponentialFamily w, Map Natural f y x
      , Translation z y , Translation w x
      , SamplePoint y ~ SamplePoint z, SamplePoint x ~ SamplePoint w )
  => ConjugatedLikelihood f y x z w where
    conjugationParameters
        :: Natural # Affine f y z x -- ^ Categorical likelihood
        -> (Double, Natural # w) -- ^ Conjugation parameters


--- Functions ---


-- Construction --

-- | Creates a 'Harmonium' from component parameters.
joinHarmonium
    :: (Manifold w, Manifold z, Manifold (f y x))
    => c # z -- ^ Visible layer biases
    -> c # f y x -- ^ ^ Interaction parameters
    -> c # w -- ^ Hidden layer Biases
    -> c # AffineHarmonium f y x z w -- ^ Harmonium
joinHarmonium nz nyx = join (join nz nyx)

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium
    :: (Manifold z, Manifold (f y x), Manifold w)
    => c # AffineHarmonium f y x z w -- ^ Harmonium
    -> (c # z, c # f y x, c # w) -- ^ Biases and interaction parameters
splitHarmonium hrm =
    let (fzx,nw) = split hrm
        (nz,nyx) = split fzx
     in (nz,nyx,nw)

-- | Build a mixture model in source coordinates.
joinSourceMixture
    :: (KnownNat k, Manifold z)
    => S.Vector (k+1) (Source # z) -- ^ Mixture components
    -> Source # Categorical k -- ^ Weights
    -> Source # Mixture z k
joinSourceMixture szs sx =
    let (sz,szs') = S.splitAt szs
        aff = join (S.head sz) (fromColumns szs')
     in join aff sx

-- | Build a mixture model in source coordinates.
splitSourceMixture
    :: (KnownNat k, Manifold z)
    => Source # Mixture z k
    -> (S.Vector (k+1) (Source # z), Source # Categorical k)
splitSourceMixture mxmdl =
    let (aff,sx) = split mxmdl
        (sz0,szs0') = split aff
     in (S.cons sz0 $ toColumns szs0' ,sx)

-- | Build a mixture model in mean coordinates.
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
     in joinHarmonium mz mzx mx

-- | Split a mixture model in mean coordinates.
splitMeanMixture
    :: ( KnownNat k, DuallyFlatExponentialFamily z )
    => Mean # Mixture z k
    -> (S.Vector (k+1) (Mean # z), Mean # Categorical k)
splitMeanMixture hrm =
    let (mz,mzx,mx) = splitHarmonium hrm
        twmzs = toRows $ transpose mzx
        wmzs = S.cons (mz - S.foldr (+) 0 twmzs) twmzs
        wghts = categoricalWeights mx
        mzs = S.zipWith (/>) wghts wmzs
     in (mzs,mx)

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
        affzx = join nz nzx
        rprms = snd $ conjugationParameters affzx
        nx = nx0 - rprms
     in joinHarmonium nz nzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitNaturalMixture
    :: forall k z . ( KnownNat k, LegendreExponentialFamily z )
    => Natural # Mixture z k -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # z), Natural # Categorical k) -- ^ (components, weights)
splitNaturalMixture hrm =
    let (nz,nzx,nx) = splitHarmonium hrm
        affzx = join nz nzx
        rprms = snd $ conjugationParameters affzx
        nx0 = nx + rprms
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (+ nz) nzs
     in (S.cons nz nzs0',nx0)


-- Manipulation --

-- | Swap the biases and 'transpose' the interaction parameters of the given 'Harmonium'.
transposeHarmonium
    :: (Bilinear f y x, Manifold z, Manifold w)
    => c # AffineHarmonium f y x z w
    -> c # AffineHarmonium f x y w z
transposeHarmonium hrm =
        let (nz,nyx,nw) = splitHarmonium hrm
         in joinHarmonium nw (transpose nyx) nz

-- Evaluation --

-- | Computes the joint expectations of a harmonium based on a sample from the
-- observable layer.
expectationStep
    :: ( ExponentialFamily z, Map Natural f x y, Bilinear f y x
       , Translation z y, Translation w x, LegendreExponentialFamily w )
    => Sample z -- ^ Model Samples
    -> Natural # AffineHarmonium f y x z w -- ^ Harmonium
    -> Mean # AffineHarmonium f y x z w -- ^ Harmonium expected sufficient statistics
expectationStep zs hrm =
    let mzs = sufficientStatistic <$> zs
        mys = anchor <$> mzs
        pstr = fst . split $ transposeHarmonium hrm
        mws = transition <$> pstr >$> mys
        mxs = anchor <$> mws
        myx = (>$<) mys mxs
     in joinHarmonium (average mzs) myx $ average mws

---- Sampling --

-- | Initialize a Gibbs chain from a set of observations.
initialPass
    :: forall f x y z w
    . ( ExponentialFamily z, Map Natural f x y, Manifold w
      , SamplePoint y ~ SamplePoint z, Translation w x, Generative Natural w
      , ExponentialFamily y, Bilinear f y x, LegendreExponentialFamily w )
    => Natural # AffineHarmonium f y x z w -- ^ Harmonium
    -> Sample z -- ^ Model Samples
    -> Random (Sample (z, w))
initialPass hrm zs = do
    let pstr = fst . split $ transposeHarmonium hrm
    ws <- mapM samplePoint $ pstr >$>* zs
    return $ zip zs ws

-- | Update a 'Sample' with Gibbs sampling.
gibbsPass
    :: ( ExponentialFamily z, Map Natural f x y, Translation z y
       , Translation w x, SamplePoint z ~ SamplePoint y, Generative Natural w
       , ExponentialFamily y, SamplePoint x ~ SamplePoint w, Bilinear f y x
       , Map Natural f y x, ExponentialFamily x, Generative Natural z )
    => Natural # AffineHarmonium f y x z w -- ^ Harmonium
    -> Sample (z, w)
    -> Random (Sample (z, w))
gibbsPass hrm zws = do
    let ws = snd <$> zws
        pstr = fst . split $ transposeHarmonium hrm
        lkl = fst $ split hrm
    zs' <- mapM samplePoint $ lkl >$>* ws
    ws' <- mapM samplePoint $ pstr >$>* zs'
    return $ zip zs' ws'

-- Conjugation --

-- | The conjugation parameters of a conjugated `Harmonium`.
harmoniumConjugationParameters
    :: ConjugatedLikelihood f y x z w
    => Natural # AffineHarmonium f y x z w -- ^ Categorical likelihood
    -> (Double, Natural # w) -- ^ Conjugation parameters
harmoniumConjugationParameters hrm =
    conjugationParameters . fst $ split hrm

-- | The conjugation parameters of a conjugated `Harmonium`.
splitConjugatedHarmonium
    :: ConjugatedLikelihood f y x z w
    => Natural # AffineHarmonium f y x z w
    -> (Natural # Affine f y z x, Natural # w) -- ^ Conjugation parameters
splitConjugatedHarmonium hrm =
    let (lkl,nw) = split hrm
        cw = snd $ conjugationParameters lkl
     in (lkl,nw + cw)

-- | The conjugation parameters of a conjugated `Harmonium`.
joinConjugatedHarmonium
    :: ConjugatedLikelihood f y x z w
    => Natural # Affine f y z x -- ^ Conjugation parameters
    -> Natural # w
    -> Natural # AffineHarmonium f y x z w -- ^ Categorical likelihood
joinConjugatedHarmonium lkl nw =
    let cw = snd $ conjugationParameters lkl
     in join lkl $ nw - cw

-- | The conjugation parameters of a conjugated `Harmonium`.
sampleConjugated
    :: forall f y x z w
     . ( ConjugatedLikelihood f y x z w, Generative Natural w
       , Generative Natural z, Map Natural f y x )
    => Int
    -> Natural # AffineHarmonium f y x z w -- ^ Categorical likelihood
    -> Random (Sample (z,w)) -- ^ Conjugation parameters
sampleConjugated n hrm = do
    let (lkl,nw) = split hrm
        nw' = nw + snd (conjugationParameters lkl)
    ws <- sample n nw'
    let mws :: [Mean # w]
        mws = sufficientStatistic <$> ws
    zs <- mapM samplePoint $ lkl >$+> mws
    return $ zip zs ws

-- | The conjugation parameters of a conjugated `Harmonium`.
conjugatedPotential
    :: ( LegendreExponentialFamily w, ConjugatedLikelihood f y x z w )
    => Natural # AffineHarmonium f y x z w -- ^ Categorical likelihood
    -> Double -- ^ Conjugation parameters
conjugatedPotential hrm = do
    let (lkl,nw) = split hrm
        (rho0,rprms) = conjugationParameters lkl
     in potential (nw + rprms) + rho0


--- Internal ---


-- Conjugation --

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableLogDensity
    :: forall f y x z w
    . ( ExponentialFamily z, ExponentialFamily y
      , LegendreExponentialFamily w, Translation w x, Translation z y
      , Map Natural f x y, Bilinear f y x )
    => Natural # AffineHarmonium f y x z w
    -> Sample z
    -> [Double]
unnormalizedHarmoniumObservableLogDensity hrm zs =
    let (pstr,nz) = split $ transposeHarmonium hrm
        mzs = sufficientStatistic <$> zs
        nrgs = zipWith (+) (dotMap nz mzs) $ potential <$> pstr >$+> mzs
     in zipWith (+) nrgs $ logBaseMeasure (Proxy @z) <$> zs

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedDensities
    :: forall f x y z w
    . ( Bilinear f y x, Translation z y
      , LegendreExponentialFamily z, ExponentialFamily y
      , LegendreExponentialFamily w, Translation w x, Map Natural f x y)
      => (Double, Natural # w) -- ^ Conjugation Parameters
      -> Natural # AffineHarmonium f y x z w
      -> Sample z
      -> [Double]
logConjugatedDensities (rho0,rprms) hrm z =
    let udns = unnormalizedHarmoniumObservableLogDensity hrm z
        nx = snd $ split hrm
     in subtract (potential (nx + rprms) + rho0) <$> udns

-- Mixtures --

mixtureLikelihoodConjugationParameters
    :: (KnownNat k, LegendreExponentialFamily z, Translation z y)
    => Natural # Affine Tensor y z (Categorical k) -- ^ Categorical likelihood
    -> (Double, Natural # Categorical k) -- ^ Conjugation parameters
mixtureLikelihoodConjugationParameters aff =
    let (nz,nyx) = split aff
        rho0 = potential nz
        rprms = S.map (\nyxi -> subtract rho0 . potential $ nz >+> nyxi) $ toColumns nyx
     in (rho0, Point rprms)

affineMixtureToMixture
    :: (KnownNat k, Manifold z, Manifold y, Translation z y)
    => Natural # AffineMixture y z k
    -> Natural # Mixture z k
affineMixtureToMixture lmxmdl =
    let (flsk,nk) = split lmxmdl
        (nls,nlk) = split flsk
        nlsk = fromColumns . S.map (0 >+>) $ toColumns nlk
     in join (join nls nlsk) nk

mixtureToAffineMixture
    :: (KnownNat k, Manifold y, Manifold z, Translation z y)
    => Mean # Mixture z k
    -> Mean # AffineMixture y z k
mixtureToAffineMixture mxmdl =
    let (flsk,mk) = split mxmdl
        (mls,mlsk) = split flsk
        mlk = fromColumns . S.map anchor $ toColumns mlsk
     in join (join mls mlk) mk


-- Linear Gaussian Harmoniums --


linearGaussianHarmoniumConjugationParameters
    :: (KnownNat n, KnownNat k)
    => Natural # Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)
    -> (Double, Natural # MultivariateNormal k) -- ^ Conjugation parameters
linearGaussianHarmoniumConjugationParameters aff =
    let (thts,tht30) = split aff
        (tht1,tht2) = splitNaturalMultivariateNormal thts
        tht3 = toMatrix tht30
        ttht3 = S.transpose tht3
        itht2 = S.pseudoInverse tht2
        rho0 = -0.25 * tht1 `S.dotProduct` (itht2 `S.matrixVectorMultiply` tht1)
            -0.5 * (log . S.determinant . negate $ 2*tht2)
        rho1 = -0.5 * ttht3 `S.matrixVectorMultiply` (itht2 `S.matrixVectorMultiply` tht1)
        rho2 = -0.25 * ttht3 `S.matrixMatrixMultiply` (itht2 `S.matrixMatrixMultiply` tht3)
     in (rho0, joinNaturalMultivariateNormal rho1 rho2)

univariateToLinearGaussianHarmonium
    :: c # AffineHarmonium Tensor NormalMean NormalMean Normal Normal
    -> c # LinearGaussianHarmonium 1 1
univariateToLinearGaussianHarmonium hrm =
    let (z,zx,x) = splitHarmonium hrm
     in joinHarmonium (breakPoint z) (breakPoint zx) (breakPoint x)

linearGaussianHarmoniumToUnivariate
    :: c # LinearGaussianHarmonium 1 1
    -> c # AffineHarmonium Tensor NormalMean NormalMean Normal Normal
linearGaussianHarmoniumToUnivariate hrm =
    let (z,zx,x) = splitHarmonium hrm
     in joinHarmonium (breakPoint z) (breakPoint zx) (breakPoint x)

univariateToLinearModel
    :: Natural # Affine Tensor NormalMean Normal NormalMean
    -> Natural # Affine Tensor (MVNMean 1) (MultivariateNormal 1) (MVNMean 1)
univariateToLinearModel aff =
    let (z,zx) = split aff
     in join (breakPoint z) (breakPoint zx)

naturalLinearGaussianHarmoniumToJoint
    :: (KnownNat n, KnownNat k)
    => Natural # LinearGaussianHarmonium n k
    -> Natural # MultivariateNormal (n+k)
naturalLinearGaussianHarmoniumToJoint hrm =
    let (z,zx,x) = splitHarmonium hrm
        zxmtx = toMatrix zx/2
        mvnz = splitNaturalMultivariateNormal z
        mvnx = splitNaturalMultivariateNormal x
        (mu,cvr) = fromLinearGaussianHarmonium0 mvnz zxmtx mvnx
     in joinNaturalMultivariateNormal mu cvr

naturalJointToLinearGaussianHarmonium
    :: (KnownNat n, KnownNat k)
    => Natural # MultivariateNormal (n+k)
    -> Natural # LinearGaussianHarmonium n k
naturalJointToLinearGaussianHarmonium mvn =
    let (mu,cvr) = splitNaturalMultivariateNormal mvn
        ((muz,cvrz),zxmtx,(mux,cvrx)) = toLinearGaussianHarmonium0 mu cvr
        zx = 2*fromMatrix zxmtx
        z = joinNaturalMultivariateNormal muz cvrz
        x = joinNaturalMultivariateNormal mux cvrx
     in joinHarmonium z zx x

meanLinearGaussianHarmoniumToJoint
    :: (KnownNat n, KnownNat k)
    => Mean # LinearGaussianHarmonium n k
    -> Mean # MultivariateNormal (n+k)
meanLinearGaussianHarmoniumToJoint hrm =
    let (z,zx,x) = splitHarmonium hrm
        zxmtx = toMatrix zx
        mvnz = splitMeanMultivariateNormal z
        mvnx = splitMeanMultivariateNormal x
        (mu,cvr) = fromLinearGaussianHarmonium0 mvnz zxmtx mvnx
     in joinMeanMultivariateNormal mu cvr

meanJointToLinearGaussianHarmonium
    :: (KnownNat n, KnownNat k)
    => Mean # MultivariateNormal (n+k)
    -> Mean # LinearGaussianHarmonium n k
meanJointToLinearGaussianHarmonium mvn =
    let (mu,cvr) = splitMeanMultivariateNormal mvn
        ((muz,cvrz),zxmtx,(mux,cvrx)) = toLinearGaussianHarmonium0 mu cvr
        zx = fromMatrix zxmtx
        z = joinMeanMultivariateNormal muz cvrz
        x = joinMeanMultivariateNormal mux cvrx
     in joinHarmonium z zx x

fromLinearGaussianHarmonium0
    :: (KnownNat n, KnownNat k)
    => (S.Vector n Double, S.Matrix n n Double)
    -> S.Matrix n k Double
    -> (S.Vector k Double, S.Matrix k k Double)
    -> (S.Vector (n+k) Double, S.Matrix (n+k) (n+k) Double)
fromLinearGaussianHarmonium0 (muz,cvrz) zxmtx (mux,cvrx) =
    let mu = muz S.++ mux
        top = S.horizontalConcat cvrz zxmtx
        btm = S.horizontalConcat (S.transpose zxmtx) cvrx
     in (mu, S.verticalConcat top btm)

toLinearGaussianHarmonium0
    :: (KnownNat n, KnownNat k)
    => S.Vector (n+k) Double
    -> S.Matrix (n+k) (n+k) Double
    -> ( (S.Vector n Double, S.Matrix n n Double)
       , S.Matrix n k Double
       , (S.Vector k Double, S.Matrix k k Double) )
toLinearGaussianHarmonium0 mu cvr =
    let (muz,mux) = S.splitAt mu
        (tops,btms) = S.splitAt $ S.toRows cvr
        (cvrzs,zxmtxs) = S.splitAt . S.toColumns $ S.fromRows tops
        cvrz = S.fromColumns cvrzs
        zxmtx = S.fromColumns zxmtxs
        cvrx = S.fromColumns . S.drop . S.toColumns $ S.fromRows btms
     in ((muz,cvrz),zxmtx,(mux,cvrx))

harmoniumLogBaseMeasure
    :: forall f y x z w . (ExponentialFamily z, ExponentialFamily w)
    => Proxy (AffineHarmonium f y x z w)
    -> SamplePoint (z,w)
    -> Double
harmoniumLogBaseMeasure _ (z,w) =
    logBaseMeasure (Proxy @z) z + logBaseMeasure (Proxy @w) w

isotropicGaussianHarmoniumToLinear
    :: (KnownNat n, KnownNat k)
    => Natural # IsotropicGaussianHarmonium n k
    -> Natural # LinearGaussianHarmonium n k
isotropicGaussianHarmoniumToLinear isohrm =
    let (lkl,prr) = split isohrm
        (iso,tns) = split lkl
        lkl' = join (isotropicNormalToFull iso) tns
     in join lkl' prr

linearGaussianHarmoniumToIsotropic
    :: (KnownNat n, KnownNat k)
    => Mean # LinearGaussianHarmonium n k
    -> Mean # IsotropicGaussianHarmonium n k
linearGaussianHarmoniumToIsotropic lnrhrm =
    let (lkl,prr) = split lnrhrm
        (lnr,tns) = split lkl
        lkl' = join (fullNormalToIsotropic lnr) tns
     in join lkl' prr



--- Instances ---


instance Manifold (AffineHarmonium f y x z w) => Statistical (AffineHarmonium f y x z w) where
    type SamplePoint (AffineHarmonium f y x z w) = SamplePoint (z,w)

instance ( ExponentialFamily z, ExponentialFamily x, Translation z y
         , Translation w x
         , ExponentialFamily w, ExponentialFamily y, Bilinear f y x )
  => ExponentialFamily (AffineHarmonium f y x z w) where
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

instance ( KnownNat k, LegendreExponentialFamily z
         , Translation y z, LegendreExponentialFamily y, SamplePoint z ~ SamplePoint y)
  => ConjugatedLikelihood Tensor z (Categorical k) y (Categorical k) where
    conjugationParameters = mixtureLikelihoodConjugationParameters

instance ConjugatedLikelihood Tensor NormalMean NormalMean Normal Normal where
    conjugationParameters aff =
        let rprms :: Natural # MultivariateNormal 1
            (rho0,rprms) = conjugationParameters $ univariateToLinearModel aff
         in (rho0,breakPoint rprms)

instance (KnownNat n, KnownNat k) => ConjugatedLikelihood Tensor (MVNMean n) (MVNMean k)
    (MultivariateNormal n) (MultivariateNormal k) where
        conjugationParameters = linearGaussianHarmoniumConjugationParameters

instance (KnownNat n, KnownNat k) => ConjugatedLikelihood Tensor (MVNMean n) (MVNMean k)
    (IsotropicNormal n) (MultivariateNormal k) where
        conjugationParameters aff =
            let (iso,tns) = split aff
                aff' = join (isotropicNormalToFull iso) tns
             in linearGaussianHarmoniumConjugationParameters aff'

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

instance ( KnownNat k, Manifold y, Manifold z
         , LegendreExponentialFamily z, Translation z y )
  => Transition Natural Mean (AffineMixture y z k) where
    transition mxmdl0 =
        let mxmdl = affineMixtureToMixture mxmdl0
            (nzs,nx) = splitNaturalMixture mxmdl
            mx = toMean nx
            mzs = S.map transition nzs
         in mixtureToAffineMixture $ joinMeanMixture mzs mx

instance ( KnownNat k, Manifold y, Manifold z, LegendreExponentialFamily z
         , Generative Natural z, Translation z y )
  => Generative Natural (AffineMixture y z k) where
      sample n = sampleConjugated n . affineMixtureToMixture

instance (KnownNat k, DuallyFlatExponentialFamily z)
  => Transition Mean Natural (Mixture z k) where
    transition mhrm =
        let (mzs,mx) = splitMeanMixture mhrm
            nx = transition mx
            nzs = S.map transition mzs
         in joinNaturalMixture nzs nx

instance (KnownNat k, LegendreExponentialFamily z, Transition Natural Source z)
  => Transition Natural Source (Mixture z k) where
    transition nhrm =
        let (nzs,nx) = splitNaturalMixture nhrm
            sx = transition nx
            szs = S.map transition nzs
         in joinSourceMixture szs sx

instance (KnownNat k, LegendreExponentialFamily z, Transition Source Natural z)
  => Transition Source Natural (Mixture z k) where
    transition shrm =
        let (szs,sx) = splitSourceMixture shrm
            nx = transition sx
            nzs = S.map transition szs
         in joinNaturalMixture nzs nx

instance Transition Natural Mean
  (AffineHarmonium Tensor NormalMean NormalMean Normal Normal) where
      transition = linearGaussianHarmoniumToUnivariate . transition . univariateToLinearGaussianHarmonium

instance Transition Mean Natural
  (AffineHarmonium Tensor NormalMean NormalMean Normal Normal) where
      transition =  linearGaussianHarmoniumToUnivariate . transition . univariateToLinearGaussianHarmonium

instance (KnownNat n, KnownNat k) => Transition Natural Mean (LinearGaussianHarmonium n k) where
      transition = meanJointToLinearGaussianHarmonium . transition
        . naturalLinearGaussianHarmoniumToJoint

instance (KnownNat n, KnownNat k) => Transition Natural Mean (IsotropicGaussianHarmonium n k) where
      transition = linearGaussianHarmoniumToIsotropic . transition . isotropicGaussianHarmoniumToLinear

instance (KnownNat n, KnownNat k) => Transition Mean Natural (LinearGaussianHarmonium n k) where
      transition = naturalJointToLinearGaussianHarmonium . transition
        . meanLinearGaussianHarmoniumToJoint

--type instance PotentialCoordinates (Mixture z k) = Natural
--
--instance (KnownNat k, LegendreExponentialFamily z) => Legendre (Mixture z k) where
--      potential = conjugatedPotential

type instance PotentialCoordinates (AffineHarmonium f y x z w) = Natural

instance ( Manifold (f y x), LegendreExponentialFamily w, ConjugatedLikelihood f y x z w )
  => Legendre (AffineHarmonium f y x z w) where
      potential = conjugatedPotential

instance ( Manifold (f y x), LegendreExponentialFamily w
         , Transition Mean Natural (AffineHarmonium f y x z w), ConjugatedLikelihood f y x z w )
  => DuallyFlat (AffineHarmonium f y x z w) where
    dualPotential mhrm =
        let nhrm = toNatural mhrm
         in mhrm <.> nhrm - potential nhrm

instance ( Bilinear f y x, ExponentialFamily y, ExponentialFamily x
         , LegendreExponentialFamily w, ConjugatedLikelihood f y x z w )
  => AbsolutelyContinuous Natural (AffineHarmonium f y x z w) where
    logDensities = exponentialFamilyLogDensities

instance ( ConjugatedLikelihood f y x z w, LegendreExponentialFamily z
         , ExponentialFamily y, LegendreExponentialFamily w
         , Map Natural f x y, Bilinear f x y )
  => ObservablyContinuous Natural (AffineHarmonium f y x z w) where
    logObservableDensities hrm zs =
        let rho0rprms = harmoniumConjugationParameters hrm
         in logConjugatedDensities rho0rprms hrm zs

instance ( LegendreExponentialFamily z, LegendreExponentialFamily w
         , ConjugatedLikelihood f y x z w, Map Natural f x y, Bilinear f x y
         , LegendreExponentialFamily (AffineHarmonium f y x z w)
         , Manifold (f y x), SamplePoint z ~ t, ExponentialFamily y)
  => LogLikelihood Natural (AffineHarmonium f y x z w) t where
    logLikelihood xs hrm =
         average $ logObservableDensities hrm xs
    logLikelihoodDifferential zs hrm =
        let pxs = expectationStep zs hrm
            qxs = transition hrm
         in pxs - qxs

instance ( Translation z y, Manifold w, Manifold (f y x) )
  => Translation (AffineHarmonium f y x z w) y where
      (>+>) hrm ny =
          let (nz,nyx,nw) = splitHarmonium hrm
           in joinHarmonium (nz >+> ny) nyx nw
      anchor hrm =
          let (nz,_,_) = splitHarmonium hrm
           in anchor nz


{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{-# HLINT ignore "Parenthesize unary negation" #-}

{- | An Exponential Family 'Harmonium' is a product exponential family with a
particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
A 'Mixture' model is a special case of harmonium.
-}
module Goal.Graphical.Models.Harmonium (
    -- * Harmoniums
    AffineHarmonium (AffineHarmonium),
    Harmonium,
    KnownHarmonium,
    KnownAffineHarmonium,

    -- ** Constuction
    splitHarmonium,
    joinHarmonium,

    -- ** Manipulation
    transposeHarmonium,

    -- ** Evaluation
    expectationStep,

    -- ** Sampling
    initialPass,
    gibbsPass,

    -- ** Conjugated Harmoniums
    ConjugatedLikelihood (conjugationParameters),
    logConjugatedDensities,
    harmoniumConjugationParameters,
    joinConjugatedHarmonium,
    splitConjugatedHarmonium,

    -- ** Bayesian Parameter Estimation
    ConjugatedPoissonGamma,

    -- ** Mixture Models
    Mixture,
    AffineMixture,
    joinNaturalMixture,
    splitNaturalMixture,
    joinMeanMixture,
    splitMeanMixture,
    joinSourceMixture,
    splitSourceMixture,

    -- ** Linear Gaussian Harmoniums
    LinearGaussianHarmonium,
    FullGaussianHarmonium,
    DiagonalGaussianHarmonium,
    IsotropicGaussianHarmonium,

    -- * Boltzmann-Gaussian Harmonium
    BoltzmannLinearModel,
    GaussianBoltzmannHarmonium,
) where

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models

import Goal.Core.Vector.Generic qualified as G
import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Proxy (Proxy (Proxy))

--- Types ---

-- | A 2-layer harmonium.
newtype AffineHarmonium f x0 z0 x z = AffineHarmonium (Affine f x0 x z0, z)

type Harmonium f x z = AffineHarmonium f x z x z

type instance Observation (AffineHarmonium f x0 z0 x z) = SamplePoint x

{- | A 'Mixture' model is simply a 'AffineHarmonium' where the latent variable is
'Categorical'.
-}
type Mixture x k = Harmonium L.Full x (Categorical k)

-- | A 'Mixture' where only a subset of the component parameters are mixed.
type AffineMixture x0 x k =
    AffineHarmonium L.Full x0 (Categorical k) x (Categorical k)

-- | A `MultivariateNormal` reintrepreted as a join distribution over two component `MultivariateNormal`s.
type GeneralizedGaussianHarmonium t n z0 z =
    AffineHarmonium L.Full (StandardNormal n) z0 (MultivariateNormal t n) z

{- | A `MultivariateNormal` reintrepreted as a join distribution over two component `MultivariateNormal`s.
type LinearGaussianHarmonium f n k =
    AffineHarmonium L.Full (StandardNormal n) (StandardNormal k) (MultivariateNormal f n) (FullNormal k)
-}
type LinearGaussianHarmonium t n k =
    GeneralizedGaussianHarmonium t n (StandardNormal k) (FullNormal k)

-- | A `LinearGaussianHarmonium` with all covariances.
type FullGaussianHarmonium n k = LinearGaussianHarmonium L.PositiveDefinite n k

-- | A `LinearGaussianHarmonium` with a diagonal covariance between oservable variables.
type IsotropicGaussianHarmonium n k = LinearGaussianHarmonium L.Scale n k

-- | A `LinearGaussianHarmonium` with a scale covariance between oservable variables.
type DiagonalGaussianHarmonium n k = LinearGaussianHarmonium L.Diagonal n k

-- | Basic typeclass synonym for Harmoniums.
type KnownAffineHarmonium f x0 z0 x z =
    ( KnownLinear f x0 z0
    , KnownLinear f z0 x0
    , LinearSubspace x x0
    , LinearSubspace z z0
    , ExponentialFamily x0
    , ExponentialFamily z0
    , ExponentialFamily x
    , ExponentialFamily z
    )

type KnownHarmonium f x z = KnownAffineHarmonium f x z x z

--- Boltzmann-Gaussian Harmonium

-- | Linear models are linear functions with additive Guassian noise.
type BoltzmannLinearModel f n k = Affine L.Full (StandardNormal n) (MultivariateNormal f n) (Replicated k Bernoulli)

type LinearModel f n z = Affine L.Full (StandardNormal n) (MultivariateNormal f n) z

type GaussianBoltzmannHarmonium t n k =
    AffineHarmonium L.Full (StandardNormal n) (Replicated k Bernoulli) (MultivariateNormal t n) (Boltzmann k)

--- Bayesian Parameter Estimation ---

type ConjugatedPoissonGamma = AffineHarmonium L.Identity Poisson GammaRate Poisson Gamma

--- Classes ---

-- | The conjugation parameters of a conjugated likelihood.
class
    (KnownAffineHarmonium f x0 z0 x z) =>
    ConjugatedLikelihood f x0 z0 x z
    where
    conjugationParameters ::
        -- | Categorical likelihood
        Natural # Affine f x0 x z0 ->
        -- | Conjugation parameters
        (Double, Natural # z)

--- Functions ---

-- Construction --

-- | Creates a 'Harmonium' from component parameters.
joinHarmonium ::
    (KnownLinear f x0 z0, Manifold x, Manifold z) =>
    -- | Visible layer biases
    c # x ->
    -- | ^ Interaction parameters
    c # Linear f x0 z0 ->
    -- | Hidden layer Biases
    c # z ->
    -- | Harmonium
    c # AffineHarmonium f x0 z0 x z
joinHarmonium nx nx0z0 = join (join nx nx0z0)

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium ::
    (KnownLinear f x0 z0, Manifold x, Manifold z) =>
    -- | Harmonium
    c # AffineHarmonium f x0 z0 x z ->
    -- | Biases and interaction parameters
    (c # x, c # Linear f x0 z0, c # z)
splitHarmonium hrm =
    let (fxz0, nz) = split hrm
        (nx, nx0z0) = split fxz0
     in (nx, nx0z0, nz)

-- | Build a mixture model in source coordinates.
joinSourceMixture ::
    (KnownNat k, Manifold x) =>
    -- | Mixture components
    S.Vector (k + 1) (Source # x) ->
    -- | Weights
    Source # Categorical k ->
    Source # Mixture x k
joinSourceMixture sxs sk =
    let (sx, sxs') = S.splitAt sxs
        aff = join (S.head sx) (fromColumns sxs')
     in join aff sk

-- | Build a mixture model in source coordinates.
splitSourceMixture ::
    (KnownNat k, Manifold x) =>
    Source # Mixture x k ->
    (S.Vector (k + 1) (Source # x), Source # Categorical k)
splitSourceMixture mxmdl =
    let (aff, sk) = split mxmdl
        (sx0, sxs0') = split aff
     in (S.cons sx0 $ toColumns sxs0', sk)

-- | Build a mixture model in mean coordinates.
joinMeanMixture ::
    (KnownNat k, Manifold x) =>
    -- | Mixture components
    S.Vector (k + 1) (Mean # x) ->
    -- | Weights
    Mean # Categorical k ->
    Mean # Mixture x k
joinMeanMixture mxs mk =
    let wghts = categoricalWeights mk
        wmxs = S.zipWith (.>) wghts mxs
        mx = S.sum wmxs
        twmxs = S.tail wmxs
        mxk = transpose . fromRows $ twmxs
     in joinHarmonium mx mxk mk

-- | Split a mixture model in mean coordinates.
splitMeanMixture ::
    (KnownNat k, LegendreExponentialFamily x) =>
    Mean # Mixture x k ->
    (S.Vector (k + 1) (Mean # x), Mean # Categorical k)
splitMeanMixture hrm =
    let (mx, mxz, mk) = splitHarmonium hrm
        twmxs = toRows $ transpose mxz
        wmxs = S.cons (mx - S.sum twmxs) twmxs
        wghts = categoricalWeights mk
        mxs = S.zipWith (/>) wghts wmxs
     in (mxs, mk)

-- | A convenience function for building a categorical harmonium/mixture model.
joinNaturalMixture ::
    forall k x.
    ( ConjugatedLikelihood L.Full x (Categorical k) x (Categorical k)
    , KnownNat k
    , LegendreExponentialFamily x
    ) =>
    -- | Mixture components
    S.Vector (k + 1) (Natural # x) ->
    -- | Weights
    Natural # Categorical k ->
    -- | Mixture Model
    Natural # Mixture x k
joinNaturalMixture nxs0 nk0 =
    let nx0 :: S.Vector 1 (Natural # x)
        (nx0, nxs0') = S.splitAt nxs0
        nx = S.head nx0
        nxs = S.map (subtract nx) nxs0'
        nxk = fromColumns nxs
        affxk = join nx nxk
        rprms = snd $ conjugationParameters affxk
        nk = nk0 - rprms
     in joinHarmonium nx nxk nk

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitNaturalMixture ::
    forall k x.
    ( ConjugatedLikelihood L.Full x (Categorical k) x (Categorical k)
    , KnownNat k
    , LegendreExponentialFamily x
    ) =>
    -- | Categorical harmonium
    Natural # Mixture x k ->
    -- | (components, weights)
    (S.Vector (k + 1) (Natural # x), Natural # Categorical k)
splitNaturalMixture hrm =
    let (nx, nxk, nk) = splitHarmonium hrm
        affxk = join nx nxk
        rprms = snd $ conjugationParameters affxk
        nk0 = nk + rprms
        nxs = toColumns nxk
        nxs0' = S.map (+ nx) nxs
     in (S.cons nx nxs0', nk0)

-- Manipulation --

-- | Swap the biases and 'transpose' the interaction parameters of the given 'Harmonium'.
transposeHarmonium ::
    (KnownLinear f x0 z0, KnownLinear f z0 x0, Manifold x, Manifold z) =>
    c # AffineHarmonium f x0 z0 x z ->
    c # AffineHarmonium f z0 x0 z x
transposeHarmonium hrm =
    let (nz, nyx, nw) = splitHarmonium hrm
     in joinHarmonium nw (transpose nyx) nz

-- Evaluation --

{- | Computes the joint expectations of a harmonium based on a sample from the
observable layer.
-}
expectationStep ::
    (KnownAffineHarmonium f x0 z0 x z, LegendreExponentialFamily z) =>
    -- | Model Samples
    Sample x ->
    -- | Harmonium
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Harmonium expected sufficient statistics
    Mean # AffineHarmonium f x0 z0 x z
{-# INLINE expectationStep #-}
expectationStep xs hrm =
    let mxs = sufficientStatistic <$> xs
        mx0s = linearProjection <$> mxs
        pstr = fst . split $ transposeHarmonium hrm
        mzs = transition <$> pstr >$> mx0s
        mz0s = linearProjection <$> mzs
        mx0z0 = (>$<) mx0s mz0s
     in joinHarmonium (average mxs) mx0z0 $ average mzs

harmoniumLogBaseMeasure ::
    forall f z0 x0 z x.
    (ExponentialFamily z, ExponentialFamily x) =>
    Proxy (AffineHarmonium f z0 x0 z x) ->
    SamplePoint (z, x) ->
    Double
{-# INLINE harmoniumLogBaseMeasure #-}
harmoniumLogBaseMeasure _ (z, x) =
    logBaseMeasure (Proxy @z) z + logBaseMeasure (Proxy @x) x

---- Sampling --

-- | Initialize a Gibbs chain from a set of observations.
initialPass ::
    forall f x0 z0 x z.
    (KnownAffineHarmonium f x0 z0 x z, Generative Natural z, LegendreExponentialFamily z) =>
    -- | Harmonium
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Model Samples
    Sample x ->
    Random (Sample (x, z))
initialPass hrm xs = do
    let pstr = fst . split $ transposeHarmonium hrm
        mxs :: [Mean # x]
        mxs = sufficientStatistic <$> xs
        mx0s = linearProjection <$> mxs
    zs <- mapM samplePoint $ pstr >$> mx0s
    return $ zip xs zs

-- | Update a 'Sample' with Gibbs sampling.
gibbsPass ::
    forall f x0 z0 x z.
    ( KnownAffineHarmonium f x0 z0 x z
    , Generative Natural z
    , Generative Natural x
    ) =>
    -- | Harmonium
    Natural # AffineHarmonium f x0 z0 x z ->
    Sample (x, z) ->
    Random (Sample (x, z))
gibbsPass hrm xzs = do
    let zs = snd <$> xzs
        mzs :: [Mean # z]
        mzs = sufficientStatistic <$> zs
        mz0s = linearProjection <$> mzs
        pstr = fst . split $ transposeHarmonium hrm
        lkl = fst $ split hrm
    xs' <- mapM samplePoint $ lkl >$> mz0s
    let mxs' :: [Mean # x]
        mxs' = sufficientStatistic <$> xs'
        mx0s' = linearProjection <$> mxs'
    zs' <- mapM samplePoint $ pstr >$> mx0s'
    return $ zip xs' zs'

-- Conjugation --

-- | The conjugation parameters of a conjugated `Harmonium`.
harmoniumConjugationParameters ::
    (ConjugatedLikelihood f x0 z0 x z) =>
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Conjugation parameters
    (Double, Natural # z)
{-# INLINE harmoniumConjugationParameters #-}
harmoniumConjugationParameters hrm =
    conjugationParameters . fst $ split hrm

-- | The conjugation parameters of a conjugated `Harmonium`.
splitConjugatedHarmonium ::
    (ConjugatedLikelihood f x0 z0 x z) =>
    Natural # AffineHarmonium f x0 z0 x z ->
    (Natural # Affine f x0 x z0, Natural # z)
{-# INLINE splitConjugatedHarmonium #-}
splitConjugatedHarmonium hrm =
    let (lkl, nw) = split hrm
        cw = snd $ conjugationParameters lkl
     in (lkl, nw + cw)

-- | The conjugation parameters of a conjugated `Harmonium`.
joinConjugatedHarmonium ::
    (ConjugatedLikelihood f x0 z0 x z) =>
    -- | Conjugation parameters
    Natural # Affine f x0 x z0 ->
    Natural # z ->
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z
{-# INLINE joinConjugatedHarmonium #-}
joinConjugatedHarmonium lkl nw =
    let cw = snd $ conjugationParameters lkl
     in join lkl $ nw - cw

-- | The conjugation parameters of a conjugated `Harmonium`.
sampleConjugated ::
    forall f x0 z0 x z.
    ( ConjugatedLikelihood f x0 z0 x z
    , Generative Natural x
    , Generative Natural z
    ) =>
    Int ->
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Conjugation parameters
    Random (Sample (x, z))
{-# INLINE sampleConjugated #-}
sampleConjugated n hrm = do
    let (lkl, nz) = splitConjugatedHarmonium hrm
    zs <- sample n nz
    let mzs :: [Mean # z]
        mzs = sufficientStatistic <$> zs
    xs <- mapM samplePoint $ lkl >$+> mzs
    return $ zip xs zs

-- | The conjugation parameters of a conjugated `Harmonium`.
conjugatedPotential ::
    (LegendreExponentialFamily z, ConjugatedLikelihood f x0 z0 x z) =>
    -- | Categorical likelihood
    Natural # AffineHarmonium f x0 z0 x z ->
    -- | Conjugation parameters
    Double
{-# INLINE conjugatedPotential #-}
conjugatedPotential hrm = do
    let (lkl, nz) = split hrm
        (rho0, rprms) = conjugationParameters lkl
     in potential (nz + rprms) + rho0

--- Internal ---

-- Conjugation --

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableLogDensities ::
    forall f x0 z0 x z.
    ( LegendreExponentialFamily z
    , KnownAffineHarmonium f x0 z0 x z
    ) =>
    Natural # AffineHarmonium f x0 z0 x z ->
    Sample x ->
    [Double]
{-# INLINE unnormalizedHarmoniumObservableLogDensities #-}
unnormalizedHarmoniumObservableLogDensities hrm xs =
    let (pstr, nx) = split $ transposeHarmonium hrm
        mxs = sufficientStatistic <$> xs
        nrgs = zipWith (+) (dotMap nx mxs) $ potential <$> pstr >$+> mxs
     in zipWith (+) nrgs $ logBaseMeasure (Proxy @x) <$> xs

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedDensities ::
    ( LegendreExponentialFamily z
    , KnownAffineHarmonium f x0 z0 x z
    ) =>
    -- | Conjugation Parameters
    (Double, Natural # z) ->
    Natural # AffineHarmonium f x0 z0 x z ->
    Sample x ->
    [Double]
{-# INLINE logConjugatedDensities #-}
logConjugatedDensities (rho0, rprms) hrm xs =
    let udns = unnormalizedHarmoniumObservableLogDensities hrm xs
        nz = snd $ split hrm
     in subtract (potential (nz + rprms) + rho0) <$> udns

-- Mixtures --

mixtureLikelihoodConjugationParameters ::
    (KnownNat k, LegendreExponentialFamily x, LinearSubspace x x0) =>
    -- | Categorical likelihood
    Natural # Affine L.Full x0 x (Categorical k) ->
    -- | Conjugation parameters
    (Double, Natural # Categorical k)
{-# INLINE mixtureLikelihoodConjugationParameters #-}
mixtureLikelihoodConjugationParameters aff =
    let (nx, nx0z0) = split aff
        rho0 = potential nx
        rprms = S.map (\nx0z0i -> subtract rho0 . potential $ nx >+> nx0z0i) $ toColumns nx0z0
     in (rho0, Point rprms)

affineMixtureToMixture ::
    (KnownNat k, Manifold x0, Manifold x, LinearSubspace x x0) =>
    Natural # AffineMixture x0 x k ->
    Natural # Mixture x k
affineMixtureToMixture lmxmdl =
    let (flsk, nk) = split lmxmdl
        (nls, nlk) = split flsk
        nlsk = fromColumns . S.map (0 >+>) $ toColumns nlk
     in join (join nls nlsk) nk

mixtureToAffineMixture ::
    (KnownNat k, Manifold x, Manifold x0, LinearSubspace x x0) =>
    Mean # Mixture x k ->
    Mean # AffineMixture x0 x k
mixtureToAffineMixture mxmdl =
    let (flsk, mk) = split mxmdl
        (mls, mlsk) = split flsk
        mlk = fromColumns . S.map linearProjection $ toColumns mlsk
     in join (join mls mlk) mk

-- Gaussian Harmoniums --

linearModelConjugationParameters ::
    forall n f x s.
    (KnownCovariance f n, GeneralizedGaussian Natural x s) =>
    Natural # LinearModel f n x ->
    (Double, Natural # MomentParameters x s)
{-# INLINE linearModelConjugationParameters #-}
linearModelConjugationParameters aff =
    let (thts, tht3) = split aff
        tht2 :: Natural # CovarianceMatrix f n
        (tht1, tht2) = splitGaussian thts
        (itht20, lndt, _) = inverseLogDeterminant . negate $ 2 .> tht2
        itht2 = -2 .> itht20
        tht21 = itht2 >.> tht1
        chi = -0.25 * (tht1 <.> tht21) - 0.5 * lndt
        rho1 = -0.5 .> (transpose tht3 >.> tht21)
        rho2 :: Natural # Linear (CovarianceRep s) x x
        rho2 = fromTensor $ -0.25 .> changeOfBasis tht3 itht2
        rho = joinGaussian rho1 rho2
     in (chi, rho)

--- Instances ---

--- Deriving ---

deriving instance
    (Manifold (Affine f x0 x z0), Manifold z) =>
    Manifold (AffineHarmonium f x0 z0 x z)
deriving instance
    (Manifold (Affine f x0 x z0), Manifold z) =>
    Product (AffineHarmonium f x0 z0 x z)

--- Harmonium ---

instance (Manifold (AffineHarmonium f z0 x0 z x)) => Statistical (AffineHarmonium f z0 x0 z x) where
    type SamplePoint (AffineHarmonium f z0 x0 z x) = SamplePoint (z, x)

type instance PotentialCoordinates (AffineHarmonium f z0 x0 z x) = Natural

instance
    (KnownAffineHarmonium f x0 z0 x z) =>
    ExponentialFamily (AffineHarmonium f x0 z0 x z)
    where
    sufficientStatistic (z, w) =
        let mz = sufficientStatistic z
            mw = sufficientStatistic w
            my = linearProjection mz
            mx = linearProjection mw
         in joinHarmonium mz (my >.< mx) mw
    averageSufficientStatistic zws =
        let (zs, ws) = unzip zws
            mzs = sufficientStatistic <$> zs
            mws = sufficientStatistic <$> ws
            mys = linearProjection <$> mzs
            mxs = linearProjection <$> mws
         in joinHarmonium (average mzs) (mys >$< mxs) (average mws)
    logBaseMeasure = harmoniumLogBaseMeasure

instance
    ( ConjugatedLikelihood f x0 z0 x z
    , Generative Natural x
    , Generative Natural z
    ) =>
    Generative Natural (AffineHarmonium f x0 z0 x z)
    where
    sample = sampleConjugated

instance
    (LegendreExponentialFamily x, ConjugatedLikelihood f z0 x0 z x) =>
    Legendre (AffineHarmonium f z0 x0 z x)
    where
    {-# INLINE potential #-}
    potential = conjugatedPotential

instance
    ( LegendreExponentialFamily x
    , Transition Mean Natural (AffineHarmonium f z0 x0 z x)
    , ConjugatedLikelihood f z0 x0 z x
    ) =>
    DuallyFlat (AffineHarmonium f z0 x0 z x)
    where
    {-# INLINE dualPotential #-}
    dualPotential mhrm =
        let nhrm = toNatural mhrm
         in mhrm <.> nhrm - potential nhrm

instance
    (KnownAffineHarmonium f x0 z0 x z) =>
    LinearSubspace (AffineHarmonium f x0 z0 x z) x
    where
    (>+>) hrm nx' =
        let (nx, nxz, nz) = splitHarmonium hrm
         in joinHarmonium (nx + nx') nxz nz
    linearProjection hrm =
        let (nx, _, _) = splitHarmonium hrm
         in nx

instance
    (KnownAffineHarmonium f x0 z0 x z) =>
    LinearSubspace (AffineHarmonium f x0 z0 x z) z
    where
    (>+>) hrm nz' =
        let (nx, nxz, nz) = splitHarmonium hrm
         in joinHarmonium nx nxz (nz + nz')
    linearProjection hrm =
        let (_, _, nz) = splitHarmonium hrm
         in nz

instance
    ( LegendreExponentialFamily z
    , Legendre (AffineHarmonium f x0 z0 x z)
    , ExponentialFamily (AffineHarmonium f x0 z0 x z)
    ) =>
    AbsolutelyContinuous Natural (AffineHarmonium f x0 z0 x z)
    where
    logDensities = exponentialFamilyLogDensities

instance
    ( LegendreExponentialFamily z
    , ConjugatedLikelihood f x0 z0 x z
    ) =>
    ObservablyContinuous Natural (AffineHarmonium f x0 z0 x z)
    where
    logObservableDensities hrm zs =
        let rho0rprms = harmoniumConjugationParameters hrm
         in logConjugatedDensities rho0rprms hrm zs

instance
    ( LegendreExponentialFamily z
    , ExponentialFamily x
    , SamplePoint x ~ t
    , ConjugatedLikelihood f x0 z0 x z
    , Transition Natural Mean (AffineHarmonium f x0 z0 x z)
    ) =>
    LogLikelihood Natural (AffineHarmonium f x0 z0 x z) t
    where
    logLikelihood xs hrm =
        average $ logObservableDensities hrm xs
    logLikelihoodDifferential zs hrm =
        let pxs = expectationStep zs hrm
            qxs = transition hrm
         in pxs - qxs

--- Mixture ---

instance
    ( KnownAffineHarmonium L.Full x0 (Categorical k) x (Categorical k)
    , LegendreExponentialFamily x
    ) =>
    ConjugatedLikelihood L.Full x0 (Categorical k) x (Categorical k)
    where
    conjugationParameters = mixtureLikelihoodConjugationParameters

instance
    ( KnownNat k
    , Manifold y
    , Manifold z
    , LegendreExponentialFamily z
    , LinearSubspace z y
    ) =>
    Transition Natural Mean (AffineMixture y z k)
    where
    transition mxmdl0 =
        let mxmdl = affineMixtureToMixture mxmdl0
            (nzs, nx) = splitNaturalMixture mxmdl
            mx = toMean nx
            mzs = S.map transition nzs
         in mixtureToAffineMixture $ joinMeanMixture mzs mx

instance
    (KnownNat k, DuallyFlatExponentialFamily x) =>
    Transition Mean Natural (Mixture x k)
    where
    transition mhrm =
        let (mxs, mk) = splitMeanMixture mhrm
            nk = transition mk
            nxs = S.map transition mxs
         in joinNaturalMixture nxs nk

instance
    (KnownNat k, LegendreExponentialFamily x, Transition Mean Source x) =>
    Transition Mean Source (Mixture x k)
    where
    transition mhrm =
        let (mxs, mk) = splitMeanMixture mhrm
            nk = transition mk
            nxs = S.map transition mxs
         in joinSourceMixture nxs nk

instance
    (KnownNat k, LegendreExponentialFamily x, Transition Natural Source x) =>
    Transition Natural Source (Mixture x k)
    where
    transition nhrm =
        let (nxs, nk) = splitNaturalMixture nhrm
            sk = transition nk
            sxs = S.map transition nxs
         in joinSourceMixture sxs sk

instance
    (KnownNat k, LegendreExponentialFamily x, Transition Source Natural x) =>
    Transition Source Natural (Mixture x k)
    where
    transition shrm =
        let (sxs, sk) = splitSourceMixture shrm
            nk = transition sk
            nxs = S.map transition sxs
         in joinNaturalMixture nxs nk

--- Guassian Harmoniums ---

instance
    (KnownCovariance f n, GeneralizedGaussian Natural x s) =>
    ConjugatedLikelihood
        L.Full
        (StandardNormal n)
        x
        (MultivariateNormal f n)
        (MomentParameters x s)
    where
    conjugationParameters = linearModelConjugationParameters

instance
    (KnownCovariance f n, KnownNat k) =>
    Transition Natural Mean (LinearGaussianHarmonium f n k)
    where
    transition nlgh =
        let (nx, nsgmaxz, nz) = splitHarmonium nlgh
            (nmuz, nsgmaz) = splitGaussian nz
            (nmux, nsgmax) = splitGaussian nx
            prsx = -nsgmax * 2
            prsz = -nsgmaz * 2
            prsxz = -nsgmaxz
            prsxinv = inverse prsx
            schrz = schurComplement prsxinv prsxz (transpose prsxz) prsz
            sgmaz = inverse schrz
            sgmaxz = negate $ dualComposition prsxinv prsxz sgmaz
            muz = sgmaz >.> nmuz + nmux <.< sgmaxz
            mux0 = prsxinv >.> nmux
            mux = (mux0 - sgmaxz >.> (mux0 <.< prsxz)) + sgmaxz >.> nmuz
            sgmax = schurComplement (negate schrz) (transpose sgmaxz) sgmaxz prsxinv
            msgmaz = sgmaz + muz >.< muz
            msgmax = sgmax + mux >.< mux
            msgmaxz = sgmaxz + mux >.< muz
            mx = join mux msgmax
            mz = join muz msgmaz
         in joinHarmonium mx msgmaxz mz

-- instance
--     (KnownCovariance f n, KnownNat k) =>
--     Transition Mean Natural (LinearGaussianHarmonium f n k)
--     where
--     transition mlgh =
--         let (mx, msgmaxz, mz) = splitHarmonium mlgh
--             (mux, msgmax) = split mx
--             (muz, msgmaz) = split mz
--             sgmaz = msgmaz - muz >.< muz
--             sgmax = msgmax - mux >.< mux
--             sgmaxz = msgmaxz - mux >.< muz
--             sgmazinv = inverse sgmaz
--             nvrx0 = inverse $ schurComplement sgmazinv (transpose sgmaxz) sgmaxz sgmax
--             nvrx = -0.5 .> nvrx0
--             nxz = dualComposition nvrx0 sgmaxz sgmazinv
--             nmux = nvrx0 >.> mux - nxz >.> muz
--             nmvn = joinGaussian nmux nvrx
--             nlkl = join nmvn nxz
--             nz = toNatural $ join muz msgmaz
--          in joinConjugatedHarmonium (breakChart nlkl) nz
--
instance
    (KnownCovariance f n, KnownNat k) =>
    Transition Mean Source (LinearGaussianHarmonium f n k)
    where
    transition mlgh =
        let (mfxz, mz) = split mlgh
            (mx, mvrxz) = split mfxz
            mmux = fst $ split mx
            mmuz = fst $ split mz
            svrxz = breakChart $ mvrxz - mmux >.< mmuz
            sx = toSource mx
            sz = toSource mz
            sfxz = join sx svrxz
         in join sfxz sz

instance
    (KnownCovariance f n, KnownNat k) =>
    Transition Source Mean (LinearGaussianHarmonium f n k)
    where
    transition slgh =
        let (sfxz, sz) = split slgh
            (sx, svrxz) = split sfxz
            smux = fst $ split sx
            smuz = fst $ split sz
            mvrxz = breakChart $ svrxz + smux >.< smuz
            mx = toMean sx
            mz = toMean sz
            mfxz = join mx mvrxz
         in join mfxz mz

--- Gaussian-Boltzmann

instance
    (KnownCovariance t n, KnownBoltzmann L.Symmetric k) =>
    Transition Natural Mean (GaussianBoltzmannHarmonium t n k)
    where
    transition gbhrm =
        let (lmdl, bltz) = splitConjugatedHarmonium gbhrm
            blss = pointSampleSpace bltz
            mblss = sufficientStatistic <$> blss
            sblss = coordinates <$> mblss
            prbs = densities bltz blss
            bltzmtx = S.lowerTriangular . S.weightedAverageOuterProduct $ zip3 prbs sblss sblss
            (bltzmu, bltzcvr) = S.triangularSplitDiagonal bltzmtx
            mbltz = join (Point bltzmu) (Point bltzcvr)
            mxs = toMean <$> lmdl >$> mblss
            sxmus = coordinates . fst . split <$> mxs
            mxz = Point . G.toVector . S.weightedAverageOuterProduct $ zip3 prbs sxmus sblss
            mx = weightedAverage $ zip prbs mxs
         in joinHarmonium mx mxz mbltz

--- Generalized Gaussian Harmonium ---

--
-- I THINK THIS MIGHT BE WRONG. Why don't I need splitConjugated?
--
instance
    ( KnownCovariance t n
    , GeneralizedGaussian Mean z s
    , GeneralizedGaussian Natural z s
    , Transition Mean Natural (MomentParameters z s)
    ) =>
    Transition Mean Natural (GeneralizedGaussianHarmonium t n z (MomentParameters z s))
    where
    transition mlgh =
        let (mx, msgmaxz, mz) = splitHarmonium mlgh
            (mux, msgmax) = splitGaussian mx
            (muz, msgmaz) = splitGaussian mz
            sgmaz = msgmaz - muz >.< muz
            sgmax = msgmax - mux >.< mux
            sgmaxz = msgmaxz - mux >.< muz
            sgmazinv = inverse sgmaz
            nvrx0 = inverse $ schurComplement sgmazinv (transpose sgmaxz) sgmaxz sgmax
            nvrx = -0.5 .> nvrx0
            nxz = dualComposition nvrx0 sgmaxz sgmazinv
            nmux = nvrx0 >.> mux - nxz >.> muz
            nmvn = joinGaussian nmux nvrx
            nlkl = join nmvn nxz
            nz = toNatural mz
         in joinConjugatedHarmonium nlkl nz

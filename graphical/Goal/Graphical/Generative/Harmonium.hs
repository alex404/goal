{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE TypeApplications,UndecidableInstances #-}
-- | An Exponential Family 'Harmonium' is a product exponential family with a
-- particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
-- A 'Mixture' model is a special case of harmonium.
module Goal.Graphical.Generative.Harmonium
    (
    -- * Harmoniums
      Harmonium
    -- ** Constuction
    , splitHarmonium
    , splitBottomHarmonium
    , joinHarmonium
    , joinBottomHarmonium
    -- ** Manipulation
    , transposeHarmonium
    -- ** Evaluation
    , unnormalizedHarmoniumObservableDensity
    , unnormalizedHarmoniumObservableLogDensity
    -- ** Mixture Models
    , Mixture
    , joinNaturalMixture
    , splitNaturalMixture
    , joinMeanMixture
    , splitMeanMixture
    -- * Conjugated Harmoniums
    , conjugatedDensity
    , logConjugatedDensity
    , conjugatedDensities
    , logConjugatedDensities
    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Conditional
import Goal.Graphical.Generative

import qualified Goal.Core.Vector.Storable as S
--import qualified Data.Vector.Storable as VS


--- Types ---


-- | A 2-layer harmonium.
data Harmonium (f :: Type -> Type -> Type) z x

-- | A 'Mixture' model is simply a 'Harmonium' where the latent variable is
-- 'Categorical'.
type Mixture z k = Harmonium Tensor z (Categorical k)


--- Classes ---


-- | The conjugation parameters of a conjugated likelihood.
class (Manifold z, Manifold x, Manifold (f z x), Map Natural (Affine f) z x)
  => ConjugatedLikelihood f z x where
    conjugationParameters
        :: Natural # Affine f z x -- ^ Categorical likelihood
        -> (Double, Natural # x) -- ^ Conjugation parameters


--- Functions ---


-- Construction --

-- | Creates a 'Harmonium' from component parameters.
joinBottomHarmonium
    :: (Manifold z, Manifold (f z x), Manifold x)
    => c # Affine f z x -- ^ ^ Interaction parameters
    -> c # x -- ^ Hidden layer Biases
    -> c # Harmonium f z x -- ^ Harmonium
joinBottomHarmonium affzx nx =
     Point $ coordinates affzx S.++ coordinates nx

-- | Creates a 'Harmonium' from component parameters.
joinHarmonium
    :: (Manifold z, Manifold (f z x), Manifold x)
    => c # z -- ^ Visible layer biases
    -> c # f z x -- ^ ^ Interaction parameters
    -> c # x -- ^ Hidden layer Biases
    -> c # Harmonium f z x -- ^ Harmonium
joinHarmonium nz nzx nx =
     Point $ coordinates nz S.++ coordinates nzx S.++ coordinates nx

-- | Splits a 'Harmonium' into component parameters.
splitBottomHarmonium
    :: (Manifold z, Manifold (f z x), Manifold x)
    => c # Harmonium f z x -- ^ Harmonium
    -> (c # Affine f z x, c # x) -- ^ Biases and interaction parameters
splitBottomHarmonium hrm =
    let (caffzx,cx) = S.splitAt $ coordinates hrm
     in (Point caffzx, Point cx)

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium
    :: (Manifold z, Manifold (f z x), Manifold x)
    => c # Harmonium f z x -- ^ Harmonium
    -> (c # z, c # f z x, c # x) -- ^ Biases and interaction parameters
splitHarmonium hrm =
    let (cz,cs') = S.splitAt $ coordinates hrm
        (czx,cx) = S.splitAt cs'
     in (Point cz, Point czx, Point cx)

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
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx = nx0 - rprms
     in joinHarmonium nz nzx nx

-- | A convenience function for deconstructing a categorical harmonium/mixture model.
splitNaturalMixture
    :: forall k z . ( KnownNat k, LegendreExponentialFamily z )
    => Natural # Mixture z k -- ^ Categorical harmonium
    -> (S.Vector (k+1) (Natural # z), Natural # Categorical k) -- ^ (components, weights)
splitNaturalMixture hrm =
    let (nz,nzx,nx) = splitHarmonium hrm
        affzx = joinAffine nz nzx
        rprms = snd $ mixtureLikelihoodConjugationParameters affzx
        nx0 = nx + rprms
        nzs = S.map Point . S.toColumns $ toMatrix nzx
        nzs0' = S.map (+ nz) nzs
     in (S.cons nz nzs0',nx0)

-- Manipulation --

-- | Swap the biases and 'transpose' the interaction parameters of the given 'Harmonium'.
transposeHarmonium
    :: Bilinear f z x
    => Natural # Harmonium f z x
    -> Natural # Harmonium f x z
transposeHarmonium dhrm =
        let (nz,nzx,nx) = splitHarmonium dhrm
         in joinHarmonium nx (transpose nzx) nz

-- Evaluation --

---- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableDensity
    :: forall f z x . ( ExponentialFamily z, LegendreExponentialFamily x
                      , Map Natural f x z, Bilinear f z x )
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableDensity hrm z =
    let (nz,nzx,nx) = splitHarmonium hrm
        mz = sufficientStatistic z
     in exp (nz <.> mz + potential (nx + mz <.< nzx) + logBaseMeasure (Proxy @ z) z)

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableLogDensity
    :: forall f z x . ( ExponentialFamily z, LegendreExponentialFamily x
                      , Map Natural f x z, Bilinear f z x)
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableLogDensity hrm z =
    let (nz,nzx,nx) = splitHarmonium hrm
        mz = sufficientStatistic z
     in nz <.> mz + potential (nx + mz <.< nzx) + logBaseMeasure (Proxy @ z) z

-- | Computes the joint expectations of a harmonium based on a sample from the
-- observable layer.
harmoniumExpectationStep
    :: ( ExponentialFamily z, Map Natural f x z
       , Bilinear f z x, LegendreExponentialFamily x )
    => Sample z -- ^ Model Samples
    -> Natural # Harmonium f z x -- ^ Harmonium
    -> Mean # Harmonium f z x -- ^ Harmonium expected sufficient statistics
harmoniumExpectationStep zs hrm =
    let mzs = sufficientStatistic <$> zs
        (_,nzx,nx) = splitHarmonium hrm
        mxs = transition <$> joinAffine nx (transpose nzx) >$> mzs
        mzx = (>$<) mzs mxs
     in joinHarmonium (average mzs) mzx $ average mxs

-- Conjugation --

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedDensities0
    :: forall f z x .  ( Bilinear f z x, Map Natural f x z
                       , LegendreExponentialFamily z, LegendreExponentialFamily x )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f z x
      -> Sample z
      -> [Double]
logConjugatedDensities0 (rho0,rprms) hrm zs =
    let (nz,nzx,nx) = splitHarmonium hrm
        cnds = potential . (nx +) <$> (zs *<$< nzx)
        dts = dotMap nz $ sufficientStatistic <$> zs
        scls = [ logBaseMeasure (Proxy @ z) z - potential (nx + rprms) - rho0 | z <- zs ]
        zipper a b c = a + b + c
     in zipWith3 zipper cnds dts scls

-- | The density over the observable variables of a mixture model.
conjugatedDensity
    :: ( ConjugatedLikelihood f z x, LegendreExponentialFamily z
       , LegendreExponentialFamily x, Bilinear f z x, Map Natural f x z )
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> Double
conjugatedDensity hrm z =
    let rho0rprms = harmoniumConjugationParameters hrm
     in exp . head $ logConjugatedDensities0 rho0rprms hrm [z]

-- | The density over the observable variables of a mixture model.
conjugatedDensities
    :: ( ConjugatedLikelihood f z x, LegendreExponentialFamily z
       , LegendreExponentialFamily x, Bilinear f z x, Map Natural f x z )
    => Natural # Harmonium f z x
    -> Sample z
    -> [Double]
conjugatedDensities hrm zs =
    let rho0rprms = harmoniumConjugationParameters hrm
     in exp <$> logConjugatedDensities0 rho0rprms hrm zs

-- | The log-density over the observable variables of a mixture model.
logConjugatedDensity
    :: ( ConjugatedLikelihood f z x, LegendreExponentialFamily z
       , LegendreExponentialFamily x, Bilinear f z x, Map Natural f x z )
    => Natural # Harmonium f z x
    -> SamplePoint z
    -> [Double]
logConjugatedDensity hrm z =
    let rho0rprms = harmoniumConjugationParameters hrm
     in logConjugatedDensities0 rho0rprms hrm [z]

-- | The log-density over the observable variables of a mixture model.
logConjugatedDensities
    :: ( ConjugatedLikelihood f z x, LegendreExponentialFamily z
       , LegendreExponentialFamily x, Bilinear f z x, Map Natural f x z )
    => Natural # Harmonium f z x
    -> Sample z
    -> Double
logConjugatedDensities hrm zs =
    let rho0rprms = harmoniumConjugationParameters hrm
     in head $ logConjugatedDensities0 rho0rprms hrm zs

-- | The conjugation parameters of a conjugated `Harmonium`.
harmoniumConjugationParameters
    :: ConjugatedLikelihood f z x
    => Natural # Harmonium f z x -- ^ Categorical likelihood
    -> (Double, Natural # x) -- ^ Conjugation parameters
harmoniumConjugationParameters hrm =
    conjugationParameters . fst $ splitBottomHarmonium hrm

-- | The conjugation parameters of a conjugated `Harmonium`.
conjugatedPrior
    :: ConjugatedLikelihood f z x
    => Natural # Harmonium f z x -- ^ Categorical likelihood
    -> Natural # x -- ^ Conjugation parameters
conjugatedPrior hrm =
    let (affzx,nx) = splitBottomHarmonium hrm
        cx = snd $ conjugationParameters affzx
     in nx + cx

-- | The conjugation parameters of a conjugated `Harmonium`.
sampleConjugated
    :: ( ConjugatedLikelihood f z x, Generative Natural z
       , ExponentialFamily x, Generative Natural x )
    => Int
    -> Natural # Harmonium f z x -- ^ Categorical likelihood
    -> Random r (Sample (z,x)) -- ^ Conjugation parameters
sampleConjugated n hrm = do
    let (affzx,nx) = splitBottomHarmonium hrm
        nx' = nx + snd (conjugationParameters affzx)
    xs <- sample n nx'
    zs <- mapM samplePoint $ affzx >$>* xs
    return $ zip zs xs

-- | The conjugation parameters of a conjugated `Harmonium`.
conjugatedPotential
    :: (LegendreExponentialFamily x, ConjugatedLikelihood f z x)
    => Natural # Harmonium f z x -- ^ Categorical likelihood
    -> Double -- ^ Conjugation parameters
conjugatedPotential hrm = do
    let (affzx,nx) = splitBottomHarmonium hrm
        (rho0,rprms) = conjugationParameters affzx
     in potential (nx + rprms) + rho0


--- Internal ---


mixtureLikelihoodConjugationParameters
    :: (KnownNat k, LegendreExponentialFamily z)
    => Natural # z <* Categorical k -- ^ Categorical likelihood
    -> (Double, Natural # Categorical k) -- ^ Conjugation parameters
mixtureLikelihoodConjugationParameters aff =
    let (nz,nzx) = splitAffine aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz + Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)

harmoniumLogBaseMeasure
    :: forall f z x . (ExponentialFamily z, ExponentialFamily x)
    => Proxy (Harmonium f z x)
    -> SamplePoint (z,x)
    -> Double
harmoniumLogBaseMeasure _ (z,x) =
    logBaseMeasure (Proxy @ z) z + logBaseMeasure (Proxy @ x) x


--- Instances ---


instance (Manifold z, Manifold x, Manifold (f z x)) => Manifold (Harmonium f z x) where
      type Dimension (Harmonium f z x) = Dimension (f z x) + Dimension z + Dimension x

instance Manifold (Harmonium f z x) => Statistical (Harmonium f z x) where
    type SamplePoint (Harmonium f z x) = SamplePoint (z,x)

instance ( ExponentialFamily z, ExponentialFamily x, Bilinear f z x )
  => ExponentialFamily (Harmonium f z x) where
      sufficientStatistic (z,x) =
          let mz = sufficientStatistic z
              mx = sufficientStatistic x
           in joinHarmonium mz (mz >.< mx) mx
      averageSufficientStatistic zxs =
          let (zs,xs) = unzip zxs
              mzs = sufficientStatistic <$> zs
              mxs = sufficientStatistic <$> xs
           in joinHarmonium (average mzs) (mzs >$< mxs) (average mxs)
      logBaseMeasure = harmoniumLogBaseMeasure

instance (KnownNat n, LegendreExponentialFamily z)
         => ConjugatedLikelihood Tensor z (Categorical n) where
    conjugationParameters = conjugationParameters

instance ( KnownNat n, LegendreExponentialFamily z, Generative Natural z, Manifold (Mixture z n ) )
         => Generative Natural (Mixture z n) where
    sample = sampleConjugated

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
      potential = conjugatedPotential

instance (KnownNat n, DuallyFlatExponentialFamily z) => DuallyFlat (Mixture z n) where
    dualPotential mhrm =
        let nhrm = toNatural mhrm
         in mhrm <.> nhrm - potential nhrm

instance (KnownNat n, LegendreExponentialFamily z)
  => AbsolutelyContinuous Natural (Mixture z n) where
    densities = exponentialFamilyDensities

--instance ( LegendreExponentialFamily z, KnownNat n, SamplePoint z ~ t )
--  => LogLikelihood Natural (Mixture z n) t where
--    logLikelihood xs hrm =
--        let rh0rx = mixtureLikelihoodConjugationParameters . fst $ splitBottomHarmonium hrm
--         in average $ logConjugatedDensities rh0rx hrm xs
--    logLikelihoodDifferential zs hrm =
--        let pxs = expectationStep zs hrm
--            qxs = transition hrm
--         in pxs - qxs

instance ( ExponentialFamily z, Map Natural f x z, Bilinear f z x
         , LegendreExponentialFamily x )
  => ExpectationMaximization Natural (Harmonium f) z x where
      expectationStep = harmoniumExpectationStep

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE TypeApplications,UndecidableInstances #-}
-- | An Exponential Family 'Harmonium' is a product exponential family with a
-- particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
-- A 'Mixture' model is a special case of harmonium.
module Goal.Graphical.Generative.Harmonium where
--    (
--    -- * Harmoniums
--      Harmonium
--    -- ** Constuction
--    , splitHarmonium
--    , splitBottomHarmonium
--    , joinHarmonium
--    , joinBottomHarmonium
--    -- ** Manipulation
--    , transposeHarmonium
--    -- ** Evaluation
--    , unnormalizedHarmoniumObservableDensity
--    , unnormalizedHarmoniumObservableLogDensity
--    -- ** Sampling
--    , initialPass
--    , gibbsPass
--    -- ** Mixture Models
--    , Mixture
--    , joinNaturalMixture
--    , splitNaturalMixture
--    , joinMeanMixture
--    , splitMeanMixture
--    -- ** Conjugated Harmoniums
--    , ConjugatedLikelihood (conjugationParameters)
--    , joinConjugatedHarmonium
--    , splitConjugatedHarmonium
--    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Generative

import qualified Goal.Core.Vector.Storable as S
--import qualified Data.Vector.Storable as VS


--- Types ---


-- | A 2-layer harmonium.
newtype Harmonium f y x z w = Harmonium (Affine f y z x, w)

deriving instance (Manifold z, Manifold (f y x), Manifold w)
  => Manifold (Harmonium f y x z w)
deriving instance (Manifold z, Manifold (f y x), Manifold w)
  => Product (Harmonium f y x z w)

-- | A 'Mixture' model is simply a 'Harmonium' where the latent variable is
-- 'Categorical'.
type Mixture z k = Harmonium Tensor z (Categorical k) z (Categorical k)
type LocationMixture z y k = Harmonium Tensor y (Categorical k) z (Categorical k)


--- Classes ---


-- | The conjugation parameters of a conjugated likelihood.
class (Manifold z, Manifold x, Manifold (f z x))
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
    -> c # Harmonium f y x z w -- ^ Harmonium
joinHarmonium nz nyx = join (join nz nyx)

-- | Splits a 'Harmonium' into component parameters.
splitHarmonium
    :: (Manifold z, Manifold (f y x), Manifold w)
    => c # Harmonium f y x z w -- ^ Harmonium
    -> (c # z, c # f y x, c # w) -- ^ Biases and interaction parameters
splitHarmonium hrm =
    let (fzx,nw) = split hrm
        (nz,nyx) = split fzx
     in (nz,nyx,nw)

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
    => Natural # Harmonium f y x z w
    -> Natural # Harmonium f x y w z
transposeHarmonium hrm =
        let (nz,nyx,nw) = splitHarmonium hrm
         in joinHarmonium nw (transpose nyx) nz

-- Evaluation --

-- | The unnormalized density of a given 'Harmonium' 'Point'.
unnormalizedHarmoniumObservableLogDensity
    :: forall f y x z w
    . ( Manifold w, ExponentialFamily z, ExponentialFamily y
      , LegendreExponentialFamily w, Translation w x
      , Map Natural f x y, Bilinear f y x, SamplePoint y ~ SamplePoint z )
    => Natural # Harmonium f y x z w
    -> SamplePoint z
    -> Double
unnormalizedHarmoniumObservableLogDensity hrm z =
    let (pstr,nz) = split $ transposeHarmonium hrm
        mz = sufficientStatistic z
     in nz <.> mz + potential (pstr >.>* z) + logBaseMeasure (Proxy @ z) z

-- | Computes the joint expectations of a harmonium based on a sample from the
-- observable layer.
harmoniumExpectationStep
    :: ( ExponentialFamily z, Map Natural f x y, Manifold w
       , Translation z y, Translation w x
       , Bilinear f y x, LegendreExponentialFamily w )
    => Sample z -- ^ Model Samples
    -> Natural # Harmonium f y x z w -- ^ Harmonium
    -> Mean # Harmonium f y x z w -- ^ Harmonium expected sufficient statistics
harmoniumExpectationStep zs hrm =
    let mzs = sufficientStatistic <$> zs
        mys = anchor <$> mzs
        pstr = fst . split $ transposeHarmonium hrm
        mws = transition <$> pstr >$> mys
        mxs = anchor <$> mws
        myx = (>$<) mys mxs
     in joinHarmonium (average mzs) myx $ average mws

---- Sampling --

initialPass
    :: ( ExponentialFamily z, Map Natural f x y, Manifold w
       , Translation z y, Translation w x, SamplePoint z ~ SamplePoint y
       , Generative Natural w, ExponentialFamily y
       , Bilinear f y x, LegendreExponentialFamily w )
    => Natural # Harmonium f y x z w -- ^ Harmonium
    -> Sample z -- ^ Model Samples
    -> Random s (Sample (z, w))
initialPass hrm zs = do
    let pstr = fst . split $ transposeHarmonium hrm
    xs <- mapM samplePoint $ pstr >$>* zs
    return $ zip zs xs

gibbsPass
    :: ( ExponentialFamily z, Map Natural f x y, Manifold w
       , Translation z y, Translation w x, SamplePoint z ~ SamplePoint y
       , Generative Natural w, ExponentialFamily y, SamplePoint x ~ SamplePoint w
       , Map Natural f y x, ExponentialFamily x, Generative Natural z
       , Bilinear f y x, LegendreExponentialFamily w )
    => Natural # Harmonium f y x z w -- ^ Harmonium
    -> Sample (z, w)
    -> Random s (Sample (z, w))
gibbsPass hrm zws = do
    let ws = snd <$> zws
        pstr = fst . split $ transposeHarmonium hrm
        lkl = fst $ split hrm
    zs' <- mapM samplePoint $ lkl >$>* ws
    ws' <- mapM samplePoint $ pstr >$>* zs'
    return $ zip zs' ws'

-- Conjugation --

--- | Computes the negative log-likelihood of a sample point of a conjugated harmonium.
logConjugatedDensities0
    :: forall f x y z w
    . ( Bilinear f z x, Map Natural f x z, SamplePoint z ~ SamplePoint y
      , LegendreExponentialFamily z, LegendreExponentialFamily x )
      => (Double, Natural # x) -- ^ Conjugation Parameters
      -> Natural # Harmonium f y x z w
      -> Sample z
      -> [Double]
logConjugatedDensities0 (rho0,rprms) hrm zs =
    let (pstr,nz) = split $ transposeHarmonium hrm
        cnds = potential <$> pstr >$>* zs
        nx = snd $ split hrm
        dts = dotMap nz $ sufficientStatistic <$> zs
        scls = [ logBaseMeasure (Proxy @ z) z - potential (nx + rprms) - rho0 | z <- zs ]
        zipper a b c = a + b + c
     in zipWith3 zipper cnds dts scls

---- | The conjugation parameters of a conjugated `Harmonium`.
--harmoniumConjugationParameters
--    :: ConjugatedLikelihood f z x
--    => Natural # Harmonium f z x -- ^ Categorical likelihood
--    -> (Double, Natural # x) -- ^ Conjugation parameters
--harmoniumConjugationParameters hrm =
--    conjugationParameters . fst $ splitBottomHarmonium hrm
--
---- | The conjugation parameters of a conjugated `Harmonium`.
--splitConjugatedHarmonium
--    :: ConjugatedLikelihood f z x
--    => Natural # Harmonium f z x -- ^ Categorical likelihood
--    -> (Natural # Affine f z x,Natural # x) -- ^ Conjugation parameters
--splitConjugatedHarmonium hrm =
--    let (affzx,nx) = splitBottomHarmonium hrm
--        cx = snd $ conjugationParameters affzx
--     in (affzx,nx + cx)
--
---- | The conjugation parameters of a conjugated `Harmonium`.
--joinConjugatedHarmonium
--    :: ConjugatedLikelihood f z x
--    => Natural # Affine f z x -- ^ Conjugation parameters
--    -> Natural # x
--    -> Natural # Harmonium f z x -- ^ Categorical likelihood
--joinConjugatedHarmonium affzx nx =
--    let cx = snd $ conjugationParameters affzx
--     in joinBottomHarmonium affzx $ nx - cx
--
---- | The conjugation parameters of a conjugated `Harmonium`.
--sampleConjugated
--    :: ( ConjugatedLikelihood f z x, Generative Natural z
--       , ExponentialFamily x, Generative Natural x )
--    => Int
--    -> Natural # Harmonium f z x -- ^ Categorical likelihood
--    -> Random r (Sample (z,x)) -- ^ Conjugation parameters
--sampleConjugated n hrm = do
--    let (affzx,nx) = splitBottomHarmonium hrm
--        nx' = nx + snd (conjugationParameters affzx)
--    xs <- sample n nx'
--    zs <- mapM samplePoint $ affzx >$>* xs
--    return $ zip zs xs
--
---- | The conjugation parameters of a conjugated `Harmonium`.
--conjugatedPotential
--    :: (LegendreExponentialFamily x, ConjugatedLikelihood f z x)
--    => Natural # Harmonium f z x -- ^ Categorical likelihood
--    -> Double -- ^ Conjugation parameters
--conjugatedPotential hrm = do
--    let (affzx,nx) = splitBottomHarmonium hrm
--        (rho0,rprms) = conjugationParameters affzx
--     in potential (nx + rprms) + rho0


--- Internal ---


mixtureLikelihoodConjugationParameters
    :: (KnownNat k, LegendreExponentialFamily z)
    => Natural # z <* Categorical k -- ^ Categorical likelihood
    -> (Double, Natural # Categorical k) -- ^ Conjugation parameters
mixtureLikelihoodConjugationParameters aff =
    let (nz,nzx) = split aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz + nzxi) $ toColumns nzx
     in (rho0, Point rprms)

harmoniumLogBaseMeasure
    :: forall f y x z w . (ExponentialFamily z, ExponentialFamily w)
    => Proxy (Harmonium f y x z w)
    -> SamplePoint (z,w)
    -> Double
harmoniumLogBaseMeasure _ (z,w) =
    logBaseMeasure (Proxy @ z) z + logBaseMeasure (Proxy @ w) w


--- Instances ---


instance Manifold (Harmonium f y x z w) => Statistical (Harmonium f y x z w) where
    type SamplePoint (Harmonium f y x z w) = SamplePoint (z,w)

instance ( ExponentialFamily z, ExponentialFamily x, Translation z y
         , Translation w x
         , ExponentialFamily w, ExponentialFamily y, Bilinear f y x )
  => ExponentialFamily (Harmonium f y x z w) where
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

instance (KnownNat k, LegendreExponentialFamily z)
         => ConjugatedLikelihood Tensor z (Categorical k) z (Categorical k) where
    conjugationParameters = mixtureLikelihoodConjugationParameters

--instance ( KnownNat k, LegendreExponentialFamily z, Generative Natural z, Manifold (Mixture z n ) )
--         => Generative Natural (Mixture z n) where
--    sample = sampleConjugated

--instance (KnownNat k, LegendreExponentialFamily z) => Transition Natural Mean (Mixture z k) where
--    transition nhrm =
--        let (nzs,nx) = splitNaturalMixture nhrm
--            mx = toMean nx
--            mzs = S.map transition nzs
--         in joinMeanMixture mzs mx
--
--instance (KnownNat k, DuallyFlatExponentialFamily z) => Transition Mean Natural (Mixture z k) where
--    transition mhrm =
--        let (mzs,mx) = splitMeanMixture mhrm
--            nx = transition mx
--            nzs = S.map transition mzs
--         in joinNaturalMixture nzs nx
--
--instance (KnownNat k, LegendreExponentialFamily z) => Legendre (Mixture z n) where
--      type PotentialCoordinates (Mixture z n) = Natural
--      potential = conjugatedPotential
--
--instance (KnownNat k, DuallyFlatExponentialFamily z) => DuallyFlat (Mixture z n) where
--    dualPotential mhrm =
--        let nhrm = toNatural mhrm
--         in mhrm <.> nhrm - potential nhrm
--
--instance (KnownNat k, LegendreExponentialFamily z)
--  => AbsolutelyContinuous Natural (Mixture z n) where
--    logDensities = exponentialFamilyLogDensities
--
--instance (KnownNat k, LegendreExponentialFamily z)
--  => ObservablyContinuous Natural (Harmonium Tensor) z (Categorical k) where
--    logObservableDensities hrm zs =
--        let rho0rprms = harmoniumConjugationParameters hrm
--         in logConjugatedDensities0 rho0rprms hrm zs
--
--instance ( LegendreExponentialFamily z, KnownNat k, SamplePoint z ~ t )
--  => LogLikelihood Natural (Mixture z n) t where
--    logLikelihood xs hrm =
--         average $ logObservableDensities hrm xs
--    logLikelihoodDifferential zs hrm =
--        let pxs = expectationStep zs hrm
--            qxs = transition hrm
--         in pxs - qxs
--
--instance ( ExponentialFamily z, Map Natural f z x )
--  => Map Natural (Harmonium f) z x where
--      (>.>) hrm = (>.>) (fst $ splitBottomHarmonium hrm)
--      (>$>) hrm = (>$>) (fst $ splitBottomHarmonium hrm)
--
--instance ( ExponentialFamily z, Map Natural f x z, Bilinear f z x
--         , LegendreExponentialFamily x )
--  => ExpectationMaximization (Harmonium f) z x where
--      expectationStep = harmoniumExpectationStep

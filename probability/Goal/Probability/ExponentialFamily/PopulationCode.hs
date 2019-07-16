{-# LANGUAGE TypeApplications #-}

-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.PopulationCode
    ( -- * Synonyms
      Neurons
    , Response
    -- * Population Encoders
    , joinNormalPopulationEncoder
    , joinVonMisesPopulationEncoder
    , splitVonMisesPopulationEncoder
    , joinVonMisesMixturePopulationEncoder
    , splitVonMisesMixturePopulationEncoder
    , sortVonMisesMixturePopulationEncoder
    , sortedVonMisesMixturePopulationIndices
    , sortVonMisesMixturePopulationOnIndices
    -- * Learning
    , mixturePopulationPartialExpectationMaximization
    -- * Utility
    , tuningCurves
    , mixturePopulationCovariance
    , mixturePopulationNoiseCorrelations
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conditional

import qualified Data.List as L


--- Population Codes ---


type Neurons k = Replicated k Poisson
type Response k = SamplePoint (Neurons k)


--- Population Encoders ---


-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
joinNormalPopulationEncoder
    :: KnownNat k
    => Either Double (Natural # Neurons k) -- ^ Global Gain or Gains
    -> S.Vector k (Natural # Normal) -- ^ Tuning Curves
    -> Natural #> Replicated k Poisson <* Normal -- ^ Population Encoder
joinNormalPopulationEncoder engns nps =
    let mtx = fromRows nps
        sz0 = case engns of
                (Left ngn) -> Point (S.replicate ngn)
                (Right ngns) -> ngns
        nz1 = sz0 <+> Point (S.map normalBias nps)
        nz = nz1 <-> Point (S.map potential nps)
     in joinAffine nz mtx

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
joinVonMisesPopulationEncoder
    :: KnownNat k
    => Either Double (Natural # Neurons k) -- ^ Global Gain or Gains
    -> S.Vector k (Natural # VonMises) -- ^ Von Mises Curves
    -> Natural #> Replicated k Poisson <* VonMises -- ^ Population Encoder
joinVonMisesPopulationEncoder engns nps =
    let mtx = fromRows nps
        nz0 = case engns of
                (Left ngn) -> Point (S.replicate ngn)
                (Right ngns) -> ngns
        nz = nz0 <-> Point (S.map potential nps)
     in joinAffine nz mtx

-- | Splits a von mises population code.
splitVonMisesPopulationEncoder
    :: KnownNat k
    => Natural #> Replicated k Poisson <* VonMises -- ^ Population Encoder
    -> (Natural # Neurons k, S.Vector k (Natural # VonMises)) -- ^ Gains and Von Mises Curves
splitVonMisesPopulationEncoder lkl =
    let (nz0,nzx) = splitAffine lkl
        nxs = toRows nzx
        nz = nz0 <+> Point (S.map potential nxs)
     in (nz,nxs)

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
joinVonMisesMixturePopulationEncoder
    :: (KnownNat k, KnownNat n)
    => Natural # Categorical n -- ^ Weights
    -> S.Vector (n+1) (Natural # Neurons k) -- ^ Gain components
    -> S.Vector k (Natural # VonMises) -- ^ Von Mises Curves
    -> Natural #> ConditionalMixture (Neurons k) n VonMises -- ^ Mixture Encoder
joinVonMisesMixturePopulationEncoder nk ngnss nps =
    let nzx = fromRows nps
        nzk = S.map (<-> Point (S.map potential nps)) ngnss
     in joinConditionalDeepHarmonium (joinMixture nzk nk) nzx

splitVonMisesMixturePopulationEncoder
    :: (KnownNat k, KnownNat n)
    => Natural #> ConditionalMixture (Neurons k) n VonMises
    -> ( Natural # Categorical n
       , S.Vector (n+1) (Natural # Neurons k)
       , S.Vector k (Natural # VonMises) )
splitVonMisesMixturePopulationEncoder mlkl =
    let (mxmdl,nzx) = splitConditionalDeepHarmonium mlkl
        (nzk,nk) = splitMixture mxmdl
        nps = toRows nzx
        ngnss = S.map (<+> Point (S.map potential nps)) nzk
     in (nk,ngnss,nps)

sortVonMisesMixturePopulationEncoder
    :: forall k n . (KnownNat k, KnownNat n)
    => Natural #> ConditionalMixture (Neurons k) n VonMises
    -> Natural #> ConditionalMixture (Neurons k) n VonMises
sortVonMisesMixturePopulationEncoder mlkl =
    let (wghts,gnss,tcs) = splitVonMisesMixturePopulationEncoder mlkl
        mus = head . listCoordinates . toSource <$> S.toList tcs
        idxs :: S.Vector k Int
        idxs = fromJust . S.fromList . map fst . L.sortOn snd $ zip [0..] mus
        tcs' = S.backpermute tcs idxs
        gnss' = S.map (Point . flip S.backpermute idxs . coordinates) gnss
     in joinVonMisesMixturePopulationEncoder wghts gnss' tcs'

sortedVonMisesMixturePopulationIndices
    :: forall k n . (KnownNat k, KnownNat n)
    => Natural #> ConditionalMixture (Neurons k) n VonMises
    -> S.Vector k Int
sortedVonMisesMixturePopulationIndices mlkl =
    let (_,_,tcs) = splitVonMisesMixturePopulationEncoder mlkl
        mus = head . listCoordinates . toSource <$> S.toList tcs
     in fromJust . S.fromList . map fst . L.sortOn snd $ zip [0..] mus

sortVonMisesMixturePopulationOnIndices
    :: forall k n . (KnownNat k, KnownNat n)
    => S.Vector k Int
    -> Natural #> ConditionalMixture (Neurons k) n VonMises
    -> Natural #> ConditionalMixture (Neurons k) n VonMises
sortVonMisesMixturePopulationOnIndices idxs mlkl =
    let (wghts,gnss,tcs) = splitVonMisesMixturePopulationEncoder mlkl
        tcs' = S.backpermute tcs idxs
        gnss' = S.map (Point . flip S.backpermute idxs . coordinates) gnss
     in joinVonMisesMixturePopulationEncoder wghts gnss' tcs'

-- | Stimulus Dependent Noise Correlations, ordered by preferred stimulus.
mixturePopulationCovariance
    :: forall k n . ( KnownNat k, KnownNat n )
    => Natural # Mixture (Neurons k) n -- ^ Mixture Encoder
    -> Source # MultivariateNormal k -- ^ Mean Parameter Correlations
{-# INLINE mixturePopulationCovariance #-}
mixturePopulationCovariance mxmdl =
    let (ngnss, nwghts) = splitMixture mxmdl
        wghts0 = coordinates $ toSource nwghts
        wghts = 1 - S.sum wghts0 : S.toList wghts0
        gnss = toMean <$> S.toList ngnss
        mgns = weightedAveragePoint $ zip wghts gnss
        mgns2 :: Mean #> Tensor (Neurons k) (Neurons k)
        mgns2 = mgns >.< mgns
        mmgns = weightedAveragePoint $ zip wghts [ gns >.< gns | gns <- gnss ]
        cvgns = mmgns <-> mgns2
        cvnrns = cvgns <+> (fromMatrix . S.diagonalMatrix $ coordinates mgns)
     in joinMultivariateNormal (coordinates mgns) (toMatrix cvnrns)

-- | Stimulus Dependent Noise Correlations, ordered by preferred stimulus.
mixturePopulationNoiseCorrelations
    :: forall k n . ( KnownNat k, KnownNat n )
    => Natural # Mixture (Neurons k) n -- ^ Mixture Encoder
    -> S.Matrix k k Double -- ^ Mean Parameter Correlations
{-# INLINE mixturePopulationNoiseCorrelations #-}
mixturePopulationNoiseCorrelations =
    multivariateNormalCorrelations . mixturePopulationCovariance


--- Expectation Maximization


-- | EM implementation for mixture models/categorical harmoniums.
mixturePopulationPartialExpectationMaximization
    :: ( KnownNat n, KnownNat k, ExponentialFamily x )
    => Sample (Neurons k, x) -- ^ Observations
    -> Natural #> ConditionalMixture (Neurons k) n x -- ^ Current Mixture GLM
    -> Natural #> ConditionalMixture (Neurons k) n x -- ^ Updated Mixture GLM
{-# INLINE mixturePopulationPartialExpectationMaximization #-}
mixturePopulationPartialExpectationMaximization zxs mlkl =
    let (hrm,tcs) = splitConditionalDeepHarmonium mlkl
        wghts = snd $ splitMixture hrm
        gnss = S.map (Point @ Mean . S.map (max 1e-20) . coordinates) $ mixturePopulationPartialMStep zxs tcs hrm
        hrm' = joinMixture (S.map transition gnss) wghts
     in joinConditionalDeepHarmonium hrm' tcs

-- | M-step implementation for deep mixture models/categorical harmoniums. Note
-- that for the sake of type signatures, this acts on transposed harmoniums
-- (i.e. the categorical variables are at the bottom of the hierarchy).
mixturePopulationPartialMStep
    :: ( KnownNat n, KnownNat k, ExponentialFamily x )
    => Sample (Neurons k, x) -- ^ Observations
    -> Natural #> Tensor (Neurons k) x
    -> Natural # Mixture (Neurons k) n
    -> S.Vector (n+1) (Mean # Neurons k)
{-# INLINE mixturePopulationPartialMStep #-}
mixturePopulationPartialMStep zxs tcs hrm =
    let (zs,xs) = unzip zxs
        aff = fst . splitBottomHarmonium $ transposeHarmonium hrm
        szs = sufficientStatistic <$> zs
        wghtss = toMean <$> aff >$> szs
        mzs = toMean <$> tcs >$>* xs
        (cmpnts0,nrms) = foldr folder (S.replicate zero, S.replicate zero) $ zip3 wghtss szs mzs
     in S.zipWith divider cmpnts0 nrms
    where folder (wghts,sz,mz) (cmpnts,nrms) =
              let ws = categoricalWeights wghts
                  cmpnts' = S.map (.> sz) ws
                  nrms' = S.map (.> mz) ws
               in (S.zipWith (<+>) cmpnts cmpnts', S.zipWith (<+>) nrms' nrms)
          divider (Point p) (Point q) = Point $ S.zipWith (/) p q

---- | E-step implementation for deep mixture models/categorical harmoniums. Note
---- that for the sake of type signatures, this acts on transposed harmoniums
---- (i.e. the categorical variables are at the bottom of the hierarchy).
--mixturePopulationPartialExpectationMaximization'
--    :: forall n k x . ( KnownNat n, KnownNat k, ExponentialFamily x )
--    => Sample (Neurons k, x) -- ^ Observations
--    -> Natural #> ConditionalMixture (Neurons k) n x
--    -> Natural #> ConditionalMixture (Neurons k) n x
--{-# INLINE mixturePopulationPartialExpectationMaximization' #-}
--mixturePopulationPartialExpectationMaximization zxs mlkl =
--    let (zs,xs) = unzip zxs
--        (hrm,aff) = splitConditionalDeepHarmonium mlkl
--        avgz :: Mean # Neurons k
--        avgz = sufficientStatisticT zs
--        hrmzc = transposeHarmonium hrm
--        affzc = fst $ splitBottomHarmonium hrmzc
--        mwghtss = toMean <$> affzc >$>* zs
--        fss = aff >$>* xs
--        (nz,nzc) = splitAffine affzc
--        mz = toMean nz
--        mzcs = S.map toMean $ toRows nzc
--        mz' =
--    where

-- Population Code Conjugation


---- | Given a set of conjugation parameters and a population code, modulates
---- the gains of the population code to best satisfy the resulting conjugation
---- equation. Note that this uses LLS, and can hang if the calculation would
---- produce negative gains.
--conjugatePopulationEncoder
--    :: (KnownNat k, ExponentialFamily m)
--    => Double -- ^ Conjugation shift
--    -> Natural # m -- ^ Conjugation parameters
--    -> Sample m -- ^ Sample points
--    -> Natural #> Neurons k <* m -- ^ Given PPC
--    -> Natural #> Neurons k <* m -- ^ Conjugated PPC
--{-# INLINE conjugatePopulationEncoder #-}
--conjugatePopulationEncoder rho0 rprms mus lkl =
--    let dpnds = conjugationCurve rho0 rprms mus
--        indpnds = independentVariables1 lkl mus
--        gns = Point . S.map log $ S.linearLeastSquares indpnds dpnds
--        (gns0,tcs) = splitAffine lkl
--     in joinAffine (gns0 <+> gns) tcs
--
---- | A gradient for conjugateing gains which won't allow them to be negative.
--populationEncoderConjugationDifferential
--    :: (KnownNat k, ExponentialFamily x)
--    => Double -- ^ Conjugation shift
--    -> Natural # x -- ^ Conjugation parameters
--    -> Sample x -- ^ Sample points
--    -> Natural #> Tensor (Neurons k) x -- ^ linear part of ppc
--    -> Natural # Neurons k -- ^ Gains
--    -> Mean # Neurons k -- ^ Conjugated PPC
--{-# INLINE populationEncoderConjugationDifferential #-}
--populationEncoderConjugationDifferential rho0 rprms xsmps tns ngns =
--    let lkl = joinAffine ngns tns
--        rcts = conjugationCurve rho0 rprms xsmps
--        fss = transition <$> lkl >$>* xsmps
--     in averagePoint $ do
--         (rct,fs) <- zip rcts fss
--         let sms = S.sum $ coordinates fs
--             dff = sms - rct
--         return $ dff .> fs

-- | Returns the tuning curves of a population code over a set of sample points.
-- This is often useful for plotting purposes.
tuningCurves
    :: (ExponentialFamily m, KnownNat k)
    => Sample m -- Sample points
    -> Natural #> Neurons k <* m -- ^ PPC
    -> [[(SamplePoint m, Double)]] -- ^ Vector of tuning curves
tuningCurves xsmps lkl =
    let tcs = L.transpose $ listCoordinates . toMean <$> (lkl >$>* xsmps)
     in zip xsmps <$> tcs

--- Internal ---

--independentVariables1
--    :: (KnownNat k, ExponentialFamily m)
--    => Natural #> Neurons k <* m
--    -> Sample m
--    -> [S.Vector k Double]
--independentVariables1 lkl mus =
--    coordinates . toMean <$> lkl >$>* mus

normalBias :: Natural # Normal -> Double
normalBias sp =
    let [mu,vr] = listCoordinates $ toSource sp
     in - mu^(2 :: Int)/(2*vr)



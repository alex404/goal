{-# LANGUAGE
    TypeOperators,
    DataKinds
#-}

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
    -- * Conjugation
    , conjugatePopulationEncoder
    , populationEncoderConjugationDifferential
    -- * Utility
    , tuningCurves
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
import Goal.Probability.ExponentialFamily.Harmonium.Conjugation
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
    -> Mean #> Natural # Replicated k Poisson <* Normal -- ^ Population Encoder
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
    -> Mean #> Natural # Replicated k Poisson <* VonMises -- ^ Population Encoder
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
    => Mean #> Natural # Replicated k Poisson <* VonMises -- ^ Population Encoder
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
    => Natural # Categorical Int n -- ^ Weights
    -> S.Vector (n+1) (Natural # Neurons k) -- ^ Gain components
    -> S.Vector k (Natural # VonMises) -- ^ Von Mises Curves
    -> Mean #> Natural # MixtureGLM (Neurons k) n VonMises -- ^ Mixture Encoder
joinVonMisesMixturePopulationEncoder nk ngnss nps =
    let nzx = fromRows nps
        nzk = S.map (<-> Point (S.map potential nps)) ngnss
     in joinBottomSubLinear (joinMixtureModel nzk nk) nzx

splitVonMisesMixturePopulationEncoder
    :: (KnownNat k, KnownNat n)
    => Mean #> Natural # MixtureGLM (Neurons k) n VonMises
    -> ( Natural # Categorical Int n
       , S.Vector (n+1) (Natural # Neurons k)
       , S.Vector k (Natural # VonMises) )
splitVonMisesMixturePopulationEncoder mlkl =
    let (mxmdl,nzx) = splitBottomSubLinear mlkl
        (nzk,nk) = splitMixtureModel mxmdl
        nps = toRows nzx
        ngnss = S.map (<+> Point (S.map potential nps)) nzk
     in (nk,ngnss,nps)


-- Population Code Conjugation


-- | Given a set of conjugation parameters and a population code, modulates
-- the gains of the population code to best satisfy the resulting conjugation
-- equation. Note that this uses LLS, and can hang if the calculation would
-- produce negative gains.
conjugatePopulationEncoder
    :: (KnownNat k, ExponentialFamily m)
    => Double -- ^ Conjugation shift
    -> Natural # m -- ^ Conjugation parameters
    -> Sample m -- ^ Sample points
    -> Mean #> Natural # Neurons k <* m -- ^ Given PPC
    -> Mean #> Natural # Neurons k <* m -- ^ Conjugated PPC
{-# INLINE conjugatePopulationEncoder #-}
conjugatePopulationEncoder rho0 rprms mus lkl =
    let dpnds = conjugationCurve rho0 rprms mus
        indpnds = independentVariables1 lkl mus
        gns = Point . S.map log $ S.linearLeastSquares indpnds dpnds
        (gns0,tcs) = splitAffine lkl
     in joinAffine (gns0 <+> gns) tcs

-- | A gradient for conjugateing gains which won't allow them to be negative.
populationEncoderConjugationDifferential
    :: (KnownNat k, ExponentialFamily x)
    => Double -- ^ Conjugation shift
    -> Natural # x -- ^ Conjugation parameters
    -> Sample x -- ^ Sample points
    -> Mean #> Natural # Tensor (Neurons k) x -- ^ linear part of ppc
    -> Natural # Neurons k -- ^ Gains
    -> CotangentPair Natural (Neurons k) -- ^ Conjugated PPC
{-# INLINE populationEncoderConjugationDifferential #-}
populationEncoderConjugationDifferential rho0 rprms xsmps tns ngns =
    let lkl = joinAffine ngns tns
        rcts = conjugationCurve rho0 rprms xsmps
        fss = dualTransition <$> lkl >$>* xsmps
     in joinTangentPair ngns . averagePoint $ do
         (rct,fs) <- zip rcts fss
         let sms = S.sum $ coordinates fs
             dff = sms - rct
         return . primalIsomorphism $ dff .> fs

-- | Returns the tuning curves of a population code over a set of sample points.
-- This is often useful for plotting purposes.
tuningCurves
    :: (ExponentialFamily m, KnownNat k)
    => Sample m -- Sample points
    -> Mean #> Natural # Neurons k <* m -- ^ PPC
    -> [[(SamplePoint m, Double)]] -- ^ Vector of tuning curves
tuningCurves xsmps lkl =
    let tcs = L.transpose $ listCoordinates . dualTransition <$> (lkl >$>* xsmps)
     in zip xsmps <$> tcs

--- Internal ---

independentVariables1
    :: (KnownNat k, ExponentialFamily m)
    => Mean #> Natural # Neurons k <* m
    -> Sample m
    -> [S.Vector k Double]
independentVariables1 lkl mus =
    coordinates . dualTransition <$> lkl >$>* mus

normalBias :: Natural # Normal -> Double
normalBias sp =
    let [mu,vr] = listCoordinates $ toSource sp
     in - mu^(2 :: Int)/(2*vr)



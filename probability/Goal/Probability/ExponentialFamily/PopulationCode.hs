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
    , normalPopulationEncoder
    , vonMisesPopulationEncoder
    , splitVonMisesPopulationEncoder
    , vonMisesMixturePopulationEncoder
    -- * Rectification
    , rectifyPopulationCode
    , populationCodeRectificationDifferential
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
import Goal.Probability.ExponentialFamily.Rectification
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
normalPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize tuning curves
    -> Either Double (Natural # Neurons k) -- ^ Global Gain or Gains
    -> S.Vector k (Natural # Normal) -- ^ Tuning Curves
    -> Mean #> Natural # Replicated k Poisson <* Normal -- ^ Population Encoder
normalPopulationEncoder nrmb engns nps =
    let mtx = fromRows nps
        sz0 = case engns of
                (Left ngn) -> Point (S.replicate ngn)
                (Right ngns) -> ngns
        nz1 = sz0 <+> Point (S.map normalBias nps)
        nz = if nrmb then  nz1 <-> Point (S.map potential nps) else nz1
     in joinAffine nz mtx

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize tuning curves
    -> Either Double (Natural # Neurons k) -- ^ Global Gain or Gains
    -> S.Vector k (Natural # VonMises) -- ^ Von Mises Curves
    -> Mean #> Natural # Replicated k Poisson <* VonMises -- ^ Population Encoder
vonMisesPopulationEncoder nrmb engns nps =
    let mtx = fromRows nps
        nz0 = case engns of
                (Left ngn) -> Point (S.replicate ngn)
                (Right ngns) -> ngns
        nz = if nrmb then nz0 <-> Point (S.map potential nps) else nz0
     in joinAffine nz mtx

-- | Splits a von mises population code.
splitVonMisesPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize gains
    -> Mean #> Natural # Replicated k Poisson <* VonMises -- ^ Population Encoder
    -> (Natural # Neurons k, S.Vector k (Natural # VonMises)) -- ^ Gains and Von Mises Curves
splitVonMisesPopulationEncoder nrmb lkl =
    let (nz0,nzx) = splitAffine lkl
        nxs = toRows nzx
        nz = if nrmb then nz0 <+> Point (S.map potential nxs) else nz0
     in (nz,nxs)

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesMixturePopulationEncoder
    :: (KnownNat k, KnownNat n)
    => Bool -- ^ Normalize tuning curves?
    -> Natural # Categorical Int n -- ^ Weights
    -> S.Vector (n+1) (Natural # Neurons k) -- ^ Gain components
    -> S.Vector k (Natural # VonMises) -- ^ Von Mises Curves
    -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises -- ^ Mixture Encoder
vonMisesMixturePopulationEncoder nrmb nk ngnss nps =
    let nzx = fromRows nps
        nzk = if nrmb then S.map (<-> Point (S.map potential nps)) ngnss else ngnss
     in joinBottomSubLinear (buildMixtureModel nzk nk) nzx


-- Population Code Rectification


-- | Given a set of rectification parameters and a population code, modulates
-- the gains of the population code to best satisfy the resulting rectification
-- equation. Note that this uses LLS, and can hang if the calculation would
-- produce negative gains.
rectifyPopulationCode
    :: (KnownNat k, ExponentialFamily m)
    => Double -- ^ Rectification shift
    -> Natural # m -- ^ Rectification parameters
    -> Sample m -- ^ Sample points
    -> Mean #> Natural # Neurons k <* m -- ^ Given PPC
    -> Mean #> Natural # Neurons k <* m -- ^ Rectified PPC
{-# INLINE rectifyPopulationCode #-}
rectifyPopulationCode rho0 rprms mus lkl =
    let dpnds = rectificationCurve rho0 rprms mus
        indpnds = independentVariables1 lkl mus
        gns = Point . S.map log $ S.linearLeastSquares indpnds dpnds
        (gns0,tcs) = splitAffine lkl
     in joinAffine (gns0 <+> gns) tcs

-- | A gradient for rectifying gains which won't allow them to be negative.
populationCodeRectificationDifferential
    :: (KnownNat k, ExponentialFamily x)
    => Double -- ^ Rectification shift
    -> Natural # x -- ^ Rectification parameters
    -> Sample x -- ^ Sample points
    -> Mean #> Natural # Tensor (Neurons k) x -- ^ linear part of ppc
    -> Natural # Neurons k -- ^ Gains
    -> CotangentPair Natural (Neurons k) -- ^ Rectified PPC
{-# INLINE populationCodeRectificationDifferential #-}
populationCodeRectificationDifferential rho0 rprms xsmps tns ngns =
    let lkl = joinAffine ngns tns
        rcts = rectificationCurve rho0 rprms xsmps
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



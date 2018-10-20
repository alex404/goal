-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.PopulationCode
    ( -- * Population Encoders
      normalPopulationEncoder
    , vonMisesPopulationEncoder
    -- * Rectification
    , populationCodeRectificationParameters
    , rectifyPopulationCode
    , rectificationCurve
    , populationCodeRectificationDifferential
    -- * Utility
    , tuningCurves
    , sumOfTuningCurves
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Data.List as L

--- Population Encoders ---

-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
normalPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize tuning curves
    -> Either Double (S.Vector k Double) -- ^ Global Gain or Gains
    -> S.Vector k (Point Source Normal) -- ^ Tuning Curves
    -> Function Mean Natural # Replicated k Poisson <* Normal -- ^ Population Encoder
normalPopulationEncoder nrmb egns sps =
    let nps = S.map toNatural sps
        mtx = S.concat $ S.map coordinates nps
        ob0 = case egns of
                (Left gn) -> S.map ((+ log gn) . normalBias) sps
                (Right gns) -> S.zipWith (+) (log gns) $ S.map normalBias sps
        ob = if nrmb then S.zipWith (-) ob0 $ S.map potential nps else ob0
     in Point $ ob S.++ mtx

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize tuning curves
    -> Either Double (S.Vector k Double) -- ^ Global Gain Gains
    -> S.Vector k (Point Source VonMises) -- ^ Von Mises Curves
    -> Function Mean Natural # Replicated k Poisson <* VonMises -- ^ Population Encoder
vonMisesPopulationEncoder nrmb egns sps =
    let ob0 = case egns of
                (Left gn) -> S.replicate $ log gn
                (Right gns) -> log gns
        nps = S.map toNatural sps
        mtx = S.concat $ S.map coordinates nps
        ob = if nrmb then S.zipWith (-) ob0 $ S.map potential nps else ob0
     in Point $ ob S.++ mtx


-- (De)Convolution of population codes

-- | Computes the rectification curve given a set of rectification parameters,
-- at the given set of points.
rectificationCurve
    :: ExponentialFamily m
    => Double -- ^ Rectification shift
    -> Natural # m -- ^ Rectification parameters
    -> Sample m -- ^ Samples points
    -> [Double] -- ^ Rectification curve at sample points
{-# INLINE rectificationCurve #-}
rectificationCurve rho0 rprms mus = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- | Given a set of rectification parameters and a population code, modulates
-- the gains of the population code to best satisfy the resulting rectification
-- equation. Note that this uses LLS, and can hang if the calculation would
-- produce negative gains.
rectifyPopulationCode
    :: (KnownNat k, ExponentialFamily m)
    => Double -- ^ Rectification shift
    -> Natural # m -- ^ Rectification parameters
    -> Sample m -- ^ Sample points
    -> Mean #> Natural # R k Poisson <* m -- ^ Given PPC
    -> Mean #> Natural # R k Poisson <* m -- ^ Rectified PPC
{-# INLINE rectifyPopulationCode #-}
rectifyPopulationCode rho0 rprms mus lkl =
    let dpnds = rectificationCurve rho0 rprms mus
        indpnds = independentVariables1 lkl mus
        gns = Point . S.map log $ S.linearLeastSquares indpnds dpnds
        (gns0,tcs) = splitAffine lkl
     in joinAffine (gns0 <+> gns) tcs

-- | A gradient for rectifying gains which won't allow them to be negative.
populationCodeRectificationDifferential
    :: (KnownNat k, ExponentialFamily m)
    => Double -- ^ Rectification shift
    -> Natural # m -- ^ Rectification parameters
    -> Sample m -- ^ Sample points
    -> Mean #> Natural # Tensor (R k Poisson)  m -- ^ linear part of ppc
    -> Natural # R k Poisson -- ^ Given PPC
    -> CotangentPair Natural (R k Poisson) -- ^ Rectified PPC
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

-- Linear Least Squares

-- | Returns the rectification parameters which best satisfy the rectification
-- equation for the given population code.
populationCodeRectificationParameters
    :: (KnownNat k, ExponentialFamily m)
    => Mean #> Natural # R k Poisson <* m -- ^ PPC
    -> Sample m -- ^ Sample points
    -> (Double, Natural # m) -- ^ Approximate rectification parameters
{-# INLINE populationCodeRectificationParameters #-}
populationCodeRectificationParameters lkl mus =
    let dpnds = sumOfTuningCurves lkl mus
        indpnds = independentVariables0 lkl mus
        (rho0,rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

-- | The sum of the tuning curves of a population over a sample.
sumOfTuningCurves
    :: (KnownNat k, ExponentialFamily m)
    => Mean #> Natural # R k Poisson <* m -- ^ PPC
    -> Sample m -- ^ Sample Points
    -> [Double] -- ^ Sum of tuning curves at sample points
{-# INLINE sumOfTuningCurves #-}
sumOfTuningCurves lkl mus = S.sum . coordinates . dualTransition <$> lkl >$>* mus

-- | Returns the tuning curves of a population code over a set of sample points.
-- This is often useful for plotting purposes.
tuningCurves
    :: (ExponentialFamily m, KnownNat k)
    => Sample m -- Sample points
    -> Mean #> Natural # R k Poisson <* m -- ^ PPC
    -> [[(SamplePoint m, Double)]] -- ^ Vector of tuning curves
tuningCurves xsmps lkl =
    let tcs = L.transpose $ listCoordinates . dualTransition <$> (lkl >$>* xsmps)
     in zip xsmps <$> tcs

--- Internal ---

independentVariables0
    :: forall k m
    . (KnownNat k, ExponentialFamily m)
    => Mean #> Natural # R k Poisson <* m
    -> Sample m
    -> [S.Vector (Dimension m + 1) Double]
independentVariables0 _ mus =
    let sss :: [Mean # m]
        sss = sufficientStatistic <$> mus
     in (S.singleton 1 S.++) . coordinates <$> sss

independentVariables1
    :: (KnownNat k, ExponentialFamily m)
    => Mean #> Natural # R k Poisson <* m
    -> Sample m
    -> [S.Vector k Double]
independentVariables1 lkl mus =
    coordinates . dualTransition <$> lkl >$>* mus

normalBias :: Point Source Normal -> Double
normalBias sp =
    let [mu,vr] = listCoordinates sp
     in - mu^(2 :: Int)/(2*vr)



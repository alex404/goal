-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.PopulationCode
    ( -- * Population Encoders
      normalPopulationEncoder
    , vonMisesPopulationEncoder
    -- * Rectification
    , populationCodeRectificationParameters
    , rectifyPopulationCode
    , rectificationCurve
    -- * Utility
    , tuningCurves
    , sumOfTuningCurves
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions


--- Population Encoders ---

-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
--normalPopulationEncoder
--    :: S.Vector k (Point Source Normal) -- ^ Tuning Curves
--    -> Double -- ^ Gain
--    -> Point (Function Mean Natural) (Replicated k Poisson <* Normal) -- ^ Population Encoder

-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
normalPopulationEncoder
    :: KnownNat k
    => Bool -- ^ Normalize tuning curves
    -> Either Double (S.Vector k Double) -- ^ Global Gain or Gains
    -> S.Vector k (Point Source Normal) -- ^ Tuning Curves
    -> Point (Function Mean Natural) (Replicated k Poisson <* Normal) -- ^ Population Encoder
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
    -> Point (Function Mean Natural) (Replicated k Poisson <* VonMises) -- ^ Population Encoder
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
    :: (KnownNat k, ExponentialFamily m)
    => Double
    -> Natural # m
    -> Sample k m
    -> B.Vector k Double
{-# INLINE rectificationCurve #-}
rectificationCurve rho0 rprms mus = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- | Given a set of rectification parameters and a population code, returns the
-- population code which best satsisfies the resulting rectification equation.
rectifyPopulationCode
    :: (1 <= j, KnownNat k, KnownNat j, ExponentialFamily m)
    => Double
    -> Natural # m
    -> Sample j m
    -> Mean ~> Natural # R k Poisson <* m
    -> Mean ~> Natural # R k Poisson <* m
{-# INLINE rectifyPopulationCode #-}
rectifyPopulationCode rho0 rprms mus lkl =
    let dpnds = rectificationCurve rho0 rprms mus
        indpnds = independentVariables1 lkl mus
        gns = Point . S.map log $ S.linearLeastSquares indpnds (G.convert dpnds)
        (gns0,tcs) = splitAffine lkl
     in joinAffine (gns0 <+> gns) tcs

-- Linear Least Squares

-- | Returns the rectification parameters which best satisfy the rectification
-- equation for the given population code.
populationCodeRectificationParameters
    :: (KnownNat k, 1 <= j, KnownNat j, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> (Double, Natural # m)
{-# INLINE populationCodeRectificationParameters #-}
populationCodeRectificationParameters lkl mus =
    let dpnds = sumOfTuningCurves lkl mus
        indpnds = independentVariables0 lkl mus
        (rho0,rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

-- | The sum of the tuning curves of a population over a sample.
sumOfTuningCurves
    :: (KnownNat j, KnownNat k, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> S.Vector j Double
{-# INLINE sumOfTuningCurves #-}
sumOfTuningCurves lkl mus = mapReplicated (S.sum . coordinates . dualTransition) $ lkl >$>* mus

-- | Returns the tuning curves of a population code over a set of sample points.
-- This is often useful for plotting purposes.
tuningCurves
    :: (ExponentialFamily m, KnownNat k, KnownNat j)
    => Sample j m
    -> Mean ~> Natural # R k Poisson <* m
    -> B.Vector k (B.Vector j (SamplePoint m, Double))
tuningCurves xsmps lkl =
    let tcs = S.toRows . S.transpose . S.fromRows . mapReplicated (coordinates . dualTransition) $ lkl >$>* xsmps
     in B.zip xsmps . G.convert <$> G.convert tcs

--- Internal ---

independentVariables0
    :: forall j k m
    . (KnownNat k, KnownNat j, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> S.Vector j (S.Vector (Dimension m + 1) Double)
independentVariables0 _ mus =
    let sss :: B.Vector j (Mean # m)
        sss = sufficientStatistic <$> mus
     in G.convert $ ((S.singleton 1 S.++) . coordinates) <$> sss

independentVariables1
    :: (KnownNat k, KnownNat j, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> S.Vector j (S.Vector k Double)
independentVariables1 lkl mus =
    mapReplicated (coordinates . dualTransition) $ lkl >$>* mus

normalBias :: Point Source Normal -> Double
normalBias sp =
    let [mu,vr] = listCoordinates sp
     in - mu^(2 :: Int)/(2*vr)



-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.PopulationCode
    ( -- * Population Encoders
      normalPopulationEncoder
    , normalPopulationEncoder'
    , vonMisesPopulationEncoder
    , vonMisesPopulationEncoder'
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

independentVariables1
    :: (KnownNat k, KnownNat j, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> S.Vector j (S.Vector k Double)
independentVariables1 lkl mus =
    mapReplicated (coordinates . dualTransition) $ lkl >$>* mus

-- Linear Least Squares

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

sumOfTuningCurves
    :: (KnownNat j, KnownNat k, ExponentialFamily m)
    => Mean ~> Natural # R k Poisson <* m
    -> Sample j m
    -> S.Vector j Double
{-# INLINE sumOfTuningCurves #-}
sumOfTuningCurves lkl mus = mapReplicated (S.sum . coordinates . dualTransition) $ lkl >$>* mus

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

tuningCurves
    :: (ExponentialFamily m, KnownNat k, KnownNat j)
    => Sample j m
    -> Mean ~> Natural # R k Poisson <* m
    -> B.Vector k (B.Vector j (SamplePoint m, Double))
tuningCurves xsmps lkl =
    let tcs = S.toRows . S.transpose . S.fromRows . mapReplicated (coordinates . dualTransition) $ lkl >$>* xsmps
     in B.zip xsmps . G.convert <$> G.convert tcs

normalBias :: Point Source Normal -> Double
normalBias sp =
    let [mu,vr] = listCoordinates sp
     in - mu^(2 :: Int)/(2*vr)

-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
normalPopulationEncoder
    :: S.Vector k (Point Source Normal) -- ^ Tuning Curves
    -> Double -- ^ Gain
    -> Point (Function Mean Natural) (Replicated k Poisson <* Normal) -- ^ Population Encoder
normalPopulationEncoder sps gn =
    let mtx = S.concat $ S.map (coordinates . toNatural) sps
        ob = S.map ((+ log gn) . normalBias) sps
     in Point $ ob S.++ mtx

-- | Builds a linear population code, which is a population code that can be
-- expressed as an affine transformation across exponential family coordinate
-- systems.
normalPopulationEncoder'
    :: KnownNat k
    => S.Vector k (Point Source Normal) -- ^ Tuning Curves
    -> S.Vector k Double -- ^ Gains
    -> Point (Function Mean Natural) (Replicated k Poisson <* Normal) -- ^ Population Encoder
normalPopulationEncoder' sps gns =
    let mtx = S.concat $ S.map (coordinates . toNatural) sps
        ob = S.zipWith (+) (log gns) $ S.map normalBias sps
     in Point $ ob S.++ mtx


-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesPopulationEncoder
    :: KnownNat k
    => S.Vector k (Point Source VonMises) -- ^ Von Mises Curves
    -> Double -- ^ VM Gain
    -> Point (Function Mean Natural) (Replicated k Poisson <* VonMises) -- ^ Population Encoder
vonMisesPopulationEncoder sps gn =
    let mtx = S.concat $ S.map (coordinates . toNatural) sps
        ob = S.replicate $ log gn
     in Point $ ob S.++ mtx

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesPopulationEncoder'
    :: KnownNat k
    => S.Vector k (Point Source VonMises) -- ^ Von Mises Curves
    -> S.Vector k Double -- ^ Gains
    -> Point (Function Mean Natural) (Replicated k Poisson <* VonMises) -- ^ Population Encoder
vonMisesPopulationEncoder' sps gns =
    let mtx = S.concat $ S.map (coordinates . toNatural) sps
        ob = log gns
     in Point $ ob S.++ mtx


{-
-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesNormalPopulationEncoder'
    :: [Source :#: VonMises] -- ^ Von Mises Curves
    -> [Source :#: Normal] -- ^ Normal Tuning Curves
    -> Double -- ^ VM Gain
    -> Double -- ^ Normal Gain
    -> Function Mean Natural :#: Affine (Replicated Poisson) (VonMises, Normal) -- ^ Population Encoder
vonMisesNormalPopulationEncoder' vmtcs ntcs vmgn ngn =
    let (nob,nmtx) = splitAffine $ normalPopulationEncoder ntcs ngn
        vmmtx = fromCoordinates (Tensor (Replicated Poisson $ length vmtcs) VonMises) . C.concat
            $ coordinates . transitionTo Natural <$> vmtcs
        rpm = Replicated Poisson $ length vmtcs + length ntcs
        tnsm = Tensor rpm (VonMises, Normal)
        ob = fromCoordinates rpm $ C.replicate (length vmtcs) (log vmgn) C.++ coordinates nob
        tns = fromHMatrix tnsm . H.diagBlock $ [toHMatrix vmmtx, toHMatrix nmtx]
     in joinAffine ob tns


-- | Based on the non-tiled replicated normal population encoder.
multivariateNormalPopulationEncoder'
    :: [[Source :#: Normal]] -- ^ Tuning Curves
    -> [Double] -- ^ Gains
    -> Function Mean Natural :#: Affine (Replicated Poisson) MultivariateNormal -- ^ Population Encoder
multivariateNormalPopulationEncoder' spss gns =
    let nln = length spss
        pln = sum $ length <$> spss
        rp = Replicated Poisson pln
        tnsm = Tensor rp $ MultivariateNormal nln
        (ob,tns0) = splitAffine $ replicatedNormalPopulationEncoder' spss gns
        cls = H.toColumns $ toHMatrix tns0
        (mucls,sdcls0) = unzip [ (mucl,sdcl) | [mucl,sdcl] <- breakEvery 2 cls ]
        zro = C.replicate pln 0
        sdcls = concat [ replicate k zro ++ [sdcl] ++ replicate (nln - 1 - k) zro | (sdcl,k) <- zip sdcls0 [0..nln - 1] ]
        tns = fromHMatrix tnsm . H.fromColumns $ mucls ++ sdcls
     in joinAffine ob tns

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesNormalPopulationEncoder
    :: [Source :#: (VonMises, Normal)] -- ^ Tuning Curves
    -> Double -- ^ Gain
    -> Function Mean Natural :#: Affine (Replicated Poisson) (VonMises, Normal) -- ^ Population Encoder
vonMisesNormalPopulationEncoder sps gn =
    let nps = transitionTo Natural <$> sps
        rp = Replicated Poisson $ length nps
        ob = fromList rp $ (+ log gn) . normalBias . snd . splitPair' <$> sps
        tns = fromCoordinates (Tensor rp . manifold $ head sps) . C.concat $ coordinates <$> nps
     in joinAffine ob tns
     -}

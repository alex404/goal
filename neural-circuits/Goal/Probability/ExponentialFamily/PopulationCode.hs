-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.PopulationCode
    ( -- * Population Encoders
      normalPopulationEncoder
    , normalPopulationEncoder'
    , vonMisesPopulationEncoder
    , vonMisesPopulationEncoder'
    -- * Utility
    , tuningCurves
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

tuningCurves
    :: (ExponentialFamily l, KnownNat k, KnownNat j)
    => B.Vector j (SamplePoint l)
    -> Function Mean Natural # Replicated k Poisson <* l
    -> B.Vector k (B.Vector j (SamplePoint l, Double))
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

-- | Builds a population where the latent manifold is a 'Replicated'
-- 'Normal' 'Manifold'.
replicatedNormalPopulationEncoder
    :: [Source :#: Replicated Normal] -- ^ Tuning Curves
    -> Double -- ^ Gain
    -> Function Mean Natural :#: Affine (Replicated Poisson) (Replicated Normal) -- ^ Population Encoder
replicatedNormalPopulationEncoder sps gn =
    let nps = transitionTo Natural <$> sps
        rp = Replicated Poisson $ length nps
        ob = fromList rp $ (+ log gn) . sum . mapReplicated normalBias <$> sps
        tns = fromCoordinates (Tensor rp . manifold $ head sps) . C.concat $ coordinates <$> nps
     in joinAffine ob tns

-- | Non-tiled
replicatedNormalPopulationEncoder'
    :: [[Source :#: Normal]] -- ^ Tuning Curves
    -> [Double] -- ^ Gains
    -> Function Mean Natural :#: Affine (Replicated Poisson) (Replicated Normal) -- ^ Population Encoder
replicatedNormalPopulationEncoder' spss gns =
    let rp = Replicated Poisson . sum $ length <$> spss
        tnsm = Tensor rp . Replicated Normal $ length spss
        (obs,tnss) = unzip [ splitAffine $ normalPopulationEncoder sps gn | (gn,sps) <- zip gns spss ]
        ob = fromCoordinates rp . C.concat $ coordinates <$> obs
        tns = fromHMatrix tnsm . H.diagBlock $ toHMatrix <$> tnss
     in joinAffine ob tns

-- | Builds a population code where the latent manifold is a 'Replicated'
-- 'Manifold' of a 'VonMises' and 'Normal' pair. This results in a population
-- code for a 2-d dimensional stimulus with a rotational dimension, e.g. a
-- pendulum.
vonMisesPopulationEncoder'
    :: [Source :#: VonMises] -- ^ Von Mises Curves
    -> [Double] -- ^ VM Gain
    -> Function Mean Natural :#: Affine (Replicated Poisson) VonMises -- ^ Population Encoder
vonMisesPopulationEncoder' vmtcs vmgns =
    let vmmtx = fromCoordinates (Tensor (Replicated Poisson $ length vmtcs) VonMises) . C.concat
            $ coordinates . transitionTo Natural <$> vmtcs
        rpm = Replicated Poisson $ length vmtcs
        ob = fromList rpm (log <$> vmgns)
     in joinAffine ob vmmtx

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

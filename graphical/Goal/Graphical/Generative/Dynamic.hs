{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Generative.Dynamic
    (
    LatentProcess
    , sampleLatentProcess
    -- ** Construction
    , joinLatentProcess
    , splitLatentProcess
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Conditional
import Goal.Graphical.Inference
import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Generic ---


-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
data LatentProcess0 (f :: Type -> Type -> Type) (g :: Type -> Type -> Type) z x

type LatentProcess f g z x = LatentProcess0 f g [z] [x]

splitLatentProcess
    :: ( Manifold (f z x), Manifold (g x x), Manifold x, Manifold z )
    => c # LatentProcess f g z x
    -> (c # x, c # Affine f z x, c # Affine g x x)
splitLatentProcess ltnt =
    let (cx,cs') = S.splitAt $ coordinates ltnt
        (cf,cg) = S.splitAt cs'
     in (Point cx,Point cf,Point cg)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinLatentProcess
    :: ( Manifold (f z x), Manifold (g x x), Manifold x, Manifold z )
    => c # x
    -> c # Affine f z x
    -> c # Affine g x x
    -> c # LatentProcess f g z x -- ^ Conditional Harmonium
joinLatentProcess cx cf cg =
    Point $ coordinates cx S.++ coordinates cf S.++ coordinates cg

latentProcessTransition
    :: ( ConjugatedLikelihood f x x, ConjugatedLikelihood g z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear f x x
       , Generative Natural x, Generative Natural z
       , Bilinear g z x, Map Natural g x z )
    => Natural # Affine f x x
    -> Natural # Affine g z x
    -> SamplePoint x
    -> Random s (SamplePoint (z,x))
latentProcessTransition trns emsn x = do
    x' <- samplePoint $ trns >.>* x
    z' <- samplePoint $ emsn >.>* x'
    return (z',x')

sampleLatentProcess
    :: ( ConjugatedLikelihood g x x, ConjugatedLikelihood f z x
       , ExponentialFamily z, ExponentialFamily x, Bilinear g x x
       , Generative Natural x, Generative Natural z
       , Bilinear f z x, Map Natural f x z )
    => Int
    -> Natural # LatentProcess f g z x
    -> Random s (Sample (z,x))
sampleLatentProcess n ltnt = do
    let (prr,emsn,trns) = splitLatentProcess ltnt
    x0 <- samplePoint prr
    z0 <- samplePoint $ emsn >.>* x0
    iterateM (n-1) (latentProcessTransition trns emsn . snd) (z0,x0)


--- Instances ---

-- Implementations

latentProcessLogDensity
    :: ( ExponentialFamily z, ExponentialFamily x, Map Natural f z x
       , Map Natural g x x, AbsolutelyContinuous Natural x
       , AbsolutelyContinuous Natural z  )
    => Natural # x
    -> Natural # Affine f z x
    -> Natural # Affine g x x
    -> Sample (z,x)
    -> Double
latentProcessLogDensity prr emsn trns zxs =
    let (zs,xs) = unzip zxs
        prrdns = logDensity prr $ head xs
        trnsdnss = zipWith logDensity (trns >$>* xs) $ tail xs
        emsndnss = zipWith logDensity (emsn >$>* xs) zs
     in sum $ prrdns : trnsdnss ++ emsndnss

conjugatedSmoothingLogDensity
    :: ( ExponentialFamily z, ExponentialFamily x, ConjugatedLikelihood g x x
       , ConjugatedLikelihood f z x, Bilinear f z x, Bilinear g x x
       , Map Natural f x z, ObservablyContinuous Natural (Harmonium f) z x )
    => Natural # Affine f z x
    -> Natural # Affine g x x
    -> Natural # x
    -> Sample z
    -> Double
conjugatedSmoothingLogDensity emsn trns prr zs =
    let flts = fst $ conjugatedSmoothing trns emsn prr zs
        hrms = joinConjugatedHarmonium emsn <$> flts
     in sum $ zipWith logObservableDensity hrms zs

-- Latent Processes

instance ( Manifold (f z x), Manifold (g x x), Manifold x, Manifold z )
  => Manifold (LatentProcess f g z x) where
      type Dimension (LatentProcess f g z x)
        = Dimension (Affine f z x) + Dimension (Affine g x x) + Dimension x

instance Manifold (LatentProcess f g z x) => Statistical (LatentProcess f g z x) where
    type SamplePoint (LatentProcess f g z x) = [SamplePoint (z,x)]

instance ( ExponentialFamily z, ExponentialFamily x, Map Natural f z x
         , Map Natural g x x, AbsolutelyContinuous Natural x
         , AbsolutelyContinuous Natural z  )
  => AbsolutelyContinuous Natural (LatentProcess f g z x) where
    logDensity ltnt zxs =
        let (prr,emsn,trns) = splitLatentProcess ltnt
         in latentProcessLogDensity prr emsn trns zxs

instance ( ExponentialFamily z, ExponentialFamily x, ConjugatedLikelihood g x x
         , ConjugatedLikelihood f z x, Bilinear f z x, Bilinear g x x
         , Map Natural f x z, ObservablyContinuous Natural (Harmonium f) z x)
  => ObservablyContinuous Natural (LatentProcess0 f g) [z] [x] where
    logObservableDensity ltnt zs =
        let (prr,emsn,trns) = splitLatentProcess ltnt
         in conjugatedSmoothingLogDensity emsn trns prr zs

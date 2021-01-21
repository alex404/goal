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

import Goal.Graphical.Inference
import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium


--- Generic ---


-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
newtype LatentProcess f g y x z w
    = LatentProcess (Harmonium f y x z w, Affine g x w x)

type instance Observation (LatentProcess f g y x z w) = Sample z

deriving instance (Manifold (Harmonium f y x z w), Manifold (Affine g x w x))
  => Manifold (LatentProcess f g y x z w)
deriving instance (Manifold (Harmonium f y x z w), Manifold (Affine g x w x))
  => Product (LatentProcess f g y x z w)

splitLatentProcess
    :: (Manifold z, Manifold w, Manifold (f y x), Manifold (g x x))
    => c # LatentProcess f g y x z w
    -> (c # w, c # Affine f y z x, c # Affine g x w x)
splitLatentProcess ltnt =
    let (hrm,trns) = split ltnt
        (emsn,prr) = split hrm
     in (prr,emsn,trns)

joinLatentProcess
    :: (Manifold z, Manifold w, Manifold (f y x), Manifold (g x x))
    => c # w
    -> c # Affine f y z x
    -> c # Affine g x w x
    -> c # LatentProcess f g y x z w
joinLatentProcess prr emsn trns =
    let hrm = join emsn prr
     in join hrm trns

latentProcessTransition
    :: ( SamplePoint w ~ SamplePoint x, ExponentialFamily z
       , Translation w x, Translation z y, Map Natural g x x
       , ExponentialFamily x, Bilinear f x x
       , Generative Natural w, Generative Natural z
       , Bilinear g z x, Map Natural f y x )
    => Natural # Affine g x w x -- ^ Transition Distribution
    -> Natural # Affine f y z x -- ^ Emission Distribution
    -> SamplePoint w
    -> Random r (SamplePoint (z,w))
latentProcessTransition trns emsn w = do
    w' <- samplePoint $ trns >.>* w
    z' <- samplePoint $ emsn >.>* w'
    return (z',w')

sampleLatentProcess
    :: ( SamplePoint w ~ SamplePoint x, ExponentialFamily z
       , Translation w x, Translation z y, Map Natural g x x
       , ExponentialFamily x, Bilinear f x x
       , Generative Natural w, Generative Natural z
       , Bilinear g z x, Map Natural f y x )
    => Int
    -> Natural # LatentProcess f g y x z w
    -> Random s (Sample (z,x))
sampleLatentProcess n ltnt = do
    let (prr,emsn,trns) = splitLatentProcess ltnt
    x0 <- samplePoint prr
    z0 <- samplePoint $ emsn >.>* x0
    iterateM (n-1) (latentProcessTransition trns emsn . snd) (z0,x0)


--- Instances ---

-- Implementations

latentProcessLogDensity
    :: ( ExponentialFamily z, ExponentialFamily x, Map Natural f y x
       , Translation z y , Map Natural g x x, AbsolutelyContinuous Natural w
       , SamplePoint w ~ SamplePoint x, AbsolutelyContinuous Natural z, Translation w x )
    => Natural # w
    -> Natural # Affine f y z x -- ^ Emission Distribution
    -> Natural # Affine g x w x -- ^ Transition Distribution
    -> Sample (z,w)
    -> Double
latentProcessLogDensity prr emsn trns zxs =
    let (zs,xs) = unzip zxs
        prrdns = logDensity prr $ head xs
        trnsdnss = zipWith logDensity (trns >$>* xs) $ tail xs
        emsndnss = zipWith logDensity (emsn >$>* xs) zs
     in sum $ prrdns : trnsdnss ++ emsndnss

conjugatedSmoothingLogDensity
    :: ( ConjugatedLikelihood g x x w w, Bilinear g x x
       , ConjugatedLikelihood f y x z w, Bilinear f y x
       , Map Natural g x x, Map Natural f x y
       , LegendreExponentialFamily z, LegendreExponentialFamily w )
    => Natural # w
    -> Natural # Affine f y z x -- ^ Emission Distribution
    -> Natural # Affine g x w x -- ^ Transition Distribution
    -> Sample z
    -> Double
conjugatedSmoothingLogDensity prr emsn trns ws =
    let smths = fst $ conjugatedSmoothing trns emsn prr ws
        hrms = joinConjugatedHarmonium emsn <$> smths
     in sum $ zipWith logObservableDensity hrms ws

-- Latent Processes

instance Manifold (LatentProcess f g y x z w) => Statistical (LatentProcess f g y x z w) where
    type SamplePoint (LatentProcess f g y x z w) = [SamplePoint (z,x)]

instance ( ExponentialFamily z, ExponentialFamily x, Map Natural f y x
         , Translation z y , Map Natural g x x, AbsolutelyContinuous Natural w
         , SamplePoint w ~ SamplePoint x, AbsolutelyContinuous Natural z, Translation w x )
  => AbsolutelyContinuous Natural (LatentProcess f g y x z w) where
    logDensity ltnt zxs =
        let (prr,emsn,trns) = splitLatentProcess ltnt
         in latentProcessLogDensity prr emsn trns zxs

instance ( ConjugatedLikelihood g x x w w, Bilinear g x x
         , ConjugatedLikelihood f y x z w, Bilinear f y x
         , Map Natural g x x, Map Natural f x y
         , LegendreExponentialFamily z, LegendreExponentialFamily w )
  => ObservablyContinuous Natural (LatentProcess f g y x z w) where
    logObservableDensity ltnt zs =
        let (prr,emsn,trns) = splitLatentProcess ltnt
         in conjugatedSmoothingLogDensity prr emsn trns zs

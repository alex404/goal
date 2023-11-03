{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Models.Dynamic (
    LatentProcess (LatentProcess),
    HiddenMarkovModel,
    KalmanFilter,
    sampleLatentProcess,

    -- ** Construction
    joinLatentProcess,
    splitLatentProcess,

    -- ** Inference
    conjugatedFiltering,
    conjugatedSmoothing,
    conjugatedSmoothing0,
) where

--- Imports  ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Core.Vector.Storable.Linear as L

import Goal.Graphical.Inference
import Goal.Graphical.Models
import Goal.Graphical.Models.Harmonium

--- Misc

import Data.List

--- Generic ---

{- | A conditional 'Harmonium', where the observable biases of the
'Harmonium' model depend on additional variables.
-}
newtype LatentProcess f g y x z w
    = LatentProcess (AffineHarmonium f y x z w, Affine g x w x)

type HiddenMarkovModel n k =
    LatentProcess L.Full L.Full (Categorical n) (Categorical n) (Categorical k) (Categorical k)

type KalmanFilter n k =
    LatentProcess L.Full L.Full (StandardNormal n) (StandardNormal k) (FullNormal n) (FullNormal k)

type instance Observation (LatentProcess f g y x z w) = Sample z

deriving instance
    (Manifold (AffineHarmonium f y x z w), Manifold (Affine g x w x)) =>
    Manifold (LatentProcess f g y x z w)
deriving instance
    (Manifold (AffineHarmonium f y x z w), Manifold (Affine g x w x)) =>
    Product (LatentProcess f g y x z w)

{- | Split a 'LatentProcess' into a prior, an emission distribution, and a
transition distribution.
-}
splitLatentProcess ::
    (KnownLinear f y x, KnownLinear g x x, Manifold z, Manifold w) =>
    c # LatentProcess f g y x z w ->
    (c # w, c # Affine f y z x, c # Affine g x w x)
splitLatentProcess ltnt =
    let (hrm, trns) = split ltnt
        (emsn, prr) = split hrm
     in (prr, emsn, trns)

{- | Construct a 'LatentProcess' from a prior, an emission distribution, and a
transition distribution.
-}
joinLatentProcess ::
    (KnownLinear f y x, KnownLinear g x x, Manifold z, Manifold w) =>
    c # w ->
    c # Affine f y z x ->
    c # Affine g x w x ->
    c # LatentProcess f g y x z w
joinLatentProcess prr emsn trns =
    let hrm = join emsn prr
     in join hrm trns

latentProcessTransition ::
    ( SamplePoint w ~ SamplePoint x
    , ExponentialFamily z
    , LinearSubspace w x
    , LinearSubspace z y
    , KnownLinear f y x
    , KnownLinear g x x
    , ExponentialFamily x
    , Generative Natural w
    , Generative Natural z
    ) =>
    -- | Emission Distribution
    Natural # Affine f y z x ->
    -- | Transition Distribution
    Natural # Affine g x w x ->
    SamplePoint w ->
    Random (SamplePoint (z, w))
latentProcessTransition emsn trns w = do
    w' <- samplePoint $ trns >.>* w
    z' <- samplePoint $ emsn >.>* w'
    return (z', w')

{- | Generate a realization of the observable and latent states from a given
latent process.
-}
sampleLatentProcess ::
    ( SamplePoint w ~ SamplePoint x
    , ExponentialFamily z
    , LinearSubspace w x
    , LinearSubspace z y
    , KnownLinear f y x
    , KnownLinear g x x
    , ExponentialFamily x
    , Generative Natural w
    , Generative Natural z
    ) =>
    Int ->
    Natural # LatentProcess f g y x z w ->
    Random (Sample (z, x))
sampleLatentProcess n ltnt = do
    let (prr, emsn, trns) = splitLatentProcess ltnt
    x0 <- samplePoint prr
    z0 <- samplePoint $ emsn >.>* x0
    iterateM (n - 1) (latentProcessTransition emsn trns . snd) (z0, x0)

-- | Filtering for latent processes based on conjugated distributions.
conjugatedFiltering ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    ) =>
    Natural # LatentProcess f g y x z w ->
    Sample z ->
    [Natural # w]
conjugatedFiltering _ [] = []
conjugatedFiltering ltnt (z : zs') =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        prr' = conjugatedBayesRule emsn prr z
     in scanl' (conjugatedForwardStep trns emsn) prr' zs'

-- | Smoothing for latent processes based on conjugated distributions.
conjugatedSmoothing ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    ) =>
    Natural # LatentProcess f g y x z w ->
    Sample z ->
    [Natural # w]
conjugatedSmoothing ltnt zs =
    let (prr, emsn, trns) = splitLatentProcess ltnt
     in fst $ conjugatedSmoothing0 prr emsn trns zs

{- | A more low-level implementation of smoothing which also returns joint
distributions over current and subsequent states.
-}
conjugatedSmoothing0 ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    ) =>
    Natural # w ->
    -- | Emission Distribution
    Natural # Affine f y z x ->
    -- | Transition Distribution
    Natural # Affine g x w x ->
    Sample z ->
    ([Natural # w], [Natural # AffineHarmonium g x x w w])
conjugatedSmoothing0 _ _ _ [] = ([], [])
conjugatedSmoothing0 prr emsn _ [z] =
    ([conjugatedBayesRule emsn prr z], [])
conjugatedSmoothing0 prr emsn trns (z : zs) =
    let pst = conjugatedBayesRule emsn prr z
        (trns', fwd) =
            splitConjugatedHarmonium
                . transposeHarmonium
                $ joinConjugatedHarmonium trns pst
        (smth : smths, hrms) = conjugatedSmoothing0 fwd emsn trns zs
        hrm = transposeHarmonium $ joinConjugatedHarmonium trns' smth
        bwd = snd $ splitConjugatedHarmonium hrm
     in (bwd : smth : smths, hrm : hrms)

--- Instances ---

-- Implementations

latentProcessLogDensity ::
    ( ExponentialFamily z
    , ExponentialFamily x
    , KnownLinear f y x
    , KnownLinear g x x
    , LinearSubspace z y
    , AbsolutelyContinuous Natural w
    , SamplePoint w ~ SamplePoint x
    , AbsolutelyContinuous Natural z
    , LinearSubspace w x
    ) =>
    Natural # w ->
    -- | Emission Distribution
    Natural # Affine f y z x ->
    -- | Transition Distribution
    Natural # Affine g x w x ->
    Sample (z, w) ->
    Double
latentProcessLogDensity prr emsn trns zxs =
    let (zs, xs) = unzip zxs
        prrdns = logDensity prr $ head xs
        trnsdnss = zipWith logDensity (trns >$>* xs) $ tail xs
        emsndnss = zipWith logDensity (emsn >$>* xs) zs
     in sum $ prrdns : trnsdnss ++ emsndnss

latentProcessMarginalLogDensity ::
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    , ExponentialFamily y
    , LegendreExponentialFamily z
    , LegendreExponentialFamily w
    ) =>
    Natural # LatentProcess f g y x z w ->
    Sample z ->
    Double
latentProcessMarginalLogDensity ltnt zs =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        prrs =
            iterate
                (snd . splitConjugatedHarmonium . transposeHarmonium . joinConjugatedHarmonium trns)
                prr
        hrms = joinConjugatedHarmonium emsn <$> prrs
     in sum $ zipWith logObservableDensity hrms zs

-- Latent Processes

instance (Manifold (LatentProcess f g y x z w)) => Statistical (LatentProcess f g y x z w) where
    type SamplePoint (LatentProcess f g y x z w) = [SamplePoint (z, x)]

instance
    ( ExponentialFamily z
    , ExponentialFamily x
    , LinearSubspace z y
    , KnownLinear f y x
    , KnownLinear g x x
    , AbsolutelyContinuous Natural w
    , SamplePoint w ~ SamplePoint x
    , AbsolutelyContinuous Natural z
    , LinearSubspace w x
    ) =>
    AbsolutelyContinuous Natural (LatentProcess f g y x z w)
    where
    logDensities ltnt zxss = do
        zxs <- zxss
        let (prr, emsn, trns) = splitLatentProcess ltnt
        return $ latentProcessLogDensity prr emsn trns zxs

instance
    ( ConjugatedLikelihood g x x w w
    , ConjugatedLikelihood f y x z w
    , ExponentialFamily y
    , LegendreExponentialFamily z
    , LegendreExponentialFamily w
    ) =>
    ObservablyContinuous Natural (LatentProcess f g y x z w)
    where
    logObservableDensities ltnt = map (latentProcessMarginalLogDensity ltnt)

-- instance ( Manifold w , Manifold (g x x)
--         , LinearSubspace z y, Bilinear c f x y, Bilinear Natural f y x )
--  => LinearSubspace (LatentProcess f g y x z w) y where
--    (>+>) ltnt y =
--        let (ehrm,trns) = split ltnt
--            (z,yx,w) = splitHarmonium ehrm
--            z' = z >+> y
--         in join (joinHarmonium z' yx w) trns
--    projection ltnt =
--        projection . snd . split . transposeHarmonium . fst $ split ltnt

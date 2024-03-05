{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Models.Dynamic (
    LatentProcess (LatentProcess),
    KnownLatentProcess,
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
newtype LatentProcess f g x0 z0 x z
    = LatentProcess (AffineHarmonium f x0 z0 x z, Affine g z0 z z0)

type KnownLatentProcess f g x0 z0 x z = (KnownAffineHarmonium f x0 z0 x z, KnownLinear g z0 z0)

type HiddenMarkovModel n k =
    LatentProcess L.Full L.Full (Categorical n) (Categorical n) (Categorical k) (Categorical k)

type KalmanFilter n k =
    LatentProcess L.Full L.Full (StandardNormal n) (StandardNormal k) (FullNormal n) (FullNormal k)

type instance Observation (LatentProcess f g x0 z0 x z) = Sample x

deriving instance
    (Manifold (AffineHarmonium f x0 z0 x z), Manifold (Affine g z0 z z0)) =>
    Manifold (LatentProcess f g x0 z0 x z)
deriving instance
    (Manifold (AffineHarmonium f x0 z0 x z), Manifold (Affine g z0 z z0)) =>
    Product (LatentProcess f g x0 z0 x z)

{- | Split a 'LatentProcess' into a prior, an emission distribution, and a
transition distribution.
-}
splitLatentProcess ::
    (KnownLinear f x0 z0, KnownLinear g z0 z0, Manifold z, Manifold x) =>
    c # LatentProcess f g x0 z0 x z ->
    (c # z, c # Affine f x0 x z0, c # Affine g z0 z z0)
splitLatentProcess ltnt =
    let (hrm, trns) = split ltnt
        (emsn, prr) = split hrm
     in (prr, emsn, trns)

{- | Construct a 'LatentProcess' from a prior, an emission distribution, and a
transition distribution.
-}
joinLatentProcess ::
    (KnownLinear f x0 z0, KnownLinear g z0 z0, Manifold z, Manifold x) =>
    c # z ->
    c # Affine f x0 x z0 ->
    c # Affine g z0 z z0 ->
    c # LatentProcess f g x0 z0 x z
joinLatentProcess prr emsn trns =
    let hrm = join emsn prr
     in join hrm trns

{- | Generate a realization of the observable and latent states from a given
latent process.
-}
sampleLatentProcess ::
    forall f g x0 z0 x z.
    ( KnownLatentProcess f g x0 z0 x z
    , Generative Natural z
    , Generative Natural x
    ) =>
    Int ->
    Natural # LatentProcess f g x0 z0 x z ->
    Random (Sample (x, z))
sampleLatentProcess n ltnt = do
    let (prr, emsn, trns) = splitLatentProcess ltnt
    z0 <- samplePoint prr
    let mz :: Mean # z
        mz = sufficientStatistic z0
    x0 <- samplePoint $ emsn >.> linearProjection mz
    iterateM (n - 1) (latentProcessTransition emsn trns . snd) (x0, z0)

-- | Filtering for latent processes based on conjugated distributions.
conjugatedFiltering ::
    ( ConjugatedLikelihood g z0 z0 z z
    , ConjugatedLikelihood f x0 z0 x z
    ) =>
    Natural # LatentProcess f g x0 z0 x z ->
    Sample x ->
    [Natural # z]
conjugatedFiltering _ [] = []
conjugatedFiltering ltnt (x : xs') =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        prr' = conjugatedBayesRule emsn prr x
     in scanl' (conjugatedForwardStep trns emsn) prr' xs'

-- | Smoothing for latent processes based on conjugated distributions.
conjugatedSmoothing ::
    ( ConjugatedLikelihood g z0 z0 z z
    , ConjugatedLikelihood f x0 z0 x z
    ) =>
    Natural # LatentProcess f g x0 z0 x z ->
    Sample x ->
    [Natural # z]
conjugatedSmoothing ltnt xs =
    let (prr, emsn, trns) = splitLatentProcess ltnt
     in fst $ conjugatedSmoothing0 prr emsn trns xs

{- | A more low-level implementation of smoothing which also returns joint
distributions over current and subsequent states.
-}
conjugatedSmoothing0 ::
    ( ConjugatedLikelihood f x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , KnownLatentProcess f g x0 z0 x z
    ) =>
    Natural # z ->
    -- | Emission Distribution
    Natural # Affine f x0 x z0 ->
    -- | Transition Distribution
    Natural # Affine g z0 z z0 ->
    Sample x ->
    ([Natural # z], [Natural # AffineHarmonium g z0 z0 z z])
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

--- Internal ---

latentProcessMarginalLogDensity ::
    ( ConjugatedLikelihood f x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , LegendreExponentialFamily z
    ) =>
    Natural # LatentProcess f g x0 z0 x z ->
    Sample x ->
    Double
latentProcessMarginalLogDensity ltnt xs =
    let (prr, emsn, trns) = splitLatentProcess ltnt
        prrs =
            iterate
                (snd . splitConjugatedHarmonium . transposeHarmonium . joinConjugatedHarmonium trns)
                prr
        hrms = joinConjugatedHarmonium emsn <$> prrs
     in sum $ zipWith logObservableDensity hrms xs

latentProcessTransition ::
    forall f g x0 z0 x z.
    ( KnownLatentProcess f g x0 z0 x z
    , Generative Natural x
    , Generative Natural z
    ) =>
    -- | Emission Distribution
    Natural # Affine f x0 x z0 ->
    -- | Transition Distribution
    Natural # Affine g z0 z z0 ->
    SamplePoint z ->
    Random (SamplePoint (x, z))
latentProcessTransition emsn trns z = do
    let mz :: Mean # z
        mz = sufficientStatistic z
        mz0 = linearProjection mz
    z' <- samplePoint $ trns >.> mz0
    x' <- samplePoint $ emsn >.> mz0
    return (x', z')

latentProcessLogDensity ::
    forall f g x0 z0 x z.
    ( KnownLatentProcess f g x0 z0 x z
    , AbsolutelyContinuous Natural z
    , AbsolutelyContinuous Natural x
    ) =>
    Natural # z ->
    -- | Emission Distribution
    Natural # Affine f x0 x z0 ->
    -- | Transition Distribution
    Natural # Affine g z0 z z0 ->
    Sample (x, z) ->
    Double
latentProcessLogDensity prr emsn trns xzs =
    let (xs, zs) = unzip xzs
        mzs :: [Mean # z]
        mzs = sufficientStatistic <$> zs
        mzs0 = linearProjection <$> mzs
        prrdns = logDensity prr $ head zs
        trnsdnss = zipWith logDensity (trns >$> mzs0) $ tail zs
        emsndnss = zipWith logDensity (emsn >$> mzs0) xs
     in sum $ prrdns : trnsdnss ++ emsndnss

--- Instances ---

instance (Manifold (LatentProcess f g x0 z0 x z)) => Statistical (LatentProcess f g x0 z0 x z) where
    type SamplePoint (LatentProcess f g x0 z0 x z) = [SamplePoint (x, z)]

instance
    ( KnownLatentProcess f g x0 z0 x z
    , AbsolutelyContinuous Natural z
    , AbsolutelyContinuous Natural x
    ) =>
    AbsolutelyContinuous Natural (LatentProcess f g x0 z0 x z)
    where
    logDensities ltnt xzss = do
        xzs <- xzss
        let (prr, emsn, trns) = splitLatentProcess ltnt
        return $ latentProcessLogDensity prr emsn trns xzs

instance
    ( ConjugatedLikelihood f x0 z0 x z
    , ConjugatedLikelihood g z0 z0 z z
    , LegendreExponentialFamily z
    ) =>
    ObservablyContinuous Natural (LatentProcess f g x0 z0 x z)
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
--    linearProjection ltnt =
--        linearProjection . snd . split . transposeHarmonium . fst $ split ltnt

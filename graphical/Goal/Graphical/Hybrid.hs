{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Hybrid
    (
     -- * Conditional Harmoniums
      ConditionalHarmonium
    , ConditionalMixture
    , LatentProcess
    -- ** Construction
    , joinConditionalHarmonium
    , splitConditionalHarmonium
    , joinLatentProcess
    , splitLatentProcess
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Conditional
import Goal.Graphical.Generative
import Goal.Graphical.Generative.Harmonium
import Goal.Graphical.Inference

import qualified Goal.Core.Vector.Storable as S


--- Generic ---



--- Types ---


-- | A generic conditional model.
data ConditionalBias (f :: Type -> Type -> Type) z y

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium g f z x = ConditionalBias g (Harmonium f z x)

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture f z k = ConditionalHarmonium f Tensor z (Categorical k) -- ^ Function

-- | Splits a conditional 'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalHarmonium
    :: ( Map Natural f z x, Manifold (g z y), Manifold x )
    => c # ConditionalHarmonium g f z x y -- ^ Conditional Harmonium
    -> (c # g z y, c # f z x, c # x) -- ^ Matrix function and upper part
splitConditionalHarmonium dhrm =
    let (cg,cs') = S.splitAt $ coordinates dhrm
        (cf,cx) = S.splitAt cs'
     in (Point cg,Point cf,Point cx)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalHarmonium
    :: ( Map Natural f z x, Manifold (g z y), Manifold x )
    => c # g z y
    -> c # f z x
    -> c # x
    -> c # ConditionalHarmonium g f z x y -- ^ Conditional Harmonium
joinConditionalHarmonium gzy fzx x =
    Point $ coordinates gzy S.++ coordinates fzx S.++ coordinates x



--- Latent Process ---

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

latentProcessMeanParameters
    :: ( Transition Natural Mean x, Transition Natural Mean (Harmonium f z x)
       , ConjugatedLikelihood f z x, Transition Natural Mean (Harmonium g x x)
       , ConjugatedLikelihood g x x )
    => Natural # LatentProcess f g z x -- ^ Conditional Harmonium
    -> Mean # LatentProcess f g z x -- ^ Conditional Harmonium
latentProcessMeanParameters ltnt =
    let (nprr,nemsn,ntrns) = splitLatentProcess ltnt
        mprr = toMean nprr
        mesn = fst . splitBottomHarmonium . toMean $ joinConjugatedHarmonium nemsn nprr
        mtrns = fst . splitBottomHarmonium . toMean $ joinConjugatedHarmonium ntrns nprr
     in joinLatentProcess mprr mesn mtrns

latentProcessNaturalParameters
    :: ( Transition Mean Natural x, Transition Mean Natural (Harmonium f z x)
       , ConjugatedLikelihood f z x, Transition Mean Natural (Harmonium g x x)
       , ConjugatedLikelihood g x x )
    => Mean # LatentProcess f g z x -- ^ Conditional Harmonium
    -> Natural # LatentProcess f g z x -- ^ Conditional Harmonium
latentProcessNaturalParameters ltnt =
    let (mprr,memsn,mtrns) = splitLatentProcess ltnt
        nprr = toNatural mprr
        nemsn = fst . splitConjugatedHarmonium . toNatural $ joinBottomHarmonium memsn mprr
        ntrns = fst . splitConjugatedHarmonium . toNatural $ joinBottomHarmonium mtrns mprr
     in joinLatentProcess nprr nemsn ntrns

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

latentProcessExpectationStep
    :: ( ConjugatedLikelihood f z x, ConjugatedLikelihood g x x
       , Propagate Natural f z x, Propagate Natural g x x
       , ExponentialFamily z, DuallyFlatExponentialFamily x, Bilinear f z x
       , Bilinear g x x, Map Natural f x z
       , DuallyFlatExponentialFamily (Harmonium f z x)
       , DuallyFlatExponentialFamily (Harmonium g x x) )
    => Natural # x
    -> Natural # Affine f z x
    -> Natural # Affine g x x
    -> [Sample z]
    -> (Mean # x, Mean # Affine f z x, Mean # Affine g x x)
latentProcessExpectationStep prr emsn trns zss =
    let (smthss,hrmss) = unzip $ conjugatedSmoothing trns emsn prr <$> zss
        mtrns = fst . splitBottomHarmonium . average
            $ toMean <$> concat hrmss
        msmths = toMean <$> concat smthss
        mprr = average $ toMean . head <$> smthss
        mzs = concat $ map sufficientStatistic <$> zss
        mehrm = joinHarmonium (average mzs) (mzs >$< msmths) (average msmths)
        memsn = fst $ splitBottomHarmonium mehrm
     in (mprr,memsn,mtrns)


--- Instances ---


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

instance ( Transition Natural Mean x, Transition Natural Mean (Harmonium f z x)
       , ConjugatedLikelihood f z x, Transition Natural Mean (Harmonium g x x)
       , ConjugatedLikelihood g x x )
       => Transition Natural Mean (LatentProcess f g z x) where
           transition = latentProcessMeanParameters

instance ( Transition Mean Natural x, Transition Mean Natural (Harmonium f z x)
         , ConjugatedLikelihood f z x, Transition Mean Natural (Harmonium g x x)
         , ConjugatedLikelihood g x x )
         => Transition Mean Natural (LatentProcess f g z x) where
             transition = latentProcessNaturalParameters

instance ( ConjugatedLikelihood f z x, ConjugatedLikelihood g x x
         , Propagate Natural f z x, Propagate Natural g x x
         , ExponentialFamily z, DuallyFlatExponentialFamily x, Bilinear f z x
         , Bilinear g x x, Map Natural f x z
         , DuallyFlatExponentialFamily (Harmonium f z x)
         , DuallyFlatExponentialFamily (Harmonium g x x) )
         => ExpectationMaximization (LatentProcess0 f g) [z] [x] where
             expectationStep zss ltnt =
                 let (nprr,nemsn,ntrns) = splitLatentProcess ltnt
                     (mprr,memsn,mtrns) = latentProcessExpectationStep nprr nemsn ntrns zss
                  in joinLatentProcess mprr memsn mtrns

-- Conditional Harmonium

instance ( Map Natural g z y, Manifold (f z x), Manifold x )
  => Manifold (ConditionalHarmonium g f z x y) where
      type Dimension (ConditionalHarmonium g f z x y)
        = Dimension (g z y) + Dimension (f z x) + Dimension x

instance ( Map Natural g z y, Map Natural f z x, Manifold (f z x), Manifold x )
     => Map Natural (ConditionalBias g) (Harmonium f z x) y where
    (>.>) chrm my =
        let (gzy,fzx,nx) = splitConditionalHarmonium chrm
            affzx = joinAffine (gzy >.> my) fzx
         in joinBottomHarmonium affzx nx
    (>$>) chrm mys =
        let (gzy,fzx,nx) = splitConditionalHarmonium chrm
            affzxs = flip joinAffine fzx <$> (gzy >$> mys)
         in flip joinBottomHarmonium nx <$> affzxs

instance ( Propagate Natural g z y, Map Natural f z x, Manifold (f z x), Manifold x )
     => Propagate Natural (ConditionalBias g) (Harmonium f z x) y where
        propagate mhrms mys chrm =
            let (mzs,mfzxs,mxs) = unzip3 $ splitHarmonium <$> mhrms
                (ngzy,nfzx,nx) = splitConditionalHarmonium chrm
                (mgzx,nzhts) = propagate mzs mys ngzy
             in ( joinConditionalHarmonium mgzx (average mfzxs) (average mxs)
                , [joinHarmonium nzht nfzx nx | nzht <- nzhts] )

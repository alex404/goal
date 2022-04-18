{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE
    TypeApplications,
    UndecidableInstances,
    NoStarIsType,
    GeneralizedNewtypeDeriving,
    StandaloneDeriving,
    ScopedTypeVariables,
    ExplicitNamespaces,
    TypeOperators,
    KindSignatures,
    DataKinds,
    RankNTypes,
    TypeFamilies,
    FlexibleContexts,
    MultiParamTypeClasses,
    ConstraintKinds,
    FlexibleInstances
#-}
-- | An Exponential Family 'Harmonium' is a product exponential family with a
-- particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
-- A 'Mixture' model is a special case of harmonium. A 'FactorAnalysis' model
-- can also be interpreted as a 'Harmonium' with a fixed latent distribution.
module Goal.Graphical.Models.Harmonium.Gaussian
    (
    -- * Factor Analysis
      linearModelObservableDistribution
    , linearModelExpectationMaximization
    , factorAnalysisUniqueness
    ) where


--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models
import Goal.Graphical.Models.Harmonium
import Goal.Graphical.Learning

import qualified Goal.Core.Vector.Storable as S


--- Factor Analysis ---


type instance Observation (FactorAnalysis n k) = S.Vector n Double

linearModelObservableDistribution
    :: forall f n k
     . ( KnownNat n, KnownNat k, Square Natural f (MVNMean n)
       , ExponentialFamily (MultivariateNormal f n)
       , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k) )
    => Natural # LinearModel f n k
    -> Natural # FullNormal n
linearModelObservableDistribution lm =
    let lgh :: Natural # LinearGaussianHarmonium f n k
        lgh = joinConjugatedHarmonium lm standardNormal
        (lkl,nz) = split lgh
        (nx,aff) = split lkl
        (nmux,nvrx) = split nx
        nx' = join nmux . fromTensor $ toTensor nvrx
        lkl' = join nx' aff
     in snd . splitConjugatedHarmonium . transposeHarmonium $ join lkl' nz

linearModelExpectationMaximization
    :: forall f n k
     . ( Square Natural f (MVNMean n), ExponentialFamily (MultivariateNormal f n)
       , KnownNat k, Square Source f (MVNMean n)
       , Transition Source Mean (MultivariateNormal f n)
       , Transition Mean Source (MultivariateNormal f n)
       , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k)
       , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean n)
       , LinearlyComposable Tensor f (MVNMean n) (MVNMean n) (MVNMean n) )
    => [S.Vector n Double]
    -> Natural # LinearModel f n k
    -> Natural # LinearModel f n k
linearModelExpectationMaximization xs lm =
    let lgh :: Natural # LinearGaussianHarmonium f n k
        lgh =  expectationMaximization xs $ joinConjugatedHarmonium lm standardNormal
     in fst $ split lgh

factorAnalysisUniqueness
    :: forall n k . (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> S.Vector n Double
factorAnalysisUniqueness fa =
    let lds = toTensor . snd . split $ toSource fa
        sgs :: Source # Diagonal (MVNMean n) (MVNMean n)
        sgs = fromTensor . toTensor . snd . split. toSource $ linearModelObservableDistribution fa
        cms :: Source # Diagonal (MVNMean n) (MVNMean n)
        cms = fromTensor . toTensor . unsafeMatrixMultiply lds $ transpose lds
     in coordinates $ (sgs - cms) / sgs



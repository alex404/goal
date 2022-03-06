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
      FactorAnalysis
    , factorAnalysisObservableDistribution
    , factorAnalysisExpectationMaximization
    , factorAnalysisUniqueness
    -- * Principle Component Analysis
    , PrincipleComponentAnalysis
    , naturalPCAToLGH
    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models
import Goal.Graphical.Models.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Factor Analysis ---


type FactorAnalysis n k = Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)


type instance Observation (FactorAnalysis n k) = S.Vector n Double

factorAnalysisObservableDistribution
    :: (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> Natural # MultivariateNormal n
factorAnalysisObservableDistribution =
     snd . splitConjugatedHarmonium . transposeHarmonium
     . naturalFactorAnalysisToLGH

factorAnalysisExpectationMaximization
    :: ( KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Natural # FactorAnalysis n k
    -> Natural # FactorAnalysis n k
factorAnalysisExpectationMaximization zs fa =
    transition . sourceFactorAnalysisMaximizationStep . expectationStep zs
        $ naturalFactorAnalysisToLGH fa

factorAnalysisUniqueness
    :: (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> S.Vector n Double
factorAnalysisUniqueness fa =
    let lds = toMatrix . snd . split $ toSource fa
        sgs = S.takeDiagonal . snd . splitMultivariateNormal . toSource
                $ factorAnalysisObservableDistribution fa
        cms = S.takeDiagonal . S.matrixMatrixMultiply lds $ S.transpose lds
     in (sgs - cms) / sgs

-- Internal --

naturalFactorAnalysisToLGH
    :: (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> Natural # LinearGaussianHarmonium n k
naturalFactorAnalysisToLGH fa =
    let (nzs,tns) = split fa
        (mus,vrs) = S.toPair . S.toColumns . S.fromRows . S.map coordinates $ splitReplicated nzs
        mvn = joinNaturalMultivariateNormal mus $ S.diagonalMatrix vrs
        fa' = join mvn tns
     in joinConjugatedHarmonium fa' $ toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1

sourceFactorAnalysisMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => Mean # LinearGaussianHarmonium n k
    -> Source # FactorAnalysis n k
sourceFactorAnalysisMaximizationStep hrm =
    let (nz,nzx,nx) = splitHarmonium hrm
        (muz,etaz) = splitMeanMultivariateNormal nz
        (mux,etax) = splitMeanMultivariateNormal nx
        outrs = toMatrix nzx - S.outerProduct muz mux
        wmtx = S.matrixMatrixMultiply outrs $ S.inverse etax
        zcvr = etaz - S.outerProduct muz muz
        vrs = S.takeDiagonal $ zcvr - S.matrixMatrixMultiply wmtx (S.transpose outrs)
        snrms = joinReplicated $ S.zipWith (curry fromTuple) muz vrs
     in join snrms $ fromMatrix wmtx


--- Principle Component Analysis ---


type PrincipleComponentAnalysis n k = Affine Tensor (MVNMean n) (IsotropicNormal n) (MVNMean k)

naturalPCAToLGH
    :: (KnownNat n, KnownNat k)
    => Natural # PrincipleComponentAnalysis n k
    -> Natural # LinearGaussianHarmonium n k
naturalPCAToLGH pca =
    let (iso,tns) = split pca
        (mus0,vr) = split iso
        mus = coordinates mus0
        sgma = S.diagonalMatrix . S.replicate . S.head $ coordinates vr
        mvn = joinNaturalMultivariateNormal mus sgma
        pca' = join mvn tns
     in joinConjugatedHarmonium pca' $ toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1

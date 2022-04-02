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
module Goal.Graphical.Models.Harmonium.Gaussian where
--    (
--    -- * Factor Analysis
--      factorAnalysisObservableDistribution
--    , factorAnalysisExpectationMaximization
--    , factorAnalysisUniqueness
--    -- * Principle Component Analysis
--    , naturalPCAToLGH
--    , pcaObservableDistribution
--    , pcaExpectationMaximization
--    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models
import Goal.Graphical.Models.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Factor Analysis ---


type instance Observation (FactorAnalysis n k) = S.Vector n Double

--factorAnalysisObservableDistribution
--    :: (KnownNat n, KnownNat k)
--    => Natural # FactorAnalysis n k
--    -> Natural # SymmetricNormal n
--factorAnalysisObservableDistribution =
--     snd . splitConjugatedHarmonium . transposeHarmonium
--     . naturalFactorAnalysisToLGH

--factorAnalysisExpectationMaximization
--    :: ( KnownNat n, KnownNat k)
--    => [S.Vector n Double]
--    -> Natural # FactorAnalysis n k
--    -> Natural # FactorAnalysis n k
--factorAnalysisExpectationMaximization zs fa =
--    transition . sourceFactorAnalysisMaximizationStep . expectationStep zs
--        $ naturalFactorAnalysisToLGH fa

--factorAnalysisUniqueness
--    :: (KnownNat n, KnownNat k)
--    => Natural # FactorAnalysis n k
--    -> S.Vector n Double
--factorAnalysisUniqueness fa =
--    let lds = toMatrix . snd . split $ toSource fa
--        sgs = S.takeDiagonal . snd . splitMultivariateNormal . toSource
--                $ factorAnalysisObservableDistribution fa
--        cms = S.takeDiagonal . S.matrixMatrixMultiply lds $ S.transpose lds
--     in (sgs - cms) / sgs

-- Internal --

--naturalFactorAnalysisToLGH
--    :: forall n k . (KnownNat n, KnownNat k)
--    => Natural # FactorAnalysis n k
--    -> Natural # FullGaussianHarmonium n k
--naturalFactorAnalysisToLGH fa =
--    let (nzs,tns) = split fa
--        (nmu,ndiag) = split nzs
--        mvn = join nmu . fromTensor $ toTensor ndiag
--        fa' = join mvn tns
--        idnt :: Source # Diagonal (MVNMean k) (MVNMean k)
--        idnt = 1
--     in joinConjugatedHarmonium fa' . toNatural . join 0 . fromTensor $ toTensor idnt

--naturalFactorAnalysisToDGH
--    :: (KnownNat n, KnownNat k)
--    => Natural # FactorAnalysis n k
--    -> Natural # DiagonalGaussianHarmonium n k
--naturalFactorAnalysisToDGH fa =
--    let idnt :: Source # Diagonal (MVNMean k) (MVNMean k)
--        idnt = 1
--     in joinConjugatedHarmonium fa $ toNatural . join 0 . fromTensor $ toTensor idnt

--sourceFactorAnalysisMaximizationStep
--    :: forall n k . (KnownNat n, KnownNat k)
--    => Mean # DiagonalGaussianHarmonium n k
--    -> Source # FactorAnalysis n k
--sourceFactorAnalysisMaximizationStep hrm =
--    let (mz,mzx,mx) = splitHarmonium hrm
--        (muz,etaz) = split mz
--        (mux,etax) = split mx
--        outrs :: Source # Tensor (MVNMean n) (MVNMean k)
--        outrs = breakPoint mzx - muz >.< mux
--        invetax :: Source # Tensor (MVNMean k) (MVNMean k)
--        invetax = breakPoint $ inverse etax
--        wmtx = outrs <#> invetax
--        zcvr = etaz - muz >.< muz
--        vrs = zcvr - wmtx <#> transpose outrs
--        snrms = join muz $ fromTensor vrs
--     in join snrms wmtx


--- Principle Component Analysis ---


--naturalPCAToIGH
--    :: (KnownNat n, KnownNat k)
--    => Natural # PrincipleComponentAnalysis n k
--    -> Natural # IsotropicGaussianHarmonium n k
--naturalPCAToIGH pca =
--     joinConjugatedHarmonium pca $ toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1
--
--naturalPCAToLGH
--    :: (KnownNat n, KnownNat k)
--    => Natural # PrincipleComponentAnalysis n k
--    -> Natural # IsotropicGaussianHarmonium n k
--naturalPCAToLGH pca =
--    let (iso,tns) = split pca
--        (mus0,vr) = split iso
--        mus = coordinates mus0
--        sgma = S.diagonalMatrix . S.replicate . S.head $ coordinates vr
--        mvn = joinNaturalMultivariateNormal mus sgma
--        pca' = join mvn tns
--     in joinConjugatedHarmonium pca' $ toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1
--
--pcaExpectationMaximization
--    :: ( KnownNat n, KnownNat k)
--    => [S.Vector n Double]
--    -> Natural # PrincipleComponentAnalysis n k
--    -> Natural # PrincipleComponentAnalysis n k
--pcaExpectationMaximization zs pca =
--    transition . sourcePCAMaximizationStep . expectationStep zs
--        $ naturalPCAToIGH pca
--
--sourcePCAMaximizationStep
--    :: forall n k . (KnownNat n, KnownNat k)
--    => Mean # IsotropicGaussianHarmonium n k
--    -> Source # PrincipleComponentAnalysis n k
--sourcePCAMaximizationStep hrm =
--    let (mz,mzx,mx) = splitHarmonium hrm
--        (muz0,etaz) = split mz
--        (mux,etax) = splitMeanMultivariateNormal mx
--        muz = coordinates muz0
--        outrs = toMatrix mzx - S.outerProduct muz mux
--        wmtx = S.matrixMatrixMultiply outrs $ S.inverse etax
--        wmtxtr = S.transpose wmtx
--        n = fromIntegral $ natVal (Proxy @n)
--        ztr = S.head (coordinates etaz) - S.dotProduct muz muz
--        vr = ztr - 2*S.trace (S.matrixMatrixMultiply wmtx $ S.transpose outrs)
--            + S.trace (S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtx etax) wmtxtr)
--        iso = join (Point muz) $ singleton vr / n
--     in join iso $ fromMatrix wmtx

--pcaExpectationMaximization'
--    :: ( KnownNat n, KnownNat k)
--    => [S.Vector n Double]
--    -> Natural # PrincipleComponentAnalysis n k
--    -> Natural # PrincipleComponentAnalysis n k
--pcaExpectationMaximization' zs pca =
--    transition . sourcePCAMaximizationStep' . expectationStep zs
--        $ naturalPCAToLGH pca

--pcaObservableDistribution
--    :: (KnownNat n, KnownNat k)
--    => Natural # PrincipleComponentAnalysis n k
--    -> Natural # SymmetricNormal n
--pcaObservableDistribution =
--     snd . splitConjugatedHarmonium . transposeHarmonium
--     . naturalPCAToLGH




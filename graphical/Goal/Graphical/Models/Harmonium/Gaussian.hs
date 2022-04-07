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
    , pcaExpectationMaximization
    --, factorAnalysisUniqueness
    ) where
    -- * Principle Component Analysis
--    , naturalPCAToLGH
--    , pcaObservableDistribution
--    , pcaExpectationMaximization

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
        (lkl,nz) = split $ lgh
        (nx,aff) = split lkl
        (nmux,nvrx) = split nx
        nx' = join nmux . fromTensor $ toTensor nvrx
        lkl' = join nx' aff
     in snd . splitConjugatedHarmonium . transposeHarmonium $ join lkl' nz

--factorAnalysisExpectationMaximization
--    :: ( KnownNat n, KnownNat k)
--    => Observations (FactorAnalysis n k)
--    -> Natural # FactorAnalysis n k
--    -> Natural # FactorAnalysis n k
--factorAnalysisExpectationMaximization xs fa =
--    linearModelMaximizationStep . expectationStep xs $ factorAnalysisToDiagonalLGH fa

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
--    let (mx,mxz,mz) = splitHarmonium hrm
--        (mux,etax) = split mx
--        (muz,etaz) = split mz
--        xcvr = etax - mux >.< mux
--        xzcvr = mxz - mux >.< muz
--        invetaz = inverse etaz
--        wmtx = unsafeMatrixMultiply xzcvr invetaz
--        vrs = xcvr - fromTensor (unsafeMatrixMultiply wmtx $ transpose xzcvr)
--        snrms = join mux vrs
--        mfa :: Mean # FactorAnalysis n k
--        mfa = join snrms wmtx
--     in breakChart mfa

linearModelMaximizationStep
    :: forall n f k
     . ( KnownNat n, KnownNat k, Square Mean f (MVNMean n)
       , LinearlyComposable f Tensor (MVNMean n) (MVNMean n) (MVNMean k) )
    => Mean # LinearGaussianHarmonium f n k
    -> Natural # LinearModel f n k
linearModelMaximizationStep hrm =
    let (mx,mxz,mz) = splitHarmonium hrm
        (mmux,mvrx) = split mx
        (mmuz,mvrz) = split mz
        xcvr = mvrx - mmux >.< mmux
        xzcvr = mxz - mmux >.< mmuz
        invmvrz = inverse mvrz
        nvrx0 = inverse $ xcvr - fromTensor (changeOfBasis (transpose xzcvr) invmvrz)
        nvrx = -0.5 .> nvrx0
        nxz = dualComposition nvrx0 xzcvr invmvrz
        nmvn = join (nvrx0 >.> mmux) nvrx
     in join nmvn nxz


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
pcaExpectationMaximization
    :: ( KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Natural # PrincipleComponentAnalysis n k
    -> Natural # PrincipleComponentAnalysis n k
pcaExpectationMaximization xs pca =
    transition . sourcePCAMaximizationStep . expectationStep xs
        $ joinConjugatedHarmonium pca standardNormal

sourcePCAMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => Mean # IsotropicGaussianHarmonium n k
    -> Source # PrincipleComponentAnalysis n k
sourcePCAMaximizationStep hrm =
    let (mx,mxz,mz) = splitHarmonium hrm
        (mux,etax) = split mx
        (muz,etaz) = split mz
        cmux = coordinates mux
        xcvr = S.head (coordinates etax) - S.dotProduct cmux cmux
        xzcvr = mxz - mux >.< muz
        invetaz = inverse etaz
        wmtx = unsafeMatrixMultiply xzcvr invetaz
        n = fromIntegral $ natVal (Proxy @n)
        vr = xcvr - bilinearTrace (changeOfBasis (transpose xzcvr) invetaz)
        iso = join mux $ singleton vr / n
        pca :: Mean # PrincipleComponentAnalysis n k
        pca = join iso wmtx
     in breakChart pca

sourceFactorAnalysisMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => Mean # DiagonalGaussianHarmonium n k
    -> Source # FactorAnalysis n k
sourceFactorAnalysisMaximizationStep hrm =
    let (mx,mxz,mz) = splitHarmonium hrm
        (mux,etax) = split mx
        (muz,etaz) = split mz
        xcvr = etax - mux >.< mux
        xzcvr = mxz - mux >.< muz
        invetaz = inverse etaz
        wmtx = unsafeMatrixMultiply xzcvr invetaz
        vrs = xcvr - fromTensor (unsafeMatrixMultiply wmtx $ transpose xzcvr)
        snrms = join mux vrs
        mfa :: Mean # FactorAnalysis n k
        mfa = join snrms wmtx
     in breakChart mfa


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




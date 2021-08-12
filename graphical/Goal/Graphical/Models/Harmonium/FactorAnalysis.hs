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
module Goal.Graphical.Models.Harmonium.FactorAnalysis
    (
    -- * Factor Analysis
      FactorAnalysis (FactorAnalysis)
    , factorAnalysisObservableDistribution
    , factorAnalysisExpectationMaximization
    ) where

--- Imports ---


import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models
import Goal.Graphical.Models.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Types ---


newtype FactorAnalysis n k =
    FactorAnalysis (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k))
    deriving (Manifold,Product)

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
    transition .sourceFactorAnalysisMaximizationStep . expectationStep zs
        $ naturalFactorAnalysisToLGH fa

-- Factor Analysis --

naturalFactorAnalysisToLGH
    :: (KnownNat n, KnownNat k)
    => Natural # FactorAnalysis n k
    -> Natural # LinearGaussianHarmonium n k
naturalFactorAnalysisToLGH fa =
    let ltnt = toNatural . joinMultivariateNormal 0 $ S.diagonalMatrix 1
        (nzs,tns) = split fa
        (mus,vrs) = S.toPair . S.toColumns . S.fromRows
            . S.map coordinates $ splitReplicated nzs
        cvr = S.diagonalMatrix vrs
        mvn = joinNaturalMultivariateNormal mus cvr
        fa' = join mvn tns
     in joinConjugatedHarmonium fa' ltnt

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


--- Instances ---


instance ( KnownNat n, KnownNat k) => Transition Natural Source (FactorAnalysis n k)where
    transition nfa =
        let (nnrms,nmtx) = split nfa
            (nmu,nsg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates
                $ splitReplicated nnrms
            invsg = -2 * nsg
            ssg = recip invsg
            smu = nmu / invsg
            snrms = joinReplicated $ S.zipWith (curry fromTuple) smu ssg
            smtx = S.matrixMatrixMultiply (S.diagonalMatrix ssg) $ toMatrix nmtx
         in join snrms $ fromMatrix smtx

instance ( KnownNat n, KnownNat k) => Transition Source Natural (FactorAnalysis n k) where
    transition sfa =
        let (snrms,smtx) = split sfa
            (smu,ssg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates
                $ splitReplicated snrms
            invsg = recip ssg
            nmu = invsg * smu
            nsg = -0.5 * invsg
            nmtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) $ toMatrix smtx
            nnrms = joinReplicated $ S.zipWith (curry fromTuple) nmu nsg
         in join nnrms $ fromMatrix nmtx

type instance PotentialCoordinates (FactorAnalysis n k) = Natural

instance ( KnownNat k, KnownNat n )
  => ObservablyContinuous Natural (FactorAnalysis n k) where
      logObservableDensities fa = logDensities (factorAnalysisObservableDistribution fa)

instance ( KnownNat k, KnownNat n )
  => Statistical (FactorAnalysis n k) where
      type SamplePoint (FactorAnalysis n k) = (S.Vector n Double, S.Vector k Double)

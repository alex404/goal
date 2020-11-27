{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.LatentVariable
    ( -- * Latent Variable Models
      LatentVariable
    , ExpectationMaximization (expectationStep)
      -- * Factor Analysis
    , FactorAnalysis
    , joinFactorAnalysis
    , splitFactorAnalysis
    , toMultivariateNormal
    , factorAnalysisExpectationMaximization
    ) where

--- Imports ---


-- Package --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G


--- Latent Variable Class ---


type LatentVariable x o l =
    ( SamplePoint x ~ HList (o : l) )

class Statistical x => ExpectationMaximization c x where
    expectationStep
        :: LatentVariable x o l
        => [o] -> c # x -> Mean # x




--- Types ---


data FactorAnalysis (n :: Nat) (k :: Nat)


--- Construction ---


joinFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => S.Vector n Double -- ^ Mean bias
    -> S.Vector n Double -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> Source # FactorAnalysis n k
joinFactorAnalysis mus vrs mtx =
    Point $ mus S.++ vrs S.++ G.toVector mtx

splitFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => Source # FactorAnalysis n k
    -> (S.Vector n Double, S.Vector n Double, S.Matrix n k Double)
splitFactorAnalysis (Point cs) =
    let (mus,cs') = S.splitAt cs
        (vrs,mtx) = S.splitAt cs'
     in (mus,vrs,G.Matrix mtx)

toMultivariateNormal
    :: (KnownNat n, KnownNat k)
    => Source # FactorAnalysis n k
    -> Source # MultivariateNormal n
toMultivariateNormal fan =
    let (mus,vrs,mtx) = splitFactorAnalysis fan
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = S.diagonalMatrix vrs
     in joinMultivariateNormal mus $ mtx1 + mtx2

factorAnalysisExpectationMaximization
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Source # FactorAnalysis n k
    -> Source # FactorAnalysis n k
factorAnalysisExpectationMaximization xs fan =
    let (_,vrs,wmtx) = splitFactorAnalysis fan
        wmtxtr = S.transpose wmtx
        vrinv = S.pseudoInverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtxtr vrinv) wmtx
        gmtx = S.pseudoInverse $ S.matrixIdentity + mlts
        xht = average xs
        nxht = S.scale (-1) xht
        rsds = [ S.add x nxht | x <- xs ]
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtxtr) vrinv
        muhts = S.matrixVectorMultiply mlts' <$> rsds
        invsgm = S.pseudoInverse . (gmtx +) . S.averageOuterProduct $ zip muhts muhts
        wmtx0 = S.averageOuterProduct $ zip rsds muhts
        wmtx' = S.matrixMatrixMultiply wmtx0 invsgm
        vrs0 = S.withMatrix (S.scale (-1)) . S.matrixMatrixMultiply wmtx $ S.transpose wmtx0
        smtx = S.averageOuterProduct $ zip rsds rsds
        vrs' = S.takeDiagonal $ smtx + vrs0
     in joinFactorAnalysis xht vrs' wmtx'


--- Instances ---


instance (KnownNat n, KnownNat k) => Manifold (FactorAnalysis n k) where
    type Dimension (FactorAnalysis n k) = 2*n + k*n

instance (KnownNat n, KnownNat k) => Statistical (FactorAnalysis n k) where
    type SamplePoint (FactorAnalysis n k) = (S.Vector n Double, S.Vector k Double)



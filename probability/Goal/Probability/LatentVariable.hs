{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

module Goal.Probability.LatentVariable
    ( -- * Factor Analysis
        FactorAnalysis
    , joinFactorAnalysis
    , splitFactorAnalysis
    , toMultivariateNormal
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


--- Types ---


data FactorAnalysis (n :: Nat) (k :: Nat)


--- Construction ---


joinFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => S.Vector n Double -- ^ Mean bias
    -> S.Vector n Double -- ^ Variances
    -> S.Matrix n k Double -- ^ Interaction Parameters
    -> Source # FactorAnalysis n k
{-# INLINE joinFactorAnalysis #-}
joinFactorAnalysis mus vrs mtx =
    Point $ mus S.++ vrs S.++ G.toVector mtx

splitFactorAnalysis
    :: (KnownNat n, KnownNat k)
    => Source # FactorAnalysis n k
    -> (S.Vector n Double, S.Vector n Double, S.Matrix n k Double)
{-# INLINE splitFactorAnalysis #-}
splitFactorAnalysis (Point cs) =
    let (mus,cs') = S.splitAt cs
        (vrs,mtx) = S.splitAt cs'
     in (mus,vrs,G.Matrix mtx)

toMultivariateNormal
    :: (KnownNat n, KnownNat k)
    => Source # FactorAnalysis n k
    -> Source # MultivariateNormal n
{-# INLINE toMultivariateNormal #-}
toMultivariateNormal fan =
    let (mus,vrs,mtx) = splitFactorAnalysis fan
        mtx1 = S.matrixMatrixMultiply mtx (S.transpose mtx)
        mtx2 = S.diagonalMatrix vrs
     in joinMultivariateNormal mus $ addMatrices mtx1 mtx2

factorAnalysisExpectationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> Source # FactorAnalysis n k
    -> [(S.Vector k Double,S.Matrix k k Double)]
{-# INLINE factorAnalysisExpectationStep #-}
factorAnalysisExpectationStep xs fan =
    let (_,vrs,wmtx) = splitFactorAnalysis fan
        wmtx' = S.transpose wmtx
        vrinv = S.inverse $ S.diagonalMatrix vrs
        mlts = S.matrixMatrixMultiply (S.matrixMatrixMultiply wmtx' vrinv) wmtx
        gmtx = S.inverse $ addMatrices S.matrixIdentity mlts
        nxht = S.scale (-1) $ average xs
        mlts' = S.matrixMatrixMultiply (S.matrixMatrixMultiply gmtx wmtx') vrinv
     in do
         x <- xs
         let muht = S.matrixVectorMultiply mlts' (S.add x nxht)
             sgmaht = addMatrices gmtx $ S.outerProduct muht muht
         return (muht,sgmaht)

factorAnalysisMaximizationStep
    :: forall n k . (KnownNat n, KnownNat k)
    => [S.Vector n Double]
    -> [(S.Vector k Double,S.Matrix k k Double)]
    -> Source # FactorAnalysis n k
{-# INLINE factorAnalysisMaximizationStep #-}
factorAnalysisMaximizationStep xs musgms =
    let xht = average xs
        nxht = S.scale (-1) xht
        (mus,sgms) = unzip musgms
        wmtx0 = sumMatrices [S.add x nxht `S.outerProduct` mu | (x,mu) <- zip xs mus ]
        wmtx1 = S.inverse $ sumMatrices sgms
        wmtx = S.matrixMatrixMultiply wmtx0 wmtx1
        smtx0 :: S.Matrix n n Double
        smtx0 = sumMatrices [ S.outerProduct x x | x <- xs ]
        smtx = S.matrixMatrixMultiply (S.transpose smtx0) smtx0
        wmtx2 = S.matrixMatrixMultiply wmtx . S.withMatrix (S.scale (-1/fromIntegral (length xs))) $ S.transpose wmtx0
        vrs = S.takeDiagonal $ addMatrices smtx wmtx2
     in joinFactorAnalysis xht vrs wmtx

addMatrices :: (KnownNat n, KnownNat k) => S.Matrix n k Double -> S.Matrix n k Double -> S.Matrix n k Double
addMatrices (G.Matrix mtx1) (G.Matrix mtx2) = G.Matrix $ S.add mtx1 mtx2

sumMatrices :: (KnownNat n, KnownNat k) => [S.Matrix n k Double] -> S.Matrix n k Double
sumMatrices = foldr1 addMatrices


--- Instances ---


instance (KnownNat n, KnownNat k) => Manifold (FactorAnalysis n k) where
    type Dimension (FactorAnalysis n k) = 2*n + k*n

instance (KnownNat n, KnownNat k) => Statistical (FactorAnalysis n k) where
    type SamplePoint (FactorAnalysis n k) = (S.Vector n Double, S.Vector k Double)



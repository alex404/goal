{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | An Exponential Family 'Harmonium' is a product exponential family with a
particular bilinear structure (<https://papers.nips.cc/paper/2672-exponential-family-harmoniums-with-an-application-to-information-retrieval Welling, et al., 2005>).
A 'Mixture' model is a special case of harmonium.
-}
module Goal.Graphical.Models.Harmonium.Boltzmann where

--- Imports ---

--- Goal

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Models

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Misc

import Data.Proxy (Proxy (Proxy))

--- Types ---

-- | Linear models are linear functions with additive Guassian noise.
type BoltzmannLinearModel f n k = Affine L.Full (StandardNormal n) (MultivariateNormal f n) (Boltzmann n)

linearBoltzmannHarmoniumConjugationParameters ::
    forall n k f.
    ( KnownNat n
    , KnownNat k
    , KnownCovariance f n
    ) =>
    Natural # LinearModel f n k ->
    -- | Conjugation parameters
    (Double, Natural # FullNormal k)
linearBoltzmannHarmoniumConjugationParameters aff =
    let (thts, tht3) = split aff
        (tht1, tht2) = splitNaturalNormal thts
        (itht20, lndt, _) = inverseLogDeterminant . negate $ 2 .> tht2
        itht2 = -2 .> itht20
        tht21 = itht2 >.> tht1
        rho0 = -0.25 * (tht1 <.> tht21) - 0.5 * lndt
        rho1 = -0.5 .> (transpose tht3 >.> tht21)
        rho2 = -0.25 .> changeOfBasis tht3 itht2
     in (rho0, joinNaturalNormal rho1 $ fromTensor rho2)

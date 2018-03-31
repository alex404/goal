-- | Provides a model of a pendulum.
module Goal.Simulation.Physics.Models.Pendulum where
--    ( -- * The Pendulum
--    Pendulum (Pendulum) )where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B

import Goal.Simulation.Physics


--- Penduli ---


-- | The pendulum. Defined by a mass and a length.
data Pendulum l m

pendulumLength0 :: (KnownNat ln, KnownNat ld) => Proxy (ln / ld) -> Position (Pendulum (ln / ld) m) x -> Rational
pendulumLength0 prxyl _ = ratVal prxyl

pendulumLength :: (KnownNat ln, KnownNat ld) => Position (Pendulum (ln / ld) m) x -> Rational
pendulumLength = pendulumLength0 Proxy

pendulumMass0 :: (KnownNat mn, KnownNat md) => Proxy (mn / md) -> Position (Pendulum l (mn / md)) x -> Rational
pendulumMass0 prxym _ = ratVal prxym

pendulumMass :: (KnownNat mn, KnownNat md) => Position (Pendulum l (mn / md)) x -> Rational
pendulumMass = pendulumMass0 Proxy

--- Instances ---

instance (KnownNat ln, KnownNat ld, KnownNat mn, KnownNat md) => Manifold (Pendulum (ln/ld) (mn/md)) where
    type Dimension (Pendulum (ln/ld) (mn/md)) = 1

instance (KnownNat ln, KnownNat ld, KnownNat mn, KnownNat md) => Riemannian Generalized (Pendulum (ln/ld) (mn/md)) where
    metric q = Point . B.singleton . realToFrac $ pendulumMass q * pendulumLength q^(2 :: Int)

instance (KnownNat ln, KnownNat ld, KnownNat mn, KnownNat md) => Conservative Gravity (Pendulum (ln/ld) (mn/md)) where
    potentialEnergy (Gravity g) q =
        let m = realToFrac $ pendulumMass q
            l = realToFrac $ pendulumLength q
            tht = B.head $ coordinates q
         in m * realToFrac g * l * (1 - cos tht)

instance (KnownNat ln, KnownNat ld, KnownNat mn, KnownNat md) => ForceField Gravity (Pendulum (ln/ld) (mn/md)) where
    force = conservativeForce

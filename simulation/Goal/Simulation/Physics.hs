{-# LANGUAGE MultiParamTypeClasses,FlexibleInstances,FlexibleContexts #-}
-- | This module provides a general interface for working with mechanical systems.
module Goal.Simulation.Physics where
--    ( -- * Configurations
--    -- ** Types
--      Generalized
--    , PhaseSpace
--    , Velocity
--    , Phase
--    , Momentum
--    , DualPhase
--    , Acceleration
--    , Force
--    , Intertia
--    -- ** Accessors
--    , position
--    , velocity
--    , momentum
--    -- * ForceFields
--    -- ** Classes
--    , ForceField (force)
--    , Conservative (potentialEnergy)
--    , vectorField
--    , mechanicalEnergy
--    -- ** Types
--    , Gravity (Gravity)
--    , earthGravity
--    , Damping (Damping)
--    -- * Util
--    , periodic
--    , revolutions
--    ) where

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry


--- Configuration Space ---


-- Charts --

-- | Some generalized coordinates of a mechanical system.
data Generalized

-- Manifolds --

type Position m x = Point Generalized m x

-- | The 'Tangent' space of a mechanical system in 'Generalized' coordinates is the space of 'GeneralizedVelocity's.
type Velocity m x = TangentVector Generalized m x

type Phase m x = TangentPair Generalized m x

-- | The 'Tangent' space of a mechanical system in 'Generalized' coordinates is the space of 'GeneralizedVelocity's.
type Momentum m x = CotangentVector Generalized m x

type DualPhase m x = CotangentPair Generalized m x

-- | The second order 'Tangent' space of a mechanical system is the space of 'GeneralizedAcceleration's.
type Acceleration m x = TangentVector Directional (TangentSpace Generalized m) x

-- | The second order 'Tangent' space of a mechanical system is the space of 'GeneralizedAcceleration's.
type Force m x = CotangentVector Differential (TangentSpace Generalized m) x

-- | The tangent 'Bundle' of a dynamical system is also known as the 'PhaseSpace'. The "state" of a
-- mechanical system is typically understood to be an element of the 'PhaseSpace'.
type PhaseSpace m = TangentBundle Generalized m

-- | The 'Riemannian' metric on a mechanical system is known as the 'GeneralizedInertia'.
type Intertia m x = CotangentTensor Directional (TangentSpace Generalized m) x

{-
instance Primal Hamiltonian where
    type Dual Hamiltonian = Lagrangian

instance Primal Lagrangian where
    type Dual Lagrangian = Hamiltonian
-}


-- Functions --

-- | Returns the 'position' coordinates of a mechanical system.
position :: Manifold m => Phase m x -> Position m x
position = projectTangentPair

kineticEnergy :: (Riemannian Generalized m, RealFloat x) => Phase m x -> x
kineticEnergy dq = 0.5 * (flat dq <.> dq)

-- Force fields --

-- | Gravitational force.
newtype Gravity = Gravity Double

earthGravity :: Gravity
-- | Gravity on the surface of the earth.
earthGravity = Gravity 9.80665

-- | A (non-conservative) damping force.
newtype Damping = Damping Double

-- | A 'ForceField' describes how to attach a force vector (an element of the dual space of
-- generalized accelerations) to every point in the phase space.  Note that a 'ForceField' is not
-- necessarily 'Conservative'.
class Manifold m => ForceField f m where
    force :: RealFloat x => f -> Phase m x -> Force m x

-- | A 'Conservative' force depends only on 'position's and can be described as the gradient of a
-- 'potentialEnergy'.
class Manifold m => Conservative f m where
    potentialEnergy :: RealFloat x => f -> Position m x -> x

conservativeForce :: (Conservative f m, RealFloat x) => f -> Phase m x -> Force m x
conservativeForce f qdq =
    let q = position qdq
     in Point . coordinates $ differential (potentialEnergy f) q

-- | The 'vectorField' function takes a 'ForceField' on a mechanical system and converts it into the
-- appropriate element of the 'Tangent' space of the 'PhaseSpace'.
vectorField :: (Riemannian Generalized m, ForceField f m, RealFloat x)
    => f
    -> Phase m x
    -> TangentVector Directional (PhaseSpace m) x
vectorField f qdq =
    let smtx = fromJust . matrixInverse . toMatrix . metric $ position qdq
     in Point $ joinV (coordinates $ detachTangentVector qdq) (matrixVectorMultiply smtx . coordinates $ force f qdq)

-- | Computes the kinetic and potential energy of the given state of a mechanical system.

{-
mechanicalEnergy
    :: Riemannian Generalized m
    => Phase m x
    -> (Double, Double)
mechanicalEnergy f qdq =
    let v = kineticEnergy $ velocity qdq
     in (v, v - potential
-}

--- Functions ---


{-
-- | Takes a coordinate vector and clips the specified rotational degrees of
-- freedom to lie within the bounds of negative pi and pi.
periodic :: Manifold m
    => [Bool] -- ^ Which degrees of freedom are rotational
    -> (Partials :#: PhaseSpace m) -- ^ The current state
    -> (Partials :#: PhaseSpace m) -- ^ The clipped state
periodic bls qdq =
    let q = position qdq
        q' = fromList (manifold q) $ clipper <$> zip bls (listCoordinates q)
        dq' = fromCoordinates (Tangent q') . coordinates $ velocity qdq
     in tangentToBundle dq'
       where clipper (True,x) = x - 2 * pi * fromIntegral (revolutions x)
             clipper (False,x) = x

-- | Counts the number of revolutions made by an unclipped rotational degree of
-- freedom.
revolutions :: Double -> Int
revolutions x
    | x >= pi = floor $ (x + pi) / (2*pi)
    | x <= -pi = ceiling $ (x - pi) / (2*pi)
    | otherwise = 0

-}

--- Instances ---


instance (ForceField f m, ForceField g m) => ForceField (f,g) m where
    force (f,g) qdq = force f qdq  <+> force g qdq

instance (ForceField f m, ForceField g m, ForceField h m) => ForceField (f,g,h) m where
    force (f,g,h) qdq = force f qdq  <+> force g qdq <+> force h qdq


{-
instance Manifold c m => ForceField (Partials :#: PhaseSpace m -> Differentials :#: GeneralizedAcceleration m) m where
    force f = f

instance Manifold m => ForceField (Differentials :#: GeneralizedAcceleration m) m where
    force = const
    -}

-- Damping --

instance Manifold m => ForceField Damping m where
    force (Damping c) qdq =
        Point $ negate . (*realToFrac c) <$> coordinates (detachTangentVector qdq)


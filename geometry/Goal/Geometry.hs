-- | The main module of @goal-geometry@. Import this module to use all the
-- features provided by this library.
module Goal.Geometry
    (
    -- * Re-Exports
    module Goal.Geometry.Manifold
    , module Goal.Geometry.Linear
    , module Goal.Geometry.Map
    , module Goal.Geometry.Map.Multilinear
    , module Goal.Geometry.Map.Multilinear.Convolutional
    , module Goal.Geometry.Map.NeuralNetwork
    , module Goal.Geometry.Differential
    , module Goal.Geometry.Differential.GradientPursuit
    ) where


-- Imports --


-- Re-exports --

import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Map
import Goal.Geometry.Map.Multilinear
import Goal.Geometry.Map.Multilinear.Convolutional
import Goal.Geometry.Map.NeuralNetwork
import Goal.Geometry.Differential
import Goal.Geometry.Differential.GradientPursuit

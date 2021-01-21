-- | The main module of @goal-geometry@. Import this module to use all the
-- features provided by this library.
module Goal.Geometry
    (
    -- * Re-Exports
    module Goal.Geometry.Manifold
    , module Goal.Geometry.Vector
    , module Goal.Geometry.Map
    , module Goal.Geometry.Map.Linear
    , module Goal.Geometry.Map.Linear.Convolutional
    , module Goal.Geometry.Map.NeuralNetwork
    , module Goal.Geometry.Differential
    , module Goal.Geometry.Differential.GradientPursuit
    ) where


-- Imports --


-- Re-exports --

import Goal.Geometry.Manifold
import Goal.Geometry.Vector
import Goal.Geometry.Map
import Goal.Geometry.Map.Linear
import Goal.Geometry.Map.Linear.Convolutional
import Goal.Geometry.Map.NeuralNetwork
import Goal.Geometry.Differential
import Goal.Geometry.Differential.GradientPursuit

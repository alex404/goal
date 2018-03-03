-- | Training multilayer perceptrons to compute predictions for filtering by
-- contrastive divergence minimization.

module Goal.Cognition.Filter.Optimal where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Cognition.Filter.Approximate

--- Backpropagation ---

discretePrediction
    :: Discrete x
    => (Element x -> Standard :#: Categorical x) -- ^ Transition distribution
    -> (Standard :#: Categorical x) -- ^ Prior Distribution
    -> (Standard :#: Categorical x) -- ^ Prediction Distribution
discretePrediction trns p =
     foldl1' (<+>) [ density p x .> trns x | x <- samples $ manifold p ]

discreteInference
    :: (Discrete x, AbsolutelyContinuous c y)
    => (Element x -> c :#: y) -- ^ Emission distribution
    -> Sample y
    -> (Standard :#: Categorical x) -- ^ Prior Distribution
    -> (Standard :#: Categorical x) -- ^ Prediction Distribution
discreteInference emsn y p =
    let xm = manifold p
        xs = samples xm
        pxs = density p <$> xs
        emsnxs = [ density (emsn x) y | x <- xs ]
        unnrms = zipWith (*) pxs emsnxs
     in fromList xm . take (dimension xm) $ (/sum unnrms) <$> unnrms

kalmanPrediction
    :: H.Matrix Double -- ^ System Dynamics
    -> H.Matrix Double -- ^ Step-Dependent Process Noise
    -> (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Standard :#: MultivariateNormal) -- ^ Prediction Distribution
kalmanPrediction xdnm =
    extendedKalmanPrediction (xdnm H.#>) (const xdnm)

kalmanPrediction1D
    :: Double -- ^ System Dynamics
    -> Double -- ^ Step-Dependent Process Noise
    -> (Standard :#: Normal) -- ^ Prior Distribution
    -> (Standard :#: Normal) -- ^ Prediction Distribution
kalmanPrediction1D xdnm =
    extendedKalmanPrediction1D (xdnm *) (const xdnm)

kalmanInference
    :: H.Matrix Double -- ^ Observation Model
    -> H.Matrix Double -- ^ Observation Covariance
    -> Coordinates -- ^ Observation
    -> (Standard :#: MultivariateNormal) -- ^ Prior Distribution
    -> (Standard :#: MultivariateNormal) -- ^ Prediction Distribution
kalmanInference zdnm = extendedKalmanInference (zdnm H.#>) (const zdnm)

kalmanInference1D
    :: Double -- ^ System Dynamics
    -> Double -- ^ Step-Dependent Process Noise
    -> Double -- ^ Observation
    -> (Standard :#: Normal) -- ^ Prior Distribution
    -> (Standard :#: Normal) -- ^ Inference Distribution
kalmanInference1D zdnm =
    extendedKalmanInference1D (zdnm *) (const zdnm)

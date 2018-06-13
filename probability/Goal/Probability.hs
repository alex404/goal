-- | The main module of goal-probability. Import this module to use all the features provided by this library.
module Goal.Probability
    ( -- * Package Exports
      module Goal.Probability.Statistical
    , module Goal.Probability.ExponentialFamily
    , module Goal.Probability.Distributions
    , module Goal.Probability.ExponentialFamily.NeuralNetwork
    , module Goal.Probability.ExponentialFamily.PopulationCode
    , module Goal.Probability.ExponentialFamily.Harmonium
    , module Goal.Probability.ExponentialFamily.Harmonium.Inference
      -- * Utility
    , resampleVector
    , noisyFunction
    , seed
    -- * External Exports
    , module System.Random.MWC
    , module System.Random.MWC.Probability
    ) where


--- Imports ---


-- Re-exports --

import System.Random.MWC (Seed,save,restore)

import System.Random.MWC.Probability hiding (initialize,sample)
--import System.Random.MWC.Distributions (uniformShuffle)

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily.NeuralNetwork
import Goal.Probability.ExponentialFamily.PopulationCode
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Inference

-- Package --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B


--- Stochastic Functions ---

-- | Creates a seed for later RandST usage.
seed :: Random s Seed
{-# INLINE seed #-}
seed = Prob save

-- | Returns a uniform sample of elements from the given vector.
resampleVector :: (KnownNat n, KnownNat k) => B.Vector n x -> Random s (B.Vector k x)
{-# INLINE resampleVector #-}
resampleVector xs = do
    ks <- B.replicateM $ uniformR (0, B.length xs-1)
    return $ B.backpermute xs ks

-- | Returns a sample from the given function with added noise.
noisyFunction
    :: (Generative c m, Num (SamplePoint m))
    => Point c m -- ^ Noise model
    -> (y -> SamplePoint m) -- ^ Function
    -> y -- ^ Input
    -> Random s (SamplePoint m) -- ^ Stochastic Output
noisyFunction m f x = do
    ns <- samplePoint m
    return $ f x + ns

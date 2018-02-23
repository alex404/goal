-- | The main module of goal-probability. Import this module to use all the features provided by this library.
module Goal.Probability
    ( -- * Package Exports
      module Goal.Probability.Statistical
    , module Goal.Probability.ExponentialFamily
    , module Goal.Probability.Distributions
    , module Goal.Probability.ExponentialFamily.NeuralNetwork
    , module Goal.Probability.ExponentialFamily.NeuralNetwork.Convolutional
    , module Goal.Probability.ExponentialFamily.Harmonium
    , module Goal.Probability.ExponentialFamily.Harmonium.Rectification
--    , module Goal.Probability.ExponentialFamily.Harmonium.Deep
      -- * Utility
    , randomElement
    , noisyFunction
    , seed
    -- * External Exports
    , module System.Random.MWC
    , module System.Random.MWC.Probability
    ) where


--- Imports ---


-- Re-exports --

import System.Random.MWC (Seed,save,restore)

import System.Random.MWC.Probability hiding (initialize)
--import System.Random.MWC.Distributions (uniformShuffle)

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily.NeuralNetwork
import Goal.Probability.ExponentialFamily.NeuralNetwork.Convolutional
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Rectification
--import Goal.Probability.ExponentialFamily.Harmonium.Deep

-- Package --

import Goal.Geometry

--- Stochastic Functions ---

-- | Creates a seed for later RandST usage.
seed :: Random s Seed
seed = Prob save

-- | Returns a random element from a list.
randomElement :: [x] -> Random s x
randomElement xs = do
    u <- uniform
    let elm = round $ fromIntegral (length xs - 1) * (u :: Double)
    return $ xs !! elm

-- | Returns a sample from the given function with added noise.
noisyFunction
    :: (Generative c m, Num (Sample m), RealFloat x)
    => Point c m x -- ^ Noise model
    -> (y -> Sample m) -- ^ Function
    -> y -- ^ Input
    -> Random s (Sample m) -- ^ Stochastic Output
noisyFunction m f x = do
    ns <- generate m
    return $ f x + ns


{-
-- | Returns a random element from a list.
randomElement' :: V.Vector x -> RandST r x
randomElement' xs = do
    let n = V.length xs
    u <- uniformR (0,n-1)
    return $ xs V.! u

-- | Shuffles a list.
shuffleList :: [x] -> RandST r [x]
shuffleList xs = do
    let v = V.fromList xs
    v' <- toRand $ uniformShuffle v
    return $ V.toList v'

-- | Returns a set of samples from the given function with additive Gaussian noise.
noisyRange
    :: Double -- ^ The min of the function input
    -> Double -- ^ The max function input
    -> Int -- ^ Number of samples to draw from the function
    -> Double -- ^ Standard deviation of the noise
    -> (Double -> Double) -- ^ Mixture function
    -> RandST s [(Double,Double)]
noisyRange mn mx n sd f = do
    let xs = range mn mx n
        d = Standard # fromList Normal [0,sd^2]
    fxs <- mapM (\x -> (+ f x) <$> generate d) xs
    return $ zip xs fxs
    -}

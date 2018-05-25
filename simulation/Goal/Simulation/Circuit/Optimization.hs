-- | A collection of 'Circuit's for computing the differentials in gradient
-- descent algorithms.

module Goal.Simulation.Circuit.Optimization
    ( gradientAscent
    , momentumAscent
    , adamAscent
    ) where


--- Imports ---

import Goal.Core
import Goal.Geometry

import Goal.Simulation.Circuit

-- | A 'Circuit' for classic gradient descent.
gradientAscent
    :: Manifold m
    => Double -- ^ Learning Rate
    -> Circuit (TangentPair c m) (Point c m) -- ^ Gradient Ascent
{-# INLINE gradientAscent #-}
gradientAscent eps = arr (gradientStep' eps)

-- | A 'Circuit' for gradient descent with momentum.
momentumAscent
    :: Manifold m
    => Double -- ^ Learning Rate
    -> (Int -> Double) -- ^ Momentum Schedule
    -> Circuit (TangentPair c m) (Point c m) -- ^ Momentum Ascent
{-# INLINE momentumAscent #-}
momentumAscent eps mu = accumulateFunction (0,Nothing) $ \pdp (k,mm) ->
            let m = fromMaybe zero mm
                (p',m') = momentumStep eps (mu k) pdp m
             in (p',(k+1,Just m'))

-- | A 'Circuit' for gradient descent with momentum based on the Adam algorithm.
adamAscent
    :: Manifold m
    => Double -- ^ Learning Rate
    -> Double -- ^ First Moment Rate
    -> Double -- ^ Second Moment Rate
    -> Double -- ^ Second Moment regularizer
    -> Circuit (TangentPair c m) (Point c m) -- ^ Momentum Ascent
{-# INLINE adamAscent #-}
adamAscent eps b1 b2 rg = accumulateFunction (1,Nothing,Nothing) $ \dp (k,mm,mv) ->
            let m = fromMaybe zero mm
                v = fromMaybe zero mv
                (p',m',v') = adamStep eps b1 b2 rg k dp m v
             in (p',(k+1,Just m',Just v'))

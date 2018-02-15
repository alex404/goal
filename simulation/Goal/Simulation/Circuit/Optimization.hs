module Goal.Simulation.Circuit.Optimization where


--- Imports ---

import Goal.Core
import Goal.Geometry

import Goal.Simulation.Circuit


vanillaGradientAscent
    :: (Manifold m, RealFloat x)
    => x -- ^ Learning Rate
    -> Circuit (TangentPair c m x) (Point c m x) -- ^ Gradient Ascent
vanillaGradientAscent eps = arr (gradientStep' eps)

momentumAscent
    :: (Manifold m, RealFloat x)
    => x -- ^ Learning Rate
    -> (Int -> x) -- ^ Momentum Schedule
    -> Circuit (TangentPair c m x) (Point c m x) -- ^ Momentum Ascent
momentumAscent eps mu = accumulateFunction (0,Nothing) $ \pdp (k,mm) ->
            let m = fromMaybe zero mm
                (p',m') = momentumStep eps (mu k) pdp m
             in (p',(k+1,Just m'))

adamAscent
    :: (Manifold m, RealFloat x)
    => x -- ^ Learning Rate
    -> x -- ^ First Moment Rate
    -> x -- ^ Second Moment Rate
    -> x -- ^ Second Moment regularizer
    -> Circuit (TangentPair c m x) (Point c m x) -- ^ Momentum Ascent
adamAscent eps b1 b2 rg = accumulateFunction (1,Nothing,Nothing) $ \dp (k,mm,mv) ->
            let m = fromMaybe zero mm
                v = fromMaybe zero mv
                (p',m',v') = adamStep eps b1 b2 rg k dp m v
             in (p',(k+1,Just m',Just v'))

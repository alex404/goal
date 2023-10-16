-- | Gradient pursuit-based optimization on manifolds.
module Goal.Geometry.Differential.GradientPursuit (
    -- * Cauchy Sequences
    cauchyLimit,
    cauchySequence,

    -- * Gradient Pursuit
    vanillaGradient,
    gradientStep,

    -- ** Algorithms
    GradientPursuit (Classic, Momentum, Adam),
    gradientPursuitStep,
    gradientSequence,
    vanillaGradientSequence,
    gradientCircuit,
    vanillaGradientCircuit,

    -- *** Defaults
    defaultMomentumPursuit,
    defaultAdamPursuit,
) where

--- Imports ---

-- Goal --

import Goal.Core

import Goal.Geometry.Differential
import Goal.Geometry.Manifold
import Goal.Geometry.Vector

import Goal.Core.Vector.Storable qualified as S

--- Cauchy Sequences ---

{- | Attempts to calculate the limit of a sequence by finding the iteration
with a sufficiently small distance from its previous iteration.
-}
cauchyLimit ::
    -- | Distance (divergence) from previous to next
    (c # x -> c # x -> Double) ->
    -- | Epsilon
    Double ->
    -- | Input sequence
    [c # x] ->
    c # x
{-# INLINE cauchyLimit #-}
cauchyLimit f eps ps = last $ cauchySequence f eps ps

-- | Attempts to calculate the limit of a sequence. Returns the list up to the limit.
cauchySequence ::
    -- | Distance (divergence) from previous to next
    (c # x -> c # x -> Double) ->
    -- | Epsilon
    Double ->
    -- | Input list
    [c # x] ->
    -- | Truncated list
    [c # x]
{-# INLINE cauchySequence #-}
cauchySequence f eps ps =
    let pps = takeWhile taker . zip ps $ tail ps
     in head ps : fmap snd pps
  where
    taker (p1, p2) = eps < f p1 p2

--- Gradient Pursuit ---

{- | Ignore the Riemannian metric, and convert a 'Point' from a 'Dual' space to
its 'Primal' space.
-}
vanillaGradient :: (Manifold x) => c #* x -> c # x
{-# INLINE vanillaGradient #-}
vanillaGradient = breakChart

{- | 'gradientStep' takes a step size, a 'Point', a tangent vector at that
point, and returns a 'Point' with coordinates that have moved in the
direction of the tangent vector.
-}
gradientStep ::
    (Manifold x) =>
    Double ->
    -- | Point
    c # x ->
    -- | Tangent Vector
    c # x ->
    -- | Stepped point
    c # x
{-# INLINE gradientStep #-}
gradientStep eps (Point xs) pd =
    Point $ xs + coordinates (eps .> pd)

-- | An ADT reprenting three basic gradient descent algorithms.
data GradientPursuit
    = Classic
    | Momentum (Int -> Double)
    | Adam Double Double Double

-- | A standard momentum schedule.
defaultMomentumPursuit :: Double -> GradientPursuit
{-# INLINE defaultMomentumPursuit #-}
defaultMomentumPursuit mxmu = Momentum fmu
  where
    fmu k = min mxmu $ 1 - 2 ** ((negate 1 -) . logBase 2 . fromIntegral $ div k 250 + 1)

-- | Standard Adam parameters.
defaultAdamPursuit :: GradientPursuit
{-# INLINE defaultAdamPursuit #-}
defaultAdamPursuit = Adam 0.9 0.999 1e-8

-- | A single step of a gradient pursuit algorithm.
gradientPursuitStep ::
    (Manifold x) =>
    -- | Learning Rate
    Double ->
    -- | Gradient pursuit algorithm
    GradientPursuit ->
    -- | Algorithm step
    Int ->
    -- | The point
    c # x ->
    -- | The derivative
    c # x ->
    -- | The velocities
    [c # x] ->
    -- | The updated point and velocities
    (c # x, [c # x])
{-# INLINE gradientPursuitStep #-}
gradientPursuitStep eps Classic _ cp dp _ = (gradientStep eps cp dp, [])
gradientPursuitStep eps (Momentum fmu) k cp dp (v : _) =
    let (p, v') = momentumStep eps (fmu k) cp dp v
     in (p, [v'])
gradientPursuitStep eps (Adam b1 b2 rg) k cp dp (m : v : _) =
    let (p, m', v') = adamStep eps b1 b2 rg k cp dp m v
     in (p, [m', v'])
gradientPursuitStep _ _ _ _ _ _ = error "Momentum list length mismatch in gradientPursuitStep"

-- | Gradient ascent based on the 'Riemannian' metric.
gradientSequence ::
    (Riemannian c x) =>
    -- | Differential calculator
    (c # x -> c #* x) ->
    -- | Step size
    Double ->
    -- | Gradient pursuit algorithm
    GradientPursuit ->
    -- | The initial point
    c # x ->
    -- | The gradient ascent
    [c # x]
{-# INLINE gradientSequence #-}
gradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0, (repeat 0, 0))
  where
    iterator (p, (vs, k)) =
        let dp = sharp p $ f p
            (p', vs') = gradientPursuitStep eps gp k p dp vs
         in (p', (vs', k + 1))

-- | Gradient ascent which ignores the 'Riemannian' metric.
vanillaGradientSequence ::
    (Manifold x) =>
    -- | Differential calculator
    (c # x -> c #* x) ->
    -- | Step size
    Double ->
    -- | Gradient pursuit algorithm
    GradientPursuit ->
    -- | The initial point
    c # x ->
    -- | The gradient ascent
    [c # x]
{-# INLINE vanillaGradientSequence #-}
vanillaGradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0, (repeat 0, 0))
  where
    iterator (p, (vs, k)) =
        let dp = vanillaGradient $ f p
            (p', vs') = gradientPursuitStep eps gp k p dp vs
         in (p', (vs', k + 1))

-- | A 'Circuit' for gradient descent.
gradientCircuit ::
    (Monad m, Manifold x) =>
    -- | Learning Rate
    Double ->
    -- | Gradient pursuit algorithm
    GradientPursuit ->
    -- | (Point, Gradient) to Updated Point
    Circuit m (c # x, c # x) (c # x)
{-# INLINE gradientCircuit #-}
gradientCircuit eps gp = accumulateFunction (repeat 0, 0) $ \(p, dp) (vs, k) -> do
    let (p', vs') = gradientPursuitStep eps gp k p dp vs
    return (p', (vs', k + 1))

-- | A 'Circuit' for gradient descent.
vanillaGradientCircuit ::
    (Monad m, Manifold x) =>
    -- | Learning Rate
    Double ->
    -- | Gradient pursuit algorithm
    GradientPursuit ->
    -- | (Point, Gradient) to Updated Point
    Circuit m (c # x, c #* x) (c # x)
{-# INLINE vanillaGradientCircuit #-}
vanillaGradientCircuit eps gp = second (arr vanillaGradient) >>> gradientCircuit eps gp

--- Internal ---

momentumStep ::
    (Manifold x) =>
    -- | The learning rate
    Double ->
    -- | The momentum decay
    Double ->
    -- | The subsequent TangentPair
    c # x ->
    -- | The subsequent TangentPair
    c # x ->
    -- | The current velocity
    c # x ->
    -- | The (subsequent point, subsequent velocity)
    (c # x, c # x)
{-# INLINE momentumStep #-}
momentumStep eps mu p fd v =
    let v' = eps .> fd + mu .> v
     in (gradientStep 1 p v', v')

adamStep ::
    (Manifold x) =>
    -- | The learning rate
    Double ->
    -- | The first momentum rate
    Double ->
    -- | The second momentum rate
    Double ->
    -- | Second moment regularizer
    Double ->
    -- | Algorithm step
    Int ->
    -- | The subsequent gradient
    c # x ->
    -- | The subsequent gradient
    c # x ->
    -- | First order velocity
    c # x ->
    -- | Second order velocity
    c # x ->
    -- | Subsequent (point, first velocity, second velocity)
    (c # x, c # x, c # x)
{-# INLINE adamStep #-}
adamStep eps b1 b2 rg k0 p fd m v =
    let k = k0 + 1
        fd' = S.map (^ (2 :: Int)) $ coordinates fd
        m' = (1 - b1) .> fd + b1 .> m
        v' = (1 - b2) .> Point fd' + b2 .> v
        mhat = (1 - b1 ^ k) /> m'
        vhat = (1 - b2 ^ k) /> v'
        fd'' = S.zipWith (/) (coordinates mhat) . S.map ((+ rg) . sqrt) $ coordinates vhat
     in (gradientStep eps p $ Point fd'', m', v')

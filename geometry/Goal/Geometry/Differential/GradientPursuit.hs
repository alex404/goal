-- | Gradient pursuit based optimization on manifolds.

module Goal.Geometry.Differential.GradientPursuit
    ( -- * Cauchy Sequences
      cauchyLimit
    , cauchySequence
    -- * Gradient Pursuit
    , vanillaGradient
    , gradientStep
    -- ** Algorithms
    , GradientPursuit (Classic,Momentum,Adam)
    , gradientPursuitStep
    , gradientSequence
    , vanillaGradientSequence
    , gradientCircuit
    -- *** Defaults
    , defaultMomentumPursuit
    , defaultAdamPursuit
    -- *** Regularizers
    , weightDecay
    ) where


--- Imports ---


-- Goal --

import Goal.Core

import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Storable as S


--- Cauchy Sequences ---


-- | Attempts to calculate the limit of a sequence. This finds the iterate with a sufficiently small
-- distance from the previous iterate.
cauchyLimit
    :: (c # x -> c # x -> Double) -- ^ Distance (divergence) from previous to next
    -> Double -- ^ Epsilon
    -> [c # x] -- ^ Input sequence
    -> c # x
{-# INLINE cauchyLimit #-}
cauchyLimit f eps ps = last $ cauchySequence f eps ps

-- | Attempts to calculate the limit of a sequence. Returns the list up to the limit.
cauchySequence
    :: (c # x -> c # x -> Double) -- ^ Distance (divergence) from previous to next
    -> Double -- ^ Epsilon
    -> [c # x] -- ^ Input list
    -> [c # x] -- ^ Truncated list
{-# INLINE cauchySequence #-}
cauchySequence f eps ps =
    let pps = takeWhile taker . zip ps $ tail ps
     in head ps : fmap snd pps
       where taker (p1,p2) = eps < f p1 p2


--- Gradient Pursuit ---

-- | Ignore the Riemannian metric.
vanillaGradient :: Manifold x => c #* x -> c # x
{-# INLINE vanillaGradient #-}
vanillaGradient = breakPoint

-- | 'gradientStep' takes a step size, the location of a 'TangentVector', the
-- 'TangentVector' itself, and returns a 'Point' with coordinates that have
-- moved in the direction of the 'TangentVector'.
gradientStep
    :: Manifold x
    => Double
    -> c # x -- ^ Point
    -> c # x -- ^ Tangent Vector
    -> c # x -- ^ Stepped point
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
    where fmu k = min mxmu $ 1 - 2**((negate 1 -) . logBase 2 . fromIntegral $ div k 250 + 1)

-- | A standard momentum schedule.
defaultAdamPursuit :: GradientPursuit
{-# INLINE defaultAdamPursuit #-}
defaultAdamPursuit = Adam 0.9 0.999 1e-8

-- | A single step of a gradient pursuit algorithm.
gradientPursuitStep
    :: Manifold x
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Int -- ^ Algorithm step
    -> c # x -- ^ The point
    -> c # x -- ^ The derivative
    -> [c # x] -- ^ The velocities
    -> (c # x, [c # x]) -- ^ The updated point and velocities
gradientPursuitStep eps Classic _ cp dp _ = (gradientStep eps cp dp,[])
gradientPursuitStep eps (Momentum fmu) k cp dp (v:_) =
    let (p,v') = momentumStep eps (fmu k) cp dp v
     in (p,[v'])
gradientPursuitStep eps (Adam b1 b2 rg) k cp dp (m:v:_) =
    let (p,m',v') = adamStep eps b1 b2 rg k cp dp m v
     in (p,[m',v'])
gradientPursuitStep _ _ _ _ _ _ = error "Momentum list length mismatch in gradientPursuitStep"

-- | Gradient ascent based on the 'Riemannian' metric.
gradientSequence
    :: Riemannian c x
    => (c # x -> c #* x)  -- ^ Gradient calculator
    -> Double -- ^ Step size
    -> GradientPursuit  -- ^ Gradient pursuit algorithm
    -> c # x -- ^ The initial point
    -> [c # x] -- ^ The gradient ascent
{-# INLINE gradientSequence #-}
gradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0,(repeat zero,0))
        where iterator (p,(vs,k)) =
                  let dp = sharp p $ f p
                      (p',vs') = gradientPursuitStep eps gp k p dp vs
                   in (p',(vs',k+1))

-- | Gradient ascent based on the 'Riemannian' metric.
vanillaGradientSequence
    :: Manifold x
    => (c # x -> c #* x)  -- ^ Gradient calculator
    -> Double -- ^ Step size
    -> GradientPursuit  -- ^ Gradient pursuit algorithm
    -> c # x -- ^ The initial point
    -> [c # x] -- ^ The gradient ascent
{-# INLINE vanillaGradientSequence #-}
vanillaGradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0,(repeat zero,0))
        where iterator (p,(vs,k)) =
                  let dp = vanillaGradient $ f p
                      (p',vs') = gradientPursuitStep eps gp k p dp vs
                   in (p',(vs',k+1))

-- | A 'Circuit' for classic gradient descent.
gradientCircuit
    :: (Monad m, Riemannian c x)
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> c # x -- ^ Initial Point
    -> Circuit m (c #* x) (c # x) -- ^ Gradient Ascent
{-# INLINE gradientCircuit #-}
gradientCircuit eps gp p0 = accumulateFunction (p0,repeat zero,0) $ \pd (p,vs,k) -> do
    let dp = sharp p pd
        (p',vs') = gradientPursuitStep eps gp k p dp vs
    return (p',(p',vs',k+1))

-- | A 'Circuit' for classic gradient descent.
weightDecay
    :: Manifold x
    => Double -- ^ Decay Rate
    -> c # x -- ^ Point
    -> c # x -- ^ Weight-decayed point
{-# INLINE weightDecay #-}
weightDecay dcy = ((1-dcy) .>)



--- Internal ---


-- | A step of the basic xomentum algorithm.
momentumStep
    :: Manifold x
    => Double -- ^ The learning rate
    -> Double -- ^ The momentum decay
    -> c # x -- ^ The subsequent TangentPair
    -> c # x -- ^ The subsequent TangentPair
    -> c # x -- ^ The current velocity
    -> (c # x, c # x) -- ^ The (subsequent point, subsequent velocity)
{-# INLINE momentumStep #-}
momentumStep eps mu p fd v =
    let v' = eps .> fd <+> mu .> v
     in (gradientStep 1 p v', v')

-- | Note that we generally assume that momentum updates ignore the Riemannian metric.
adamStep
    :: Manifold x
    => Double -- ^ The learning rate
    -> Double -- ^ The first momentum rate
    -> Double -- ^ The second momentum rate
    -> Double -- ^ Second moment regularizer
    -> Int -- ^ Algorithm step
    -> c # x -- ^ The subsequent gradient
    -> c # x -- ^ The subsequent gradient
    -> c # x -- ^ First order velocity
    -> c # x -- ^ Second order velocity
    -> (c # x, c # x, c # x) -- ^ Subsequent (point, first velocity, second velocity)
{-# INLINE adamStep #-}
adamStep eps b1 b2 rg k0 p fd m v =
    let k = k0+1
        fd' = S.map (^(2 :: Int)) $ coordinates fd
        m' = (1-b1) .> fd <+> b1 .> m
        v' = (1-b2) .> Point fd' <+> b2 .> v
        mhat = (1-b1^k) /> m'
        vhat = (1-b2^k) /> v'
        fd'' = S.zipWith (/) (coordinates mhat) . S.map ((+ rg) . sqrt) $ coordinates vhat
     in (gradientStep eps p $ Point fd'', m',v')


---- | Gradient ascent based on the 'Riemannian' metric.
--gradientSequence
--    :: Riemannian c x
--    => Double -- ^ Step size
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent
--{-# INLINE gradientSequence #-}
--gradientSequence eps f = iterate (gradientStep' eps . sharp . f)
--
---- | Gradient ascent which ignores 'Riemannian' metric.
--vanillaGradientSequence
--    :: Manifold x
--    => Double -- ^ Step size
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent
--{-# INLINE vanillaGradientSequence #-}
--vanillaGradientSequence eps f = iterate (gradientStep' eps . breakPoint . f)
--
-- Momentum --
--
--
---- | Momentum ascent.
--momentumSequence :: Riemannian c x
--    => Double -- ^ Learning rate
--    -> (Int -> Double) -- ^ Momentum decay function
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent with momentum
--{-# INLINE momentumSequence #-}
--momentumSequence eps mu f p0 =
--    let v0 = zero
--        fd = sharp . f
--        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
--     in ps
--
---- | Vanilla Momentum ascent.
--vanillaMomentumSequence :: Manifold x
--    => Double -- ^ Learning rate
--    -> (Int -> Double) -- ^ Momentum decay function
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent with momentum
--{-# INLINE vanillaMomentumSequence #-}
--vanillaMomentumSequence eps mu f p0 =
--    let v0 = zero
--        fd = breakPoint . f
--        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
--     in ps
--
---- | Adam ascent.
--adamSequence :: Riemannian c x
--    => Double -- ^ The learning rate
--    -> Double -- ^ The first momentum rate
--    -> Double -- ^ The second momentum rate
--    -> Double -- ^ Second moment regularizer
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent with momentum
--{-# INLINE adamSequence #-}
--adamSequence eps b1 b2 rg f p0 =
--    let m0 = zero
--        v0 = zero
--        fd = sharp . f
--        (ps,_,_,_) = unzip4 $ iterate
--            (\(p,m,v,k) -> let (p',m',v') = adamStep eps b1 b2 rg k (fd p) m v in (p',m',v',k+1)) (p0,m0,v0,1)
--     in ps
--
---- | Vanilla Adam ascent.
--vanillaAdamSequence :: Manifold x
--    => Double -- ^ The learning rate
--    -> Double -- ^ The first momentum rate
--    -> Double -- ^ The second momentum rate
--    -> Double -- ^ Second moment regularizer
--    -> (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ The initial point
--    -> [c # x] -- ^ The gradient ascent with momentum
--{-# INLINE vanillaAdamSequence #-}
--vanillaAdamSequence eps b1 b2 rg f p0 =
--    let m0 = zero
--        v0 = zero
--        fd = breakPoint . f
--        (ps,_,_,_) = unzip4 $ iterate
--            (\(p,m,v,k) -> let (p',m',v') = adamStep eps b1 b2 rg k (fd p) m v in (p',m',v',k+1)) (p0,m0,v0,1)
--     in ps
--
-- Newton --

-- | A step of the Newton algorithm for nonlinear optimization.
--newtonStep
--    :: Manifold x
--    => c # x
--    -> CotangentPair c x -- ^ Derivatives
--    -> CotangentTensor c x -- ^ Hessian
--    -> c # x -- ^ Step
--{-# INLINE newtonStep #-}
--newtonStep p df ddf = gradientStep (-1) p $ inverse ddf >.> df
--
---- | An infinite list of iterations of the Newton algorithm for nonlinear optimization.
--newtonSequence
--    :: Manifold x
--    => (c # x -> CotangentPair c x)  -- ^ Gradient calculator
--    -> c # x -- ^ Initial point
--    -> [c # x] -- ^ Newton sequence
--{-# INLINE newtonSequence #-}
--newtonSequence f = iterate iterator
--    where iterator p = newtonStep p (differential f p) (hessian f p)


-- Gauss Newton --

{-
-- | A step of the Gauss-Newton algorithm for nonlinear optimization.
gaussNewtonStep
    :: (Manifold x, RealFrac x)
    => x -- ^ Damping factor
    -> Vector v k x -- ^ Residuals
    -> [CotangentVector c x] -- ^ Residual differentials
    -> c # x -- ^ Parameter estimates
gaussNewtonStep eps rs grds =
    gradientStep (-eps) $ linearLeastSquares0 (fromRows (Euclidean $ length grds) grds) rs

-- | An infinite list of iterations of the Gauss-Newton algorithm for nonlinear optimization.
gaussNewtonSequence :: (Manifold x, RealFrac x)
    => x -- ^ Damping Factor
    -> (c # x -> [x]) -- ^ Residual Function
    -> (c # x -> [Differentials :#: Tangent c x]) -- ^ Residual Differential
    -> (c # x) -- ^ Initial guess
    -> [c # x] -- ^ Gauss-Newton Sequence
gaussNewtonSequence dmp rsf rsf' = iterate iterator
  where iterator p = gaussNewtonStep dmp (rsf p) (rsf' p)
  -}

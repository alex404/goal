-- | Gradient pursuit based optimization on manifolds.

module Goal.Geometry.Differential.GradientPursuit
    ( -- * Cauchy Sequences
      cauchyLimit
    , cauchySequence
    -- * Gradient Pursuit
    , GradientPursuit (Classic,Momentum,Adam)
    , gradientPursuitStep
    , gradientSequence
    , vanillaGradientSequence
    , gradientCircuit
    -- ** Defaults
    , defaultMomentumPursuit
    , defaultAdamPursuit
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
    :: (Point c m -> Point c m -> Double) -- ^ Distance (divergence) from previous to next
    -> Double -- ^ Epsilon
    -> [Point c m] -- ^ Input sequence
    -> Point c m
{-# INLINE cauchyLimit #-}
cauchyLimit f eps ps = last $ cauchySequence f eps ps

-- | Attempts to calculate the limit of a sequence. Returns the list up to the limit.
cauchySequence
    :: (Point c m -> Point c m -> Double) -- ^ Distance (divergence) from previous to next
    -> Double -- ^ Epsilon
    -> [Point c m] -- ^ Input list
    -> [Point c m] -- ^ Truncated list
{-# INLINE cauchySequence #-}
cauchySequence f eps ps =
    let pps = takeWhile taker . zip ps $ tail ps
     in head ps : fmap snd pps
       where taker (p1,p2) = eps < f p1 p2


--- Gradient Pursuit ---


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
    :: Manifold m
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Int -- ^ Algorithm step
    -> TangentPair c m -- ^ The subsequent TangentPair
    -> [TangentVector c m] -- ^ The velocities
    -> (Point c m, [TangentVector c m]) -- ^ The updated point and velocities
gradientPursuitStep eps Classic _ dp _ = (gradientStep' eps dp,[])
gradientPursuitStep eps (Momentum fmu) k dp (v:_) =
    let (p,v') = momentumStep eps (fmu k) dp v
     in (p,[v'])
gradientPursuitStep eps (Adam b1 b2 rg) k dp (m:v:_) =
    let (p,m',v') = adamStep eps b1 b2 rg k dp m v
     in (p,[m',v'])
gradientPursuitStep _ _ _ _ _ = error "Momentum list length mismatch in gradientPursuitStep"


-- | Gradient ascent based on the 'Riemannian' metric.
gradientSequence
    :: Riemannian c m
    => (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
    -> Double -- ^ Step size
    -> GradientPursuit  -- ^ Gradient pursuit algorithm
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent
{-# INLINE gradientSequence #-}
gradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0,(repeat zero,0))
        where iterator (p,(vs,k)) =
                  let dp = sharp $ f p
                      (p',vs') = gradientPursuitStep eps gp k dp vs
                   in (p',(vs',k+1))

-- | Gradient ascent based on the 'Riemannian' metric.
vanillaGradientSequence
    :: Manifold m
    => (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
    -> Double -- ^ Step size
    -> GradientPursuit  -- ^ Gradient pursuit algorithm
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent
{-# INLINE vanillaGradientSequence #-}
vanillaGradientSequence f eps gp p0 =
    fst <$> iterate iterator (p0,(repeat zero,0))
        where iterator (p,(vs,k)) =
                  let dp = vanillaGradient' $ f p
                      (p',vs') = gradientPursuitStep eps gp k dp vs
                   in (p',(vs',k+1))

-- | A 'Circuit' for classic gradient descent.
gradientCircuit
    :: (Monad m, Manifold z)
    => Double -- ^ Learning Rate
    -> GradientPursuit -- ^ Gradient pursuit algorithm
    -> Circuit m (TangentPair c z) (Point c z) -- ^ Gradient Ascent
{-# INLINE gradientCircuit #-}
gradientCircuit eps gp = accumulateFunction (repeat zero,0) $ \pdp (vs,k) -> do
    let (p',vs') = gradientPursuitStep eps gp k pdp vs
    return (p',(vs',k+1))


--- Internal ---


-- | A step of the basic momentum algorithm.
momentumStep
    :: Manifold m
    => Double -- ^ The learning rate
    -> Double -- ^ The momentum decay
    -> TangentPair c m -- ^ The subsequent TangentPair
    -> TangentVector c m -- ^ The current velocity
    -> (Point c m, TangentVector c m) -- ^ The (subsequent point, subsequent velocity)
{-# INLINE momentumStep #-}
momentumStep eps mu pfd v =
    let (p,fd) = splitTangentPair pfd
        v' = eps .> fd <+> mu .> v
     in (gradientStep 1 p v', v')

-- | Note that we generally assume that momentum updates ignore the Riemannian metric.
adamStep
    :: Manifold m
    => Double -- ^ The learning rate
    -> Double -- ^ The first momentum rate
    -> Double -- ^ The second momentum rate
    -> Double -- ^ Second moment regularizer
    -> Int -- ^ Algorithm step
    -> TangentPair c m -- ^ The subsequent gradient
    -> TangentVector c m -- ^ First order velocity
    -> TangentVector c m -- ^ Second order velocity
    -> (Point c m, TangentVector c m, TangentVector c m) -- ^ Subsequent (point, first velocity, second velocity)
{-# INLINE adamStep #-}
adamStep eps b1 b2 rg k0 pfd m v =
    let k = k0+1
        (p,fd) = splitTangentPair pfd
        fd' = S.map (^(2 :: Int)) $ coordinates fd
        m' = (1-b1) .> fd <+> b1 .> m
        v' = (1-b2) .> Point fd' <+> b2 .> v
        mhat = (1-b1^k) /> m'
        vhat = (1-b2^k) /> v'
        fd'' = S.zipWith (/) (coordinates mhat) . S.map ((+ rg) . sqrt) $ coordinates vhat
     in (gradientStep eps p $ Point fd'', m',v')


---- | Gradient ascent based on the 'Riemannian' metric.
--gradientSequence
--    :: Riemannian c m
--    => Double -- ^ Step size
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent
--{-# INLINE gradientSequence #-}
--gradientSequence eps f = iterate (gradientStep' eps . sharp . f)
--
---- | Gradient ascent which ignores 'Riemannian' metric.
--vanillaGradientSequence
--    :: Manifold m
--    => Double -- ^ Step size
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent
--{-# INLINE vanillaGradientSequence #-}
--vanillaGradientSequence eps f = iterate (gradientStep' eps . breakPoint . f)
--
-- Momentum --
--
--
---- | Momentum ascent.
--momentumSequence :: Riemannian c m
--    => Double -- ^ Learning rate
--    -> (Int -> Double) -- ^ Momentum decay function
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent with momentum
--{-# INLINE momentumSequence #-}
--momentumSequence eps mu f p0 =
--    let v0 = zero
--        fd = sharp . f
--        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
--     in ps
--
---- | Vanilla Momentum ascent.
--vanillaMomentumSequence :: Manifold m
--    => Double -- ^ Learning rate
--    -> (Int -> Double) -- ^ Momentum decay function
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent with momentum
--{-# INLINE vanillaMomentumSequence #-}
--vanillaMomentumSequence eps mu f p0 =
--    let v0 = zero
--        fd = breakPoint . f
--        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
--     in ps
--
---- | Adam ascent.
--adamSequence :: Riemannian c m
--    => Double -- ^ The learning rate
--    -> Double -- ^ The first momentum rate
--    -> Double -- ^ The second momentum rate
--    -> Double -- ^ Second moment regularizer
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent with momentum
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
--vanillaAdamSequence :: Manifold m
--    => Double -- ^ The learning rate
--    -> Double -- ^ The first momentum rate
--    -> Double -- ^ The second momentum rate
--    -> Double -- ^ Second moment regularizer
--    -> (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ The initial point
--    -> [Point c m] -- ^ The gradient ascent with momentum
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
--    :: Manifold m
--    => Point c m
--    -> CotangentPair c m -- ^ Derivatives
--    -> CotangentTensor c m -- ^ Hessian
--    -> Point c m -- ^ Step
--{-# INLINE newtonStep #-}
--newtonStep p df ddf = gradientStep (-1) p $ inverse ddf >.> df
--
---- | An infinite list of iterations of the Newton algorithm for nonlinear optimization.
--newtonSequence
--    :: Manifold m
--    => (Point c m -> CotangentPair c m)  -- ^ Gradient calculator
--    -> Point c m -- ^ Initial point
--    -> [Point c m] -- ^ Newton sequence
--{-# INLINE newtonSequence #-}
--newtonSequence f = iterate iterator
--    where iterator p = newtonStep p (differential f p) (hessian f p)


-- Gauss Newton --

{-
-- | A step of the Gauss-Newton algorithm for nonlinear optimization.
gaussNewtonStep
    :: (Manifold m, RealFrac x)
    => x -- ^ Damping factor
    -> Vector v k x -- ^ Residuals
    -> [CotangentVector c m] -- ^ Residual differentials
    -> Point c m -- ^ Parameter estimates
gaussNewtonStep eps rs grds =
    gradientStep (-eps) $ linearLeastSquares0 (fromRows (Euclidean $ length grds) grds) rs

-- | An infinite list of iterations of the Gauss-Newton algorithm for nonlinear optimization.
gaussNewtonSequence :: (Manifold m, RealFrac x)
    => x -- ^ Damping Factor
    -> (Point c m -> [x]) -- ^ Residual Function
    -> (Point c m -> [Differentials :#: Tangent c m]) -- ^ Residual Differential
    -> (Point c m) -- ^ Initial guess
    -> [Point c m] -- ^ Gauss-Newton Sequence
gaussNewtonSequence dmp rsf rsf' = iterate iterator
  where iterator p = gaussNewtonStep dmp (rsf p) (rsf' p)
  -}

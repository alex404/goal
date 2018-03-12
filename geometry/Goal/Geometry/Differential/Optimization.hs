-- | Provides a few general tools and algorithms for numerical optimization.

module Goal.Geometry.Differential.Optimization (
    -- * Cauchy Sequences
      cauchyLimit
    , cauchySequence
    -- * Gradient Pursuit
    , gradientSequence
    , vanillaGradientSequence
    -- ** Momentum
    , momentumStep
    , defaultMomentumSchedule
    , momentumSequence
    , vanillaMomentumSequence
    -- ** Adam
    , adamStep
    , adamSequence
    , vanillaAdamSequence
    -- * Least Squares
    , linearLeastSquares
    , linearLeastSquares0
    -- ** Newton
    , newtonStep
    , newtonSequence
    ) where


--- Imports ---


import Data.List (unzip4)

-- Goal --

import Goal.Core

import Goal.Geometry.Manifold
import Goal.Geometry.Linear
import Goal.Geometry.Map
import Goal.Geometry.Map.Multilinear
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S

--- Cauchy Sequences ---


-- | Attempts to calculate the limit of a sequence. This finds the iterate with a sufficiently small
-- distance from the previous iterate.
cauchyLimit
    :: Ord x
    => (Point c m -> Point c m -> x) -- ^ Distance (divergence) from previous to next
    -> x -- ^ Epsilon
    -> [Point c m] -- ^ Input sequence
    -> Point c m
{-# INLINE cauchyLimit #-}
cauchyLimit f eps ps = last $ cauchySequence f eps ps

-- | Attempts to calculate the limit of a sequence. Returns the list up to the limit.
cauchySequence
    :: Ord x
    => (Point c m -> Point c m -> x) -- ^ Distance (divergence) from previous to next
    -> x -- ^ Epsilon
    -> [Point c m] -- ^ Input list
    -> [Point c m] -- ^ Truncated list
{-# INLINE cauchySequence #-}
cauchySequence f eps ps =
    let pps = takeWhile taker . zip ps $ tail ps
     in head ps : fmap snd pps
       where taker (p1,p2) = eps < f p1 p2


--- Gradient Pursuit ---


-- | Gradient ascent based on the 'Riemannian' metric.
gradientSequence
    :: Riemannian c m
    => Double -- ^ Step size
    -> (forall x. RealFloat x =>  B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent
{-# INLINE gradientSequence #-}
gradientSequence eps f = iterate (gradientStep' eps . sharp . differential' f)

-- | Gradient ascent which ignores 'Riemannian' metric.
vanillaGradientSequence
    :: Manifold m
    => Double -- ^ Step size
    -> (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent
{-# INLINE vanillaGradientSequence #-}
vanillaGradientSequence eps f = iterate (gradientStep' eps . breakChart . differential' f)

-- Momentum --

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

defaultMomentumSchedule :: RealFloat x => x -> Int -> x
{-# INLINE defaultMomentumSchedule #-}
defaultMomentumSchedule mxmu k = min mxmu $ 1 - 2**((negate 1 -) . logBase 2 . fromIntegral $ div k 250 + 1)

-- | Momentum ascent.
momentumSequence :: Riemannian c m
    => Double -- ^ Learning rate
    -> (Int -> Double) -- ^ Momentum decay function
    -> (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent with momentum
{-# INLINE momentumSequence #-}
momentumSequence eps mu f p0 =
    let v0 = zero
        fd = sharp . differential' f
        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
     in ps

-- | Vanilla Momentum ascent.
vanillaMomentumSequence :: Manifold m
    => Double -- ^ Learning rate
    -> (Int -> Double) -- ^ Momentum decay function
    -> (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent with momentum
{-# INLINE vanillaMomentumSequence #-}
vanillaMomentumSequence eps mu f p0 =
    let v0 = zero
        fd = breakChart . differential' f
        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
     in ps

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
adamStep eps b1 b2 rg k pfd m v =
    let (p,fd) = splitTangentPair pfd
        fd' = G.map (^(2 :: Int)) $ coordinates fd
        m' = (1-b1) .> fd <+> b1 .> m
        v' = (1-b2) .> Point fd' <+> b2 .> v
        mhat = (1-b1^k) /> m'
        vhat = (1-b2^k) /> v'
        fd'' = G.zipWith (/) (coordinates mhat) . G.map ((+ rg) . sqrt) $ coordinates vhat
     in (gradientStep eps p $ Point fd'', m',v')

-- | Adam ascent.
adamSequence :: Riemannian c m
    => Double -- ^ The learning rate
    -> Double -- ^ The first momentum rate
    -> Double -- ^ The second momentum rate
    -> Double -- ^ Second moment regularizer
    -> (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent with momentum
{-# INLINE adamSequence #-}
adamSequence eps b1 b2 rg f p0 =
    let m0 = zero
        v0 = zero
        fd = sharp . differential' f
        (ps,_,_,_) = unzip4 $ iterate
            (\(p,m,v,k) -> let (p',m',v') = adamStep eps b1 b2 rg k (fd p) m v in (p',m',v',k+1)) (p0,m0,v0,1)
     in ps

-- | Vanilla Adam ascent.
vanillaAdamSequence :: Manifold m
    => Double -- ^ The learning rate
    -> Double -- ^ The first momentum rate
    -> Double -- ^ The second momentum rate
    -> Double -- ^ Second moment regularizer
    -> (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ The initial point
    -> [Point c m] -- ^ The gradient ascent with momentum
{-# INLINE vanillaAdamSequence #-}
vanillaAdamSequence eps b1 b2 rg f p0 =
    let m0 = zero
        v0 = zero
        fd = breakChart . differential' f
        (ps,_,_,_) = unzip4 $ iterate
            (\(p,m,v,k) -> let (p',m',v') = adamStep eps b1 b2 rg k (fd p) m v in (p',m',v',k+1)) (p0,m0,v0,1)
     in ps


--- Least Squares ---

-- | Linear least squares estimation.
linearLeastSquares
    :: (Manifold m, KnownNat k, 1 <= k)
    => S.Vector k (Point c m) -- ^ Independent variable observations
    -> S.Vector k Double -- ^ Dependent variable observations
    -> Point (Dual c) m -- ^ Parameter estimates
{-# INLINE linearLeastSquares #-}
linearLeastSquares xs ys =
    let mtx = S.fromRows $ S.map coordinates xs
     in linearLeastSquares0 mtx ys

-- | Linear least squares estimation, where the design matrix is provided directly.
linearLeastSquares0
    :: (Manifold m, KnownNat k)
    => S.Matrix k (Dimension m) Double -- ^ Independent variable observations
    -> S.Vector k Double -- ^ Dependent variable observations
    -> Point c m -- ^ Parameter estimates
{-# INLINE linearLeastSquares0 #-}
linearLeastSquares0 mtx ys =
    let tmtx = S.transpose mtx
        prj = S.matrixMatrixMultiply (S.inverse $ S.matrixMatrixMultiply tmtx mtx) tmtx
     in Point $ S.matrixVectorMultiply prj ys

-- Newton --

-- | A step of the Newton algorithm for nonlinear optimization.
newtonStep
    :: Manifold m
    => Point c m
    -> CotangentVector c m -- ^ Derivatives
    -> CotangentTensor c m -- ^ Hessian
    -> Point c m -- ^ Step
{-# INLINE newtonStep #-}
newtonStep p df ddf = gradientStep (-1) p $ inverse ddf >.> df

-- | An infinite list of iterations of the Newton algorithm for nonlinear optimization.
newtonSequence
    :: Manifold m
    => (forall x. RealFloat x => B.Vector (Dimension m) x -> x)  -- ^ Function to minimize
    -> Point c m -- ^ Initial point
    -> [Point c m] -- ^ Newton sequence
{-# INLINE newtonSequence #-}
newtonSequence f = iterate iterator
    where iterator p = newtonStep p (differential f p) (hessian f p)


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

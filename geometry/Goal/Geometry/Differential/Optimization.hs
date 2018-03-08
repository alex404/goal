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

--- Cauchy Sequences ---


-- | Attempts to calculate the limit of a sequence. This finds the iterate with a sufficiently small
-- distance from the previous iterate.
cauchyLimit
    :: Ord x
    => (Point v c m x -> Point v c m x -> x) -- ^ Distance (divergence) from previous to next
    -> x -- ^ Epsilon
    -> [Point v c m x] -- ^ Input sequence
    -> Point v c m x
{-# INLINE cauchyLimit #-}
cauchyLimit f eps ps = last $ cauchySequence f eps ps

-- | Attempts to calculate the limit of a sequence. Returns the list up to the limit.
cauchySequence
    :: Ord x
    => (Point v c m x -> Point v c m x -> x) -- ^ Distance (divergence) from previous to next
    -> x -- ^ Epsilon
    -> [Point v c m x] -- ^ Input list
    -> [Point v c m x] -- ^ Truncated list
{-# INLINE cauchySequence #-}
cauchySequence f eps ps =
    let pps = takeWhile taker . zip ps $ tail ps
     in head ps : fmap snd pps
       where taker (p1,p2) = eps < f p1 p2


--- Gradient Pursuit ---


-- | Gradient ascent based on the 'Riemannian' metric.
gradientSequence
    :: (Riemannian v c m x, GVector v x, RealFloat x)
    => x -- ^ Step size
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent
{-# INLINE gradientSequence #-}
gradientSequence eps f = iterate (gradientStep' eps . sharp . differential' f)

-- | Gradient ascent which ignores 'Riemannian' metric.
vanillaGradientSequence
    :: (Manifold m, GVector v x, RealFloat x)
    => x -- ^ Step size
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent
{-# INLINE vanillaGradientSequence #-}
vanillaGradientSequence eps f = iterate (gradientStep' eps . breakChart . differential' f)

-- Momentum --

-- | A step of the basic momentum algorithm.
momentumStep
    :: (Manifold m, GVector v x, RealFloat x)
    => x -- ^ The learning rate
    -> x -- ^ The momentum decay
    -> TangentPair v c m x -- ^ The subsequent TangentPair
    -> TangentVector v c m x -- ^ The current velocity
    -> (Point v c m x, TangentVector v c m x) -- ^ The (subsequent point, subsequent velocity)
{-# INLINE momentumStep #-}
momentumStep eps mu pfd v =
    let (p,fd) = splitTangentPair pfd
        v' = eps .> fd <+> mu .> v
     in (gradientStep 1 p v', v')

defaultMomentumSchedule :: RealFloat x => x -> Int -> x
{-# INLINE defaultMomentumSchedule #-}
defaultMomentumSchedule mxmu k = min mxmu $ 1 - 2**((negate 1 -) . logBase 2 . fromIntegral $ div k 250 + 1)

-- | Momentum ascent.
momentumSequence :: (Riemannian v c m x, GVector v x, RealFloat x)
    => x -- ^ Learning rate
    -> (Int -> x) -- ^ Momentum decay function
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent with momentum
{-# INLINE momentumSequence #-}
momentumSequence eps mu f p0 =
    let v0 = zero
        fd = sharp . differential' f
        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
     in ps

-- | Vanilla Momentum ascent.
vanillaMomentumSequence :: (Manifold m, GVector v x, RealFloat x)
    => x -- ^ Learning rate
    -> (Int -> x) -- ^ Momentum decay function
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent with momentum
{-# INLINE vanillaMomentumSequence #-}
vanillaMomentumSequence eps mu f p0 =
    let v0 = zero
        fd = breakChart . differential' f
        (ps,_,_) = unzip3 $ iterate (\(p,v,k) -> let (p',v') = momentumStep eps (mu k) (fd p) v in (p',v',k+1)) (p0,v0,0)
     in ps

-- | Note that we generally assume that momentum updates ignore the Riemannian metric.
adamStep
    :: (Manifold m, GVector v x, RealFloat x)
    => x -- ^ The learning rate
    -> x -- ^ The first momentum rate
    -> x -- ^ The second momentum rate
    -> x -- ^ Second moment regularizer
    -> Int -- ^ Algorithm step
    -> TangentPair v c m x -- ^ The subsequent gradient
    -> TangentVector v c m x -- ^ First order velocity
    -> TangentVector v c m x -- ^ Second order velocity
    -> (Point v c m x, TangentVector v c m x, TangentVector v c m x) -- ^ Subsequent (point, first velocity, second velocity)
{-# INLINE adamStep #-}
adamStep eps b1 b2 rg k pfd m v =
    let (p,fd) = splitTangentPair pfd
        fd' = G.map (^(2 :: Int)) fd
        m' = (1-b1) .> fd <+> b1 .> m
        v' = (1-b2) .> fd' <+> b2 .> v
        mhat = (1-b1^k) /> m'
        vhat = (1-b2^k) /> v'
        fd'' = G.zipWith (/) (coordinates mhat) $ (+ rg) . sqrt <$> coordinates vhat
     in (gradientStep eps p $ Point fd'', m',v')

-- | Adam ascent.
adamSequence :: (Riemannian v c m x, GVector v x, RealFloat x)
    => x -- ^ The learning rate
    -> x -- ^ The first momentum rate
    -> x -- ^ The second momentum rate
    -> x -- ^ Second moment regularizer
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent with momentum
{-# INLINE adamSequence #-}
adamSequence eps b1 b2 rg f p0 =
    let m0 = zero
        v0 = zero
        fd = sharp . differential' f
        (ps,_,_,_) = unzip4 $ iterate
            (\(p,m,v,k) -> let (p',m',v') = adamStep eps b1 b2 rg k (fd p) m v in (p',m',v',k+1)) (p0,m0,v0,1)
     in ps

-- | Vanilla Adam ascent.
vanillaAdamSequence :: (Manifold m, GVector v x, RealFloat x)
    => x -- ^ The learning rate
    -> x -- ^ The first momentum rate
    -> x -- ^ The second momentum rate
    -> x -- ^ Second moment regularizer
    -> (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ The initial point
    -> [Point v c m x] -- ^ The gradient ascent with momentum
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
    :: (Manifold m, KnownNat k, 1 <= k, GVector v x, RealFloat x)
    => Vector v k (Point v c m x) -- ^ Independent variable observations
    -> Vector v k x -- ^ Dependent variable observations
    -> Point v (Dual c) m x -- ^ Parameter estimates
{-# INLINE linearLeastSquares #-}
linearLeastSquares xs ys =
    let mtx = fromMatrix . fromRows $ coordinates <$> xs
     in linearLeastSquares0 mtx ys

-- | Linear least squares estimation, where the design matrix is provided directly.
linearLeastSquares0
    :: (Manifold m, KnownNat k, Storable x)
    => Point v (Function c Cartesian) (Product (Euclidean k) m) x -- ^ Design matrix
    -> Vector SVector k x -- ^ Independent variables
    -> Point v c m x -- ^ Parameter estimates
{-# INLINE linearLeastSquares0 #-}
linearLeastSquares0 mtx ys =
    let tmtx = transpose mtx
        prj = (fromJust . inverse $ tmtx <#> mtx) <#> tmtx
     in prj >.> Point ys

-- Newton --

-- | A step of the Newton algorithm for nonlinear optimization.
newtonStep
    :: (Manifold m, GVector v x, RealFloat x)
    => Point v c m x
    -> CotangentVector v c m x -- ^ Derivatives
    -> CotangentTensor v c m x -- ^ Hessian
    -> Point v c m x -- ^ Step
{-# INLINE newtonStep #-}
newtonStep p df ddf = gradientStep (-1) p $ fromJust (inverse ddf) >.> df

-- | An infinite list of iterations of the Newton algorithm for nonlinear optimization.
newtonSequence
    :: (Manifold m, GVector v x)
    => (forall z. RealFloat z => BPoint c m z -> z)  -- ^ Function to minimize
    -> Point v c m x -- ^ Initial point
    -> [Point v c m x] -- ^ Newton sequence
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
    -> [CotangentVector v c m x] -- ^ Residual differentials
    -> Point v c m x -- ^ Parameter estimates
gaussNewtonStep eps rs grds =
    gradientStep (-eps) $ linearLeastSquares0 (fromRows (Euclidean $ length grds) grds) rs

-- | An infinite list of iterations of the Gauss-Newton algorithm for nonlinear optimization.
gaussNewtonSequence :: (Manifold m, RealFrac x)
    => x -- ^ Damping Factor
    -> (Point v c m x -> [x]) -- ^ Residual Function
    -> (Point v c m x -> [Differentials :#: Tangent c m]) -- ^ Residual Differential
    -> (Point v c m x) -- ^ Initial guess
    -> [Point v c m x] -- ^ Gauss-Newton Sequence
gaussNewtonSequence dmp rsf rsf' = iterate iterator
  where iterator p = gaussNewtonStep dmp (rsf p) (rsf' p)
  -}

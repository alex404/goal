{-# LANGUAGE
    RankNTypes,
    TypeOperators,
    TypeApplications,
    FlexibleContexts,
    ScopedTypeVariables
#-}
-- | The main module of goal-probability. Import this module to use all the features provided by this library.
module Goal.Probability
    ( -- * Package Exports
      module Goal.Probability.Statistical
    , module Goal.Probability.ExponentialFamily
    , module Goal.Probability.Distributions
    , module Goal.Probability.ExponentialFamily.NeuralNetwork
    , module Goal.Probability.ExponentialFamily.PopulationCode
    , module Goal.Probability.ExponentialFamily.Rectification
    , module Goal.Probability.ExponentialFamily.Harmonium
    , module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    , module Goal.Probability.ExponentialFamily.Harmonium.Optimization
      -- * Utility
    , resampleVector
    , subsampleVector
    , shuffleList
    , noisyFunction
    , estimateCorrelations
    , seed
    -- * External Exports
    , module System.Random.MWC
    , module System.Random.MWC.Probability
    ) where


--- Imports ---


-- Re-exports --

import System.Random.MWC (Seed,save,restore)
import qualified System.Random.MWC as MWC

import System.Random.MWC.Probability hiding (initialize,sample)
import System.Random.MWC.Distributions (uniformShuffle)

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily.NeuralNetwork
import Goal.Probability.ExponentialFamily.PopulationCode
import Goal.Probability.ExponentialFamily.Rectification
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Conditional
import Goal.Probability.ExponentialFamily.Harmonium.Optimization


-- Package --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic.Mutable as M
import qualified Goal.Core.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable.Base as MV
import qualified Data.Vector as V

--- Stochastic Functions ---

-- | Creates a seed for later RandST usage.
seed :: Random s Seed
{-# INLINE seed #-}
seed = Prob save

shuffleList :: [a] -> Random r [a]
shuffleList xs = fmap V.toList . Prob $ uniformShuffle (V.fromList xs)

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

subsampleVector
    :: forall k m x r . (KnownNat k, KnownNat m)
    => B.Vector (k + m) x
    -> Random r (B.Vector k x)
{-# INLINE subsampleVector #-}
subsampleVector v = Prob $ \gn -> do
    let k = natValInt (Proxy :: Proxy k)
    mv <- G.thaw v
    randomSubSample0 k mv gn
    v' <- G.unsafeFreeze mv
    let foo :: (B.Vector k x, B.Vector m x)
        foo = B.splitAt v'
    return $ fst foo

randomSubSample0
    :: (KnownNat n, PrimMonad m, MV.MVector v a)
    => Int -> G.MVector v n (PrimState m) a -> Gen (PrimState m) -> m ()
{-# INLINE randomSubSample0 #-}
randomSubSample0 k v gn = looper 0
    where n = M.length v
          looper i
            | i == k = return ()
            | otherwise = do
                j <- MWC.uniformR (i,n-1) gn
                M.unsafeSwap v i j
                looper (i+1)

pop :: Int -> [x] -> (x,[x])
pop idx xs = (x,lft ++ rgt)
  where (lft, (x:rgt)) = splitAt idx xs

estimateCorrelations
    :: forall k x v . (G.VectorClass v x, G.VectorClass v Double, KnownNat k, Real x)
    => [G.Vector v k x]
    -> [Double]
{-# INLINE estimateCorrelations #-}
estimateCorrelations zs =
    let mnrm :: Source # MultivariateNormal k
        mnrm = mle $ G.convert . G.map realToFrac <$> zs
        k = natValInt $ Proxy @ k
        lwrs = drop k $ listCoordinates mnrm
        trngs = subtract 1 . S.triangularNumber <$> [1..k]
        folder trng (vrs',cvrs') = let (vr,cvrs'') = pop trng cvrs' in (vr:vrs',cvrs'')
        (vrs,cvrs) = foldr folder ([],lwrs) trngs
        sds = sqrt <$> vrs
        dvs = concat $ do
            i <- [2..k]
            let sbsds = take i sds
                (rwsds,clsd) = splitAt (i-1) sbsds
            return $ (head clsd *) <$> rwsds
     in zipWith (/) cvrs dvs



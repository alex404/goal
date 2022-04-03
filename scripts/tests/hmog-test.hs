#! stack runghc

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S

--- Globals ---

strngx :: Source # MVNCovariance (MVNMean 3) (MVNMean 3)
strngx = fromTuple (2,-0.5,2,-0.5,-0.3,1.5)

snrmx :: Source # FullNormal 3
snrmx = join 1 strngx

strngz :: Source # MVNCovariance (MVNMean 2) (MVNMean 2)
strngz = fromTuple (2,0.5,3)

snrmz :: Source # FullNormal 2
snrmz = join (-1) strngz

stnsxz :: Source # Tensor (MVNMean 3) (MVNMean 2)
stnsxz = fromTuple (-0.1,1,0.5,0.2,-0.5,-1)

slgh :: Source # FullGaussianHarmonium 3 2
slgh = join (join snrmx stnsxz) snrmz

nsmps :: Int
nsmps = 100000

-- | Converts a point on a 'Tensor manifold into a Matrix.
toMatrix :: (Manifold x, Manifold y) => c # Tensor y x -> S.Matrix (Dimension y) (Dimension x) Double
{-# INLINE toMatrix #-}
toMatrix (Point xs) = G.Matrix xs

-- | Converts a Matrix into a 'Point' on a 'Tensor 'Manifold'.
fromMatrix :: S.Matrix (Dimension y) (Dimension x) Double -> c # Tensor y x
{-# INLINE fromMatrix #-}
fromMatrix (G.Matrix xs) = Point xs


printBig bgmtx = do
    S.mapM (print . S.toList) . S.toRows $ bgmtx

printBlocks tl tr br = do
    mapM_ print $ zipWith (++)
        (map listCoordinates . S.toList . toRows $ toTensor tl)
        (map listCoordinates . S.toList . toRows $ tr )
    mapM_ print $ zipWith (++)
        (map listCoordinates . S.toList . toRows $ transpose tr)
        (map listCoordinates . S.toList . toRows $ toTensor br)


--- Main ---


main :: IO ()
main = do
    let nlgh :: Natural # FullGaussianHarmonium 3 2
        nlgh = toNatural slgh

        mlgh :: Mean # FullGaussianHarmonium 3 2
        mlgh = toMean nlgh

    putStrLn "Isotransform (Source - Natural) error:"
    print . euclideanDistance slgh . toSource $ toNatural slgh

    putStrLn "Isotransform (Source - Mean) error:"
    print . euclideanDistance slgh . toSource $ toMean slgh

    xzs <- realize $ sample nsmps nlgh
    let (xs,zs) = unzip xzs
    let nprr = snd $ splitConjugatedHarmonium nlgh
    --zs' <- realize $ sample nsmps nprr

    let mlgh' :: Mean # FullGaussianHarmonium 3 2
        mlgh' = averageSufficientStatistic xzs

    let (maff,mprr) = split mlgh
        (maff',mprr') = split mlgh'

    putStrLn "Sampling error:"
    --print . euclideanDistance mprr' $ averageSufficientStatistic zs'
    --print $ euclideanDistance mlgh mlgh'

    let cvr0 :: S.Matrix 5 5 Double
        cvr0 =
            let top = S.horizontalConcat (toMatrix $ toTensor strngx) (toMatrix stnsxz)
                btm = S.horizontalConcat (S.transpose $ toMatrix stnsxz) (toMatrix $ toTensor strngz)
             in S.verticalConcat top btm

    let cvr = fromTensor $ fromMatrix cvr0
        smvn :: Source # FullNormal 5
        smvn = join (fromTuple (1,1,1,-1,-1)) cvr


    xzs0' <- realize $ sample nsmps smvn
    let xs = fst <$> xzs
        xzs' = S.splitAt <$> xzs0'
        (xs',zs') = unzip xzs'
        xmrg0,xmrg',xmrg'' :: Mean # FullNormal 3
        xmrg0 = toMean snrmx
        xmrg' = averageSufficientStatistic xs
        xmrg'' = averageSufficientStatistic xs'
        xmrg = snd . split $ transposeHarmonium mlgh

    print $ euclideanDistance xmrg0 xmrg
    print . euclideanDistance (averageSufficientStatistic zs) . transition . snd $ splitConjugatedHarmonium nlgh
    print . euclideanDistance (averageSufficientStatistic zs') . transition . snd $ splitConjugatedHarmonium nlgh
    --print $ euclideanDistance smvn $ mle xzs0'
    --print $ euclideanDistance xmrg' xmrg
    --print $ euclideanDistance xmrg'' xmrg
    --print $ euclideanDistance xmrg'' xmrg'
    --print $ euclideanDistance mlgh $ averageSufficientStatistic xzs
    --print $ euclideanDistance mlgh $ averageSufficientStatistic xzs'

    ----print bigMatrix
    --putStrLn "Simple Inverse"
    --printBig $ S.inverse bigMatrix
    --putStrLn "Advanced Inverse"
    --printBlocks nvrx0 nvrxz0 nvrz0
    --print nvrz


#! stack runghc

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE ScopedTypeVariables,DataKinds,TypeOperators,TypeFamilies #-}
--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S

--- Globals ---

type XCovar = Diagonal
--type XCovar = Diagonal
--type XCovar = MVNCovariance
type LGH = LinearGaussianHarmonium XCovar 3 2

strngx :: Source # XCovar (MVNMean 3) (MVNMean 3)
--strngx = fromTuple (2,-0.5,2,-0.5,-0.3,1.5)
strngx = fromTensor $ toTensor strngx0

strngx0 :: Source # Diagonal (MVNMean 3) (MVNMean 3)
strngx0 = fromTuple (2,2,1.5)

snrmx :: Source # MultivariateNormal XCovar 3
snrmx = join 1 strngx

strngz :: Source # MVNCovariance (MVNMean 2) (MVNMean 2)
strngz = fromTuple (2,0.5,3)
--strngz = fromTuple (1,0,1)

snrmz :: Source # FullNormal 2
snrmz = join (-1) strngz
--snrmz = join 0 strngz

stnsxz :: Source # Tensor (MVNMean 3) (MVNMean 2)
stnsxz = fromTuple (-0.1,1,0.5,0.2,-0.5,-1)
--stnsxz = 0
--stnsxz = fromTuple (1,0,1,0,0,1)

slgh :: Source # LGH
slgh = join (join snrmx stnsxz) snrmz

vr0 :: Double
vr0 = 2

prr0,prr1,prr2 :: Source # FullNormal 2
prr0 = fromTuple (0,0,vr0,0,vr0)
prr1 = fromTuple (0.0,0,vr0,-0.01,vr0)
prr2 = fromTuple (0,0,vr0,0.01,vr0)

nprr0,nprr1,nprr2 :: Natural # FullNormal 2
nprr0 = transition prr0
nprr1 = transition prr1
nprr2 = transition prr2

mog0 :: Natural # Mixture (FullNormal 2) 2
mog0 = joinNaturalMixture (S.fromTuple (nprr0,nprr1,nprr2)) 0

hmog :: Natural # HierarchicalMixtureOfGaussians XCovar 3 2 2
hmog = join (fst . split $ toNatural slgh) mog0

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

toDiagonalLGH
    ::  Natural # FullGaussianHarmonium 3 2
    ->  Natural # DiagonalGaussianHarmonium 3 2
toDiagonalLGH nlgh =
    let (nzx,nx) = split $ transposeHarmonium nlgh
        (nmu,nvr) = split nx
        nvr' :: Natural # Diagonal (MVNMean 3) (MVNMean 3)
        nvr' = fromTensor $ toTensor nvr
     in transposeHarmonium . join nzx $ join nmu nvr'

comparison p q = do
    let xs = listCoordinates p
        ys = listCoordinates q
        strs = [ concat [ "(", showFFloat (Just 4) x ", "
                        , showFFloat (Just 4) y ", "
                        , showFFloat (Just 4) (x-y) ", " ]
                        | (x,y) <- zip xs ys ]
    mapM_ putStrLn strs


--- Main ---


main :: IO ()
main = do
    let nlgh :: Natural # LGH
        nlgh = toNatural slgh

        mlgh :: Mean # LGH
        mlgh = toMean nlgh

    let (nzx,nx) = split $ transposeHarmonium nlgh
        (nmu,nvr) = split nx
        nvr' :: Natural # Diagonal (MVNMean 3) (MVNMean 3)
        nvr' = fromTensor $ toTensor nvr

        nlgh' :: Natural # LGH
        nlgh' = transposeHarmonium . join nzx . join nmu . fromTensor $ toTensor nvr'

        stnd :: Source # FullNormal 2
        stnd = standardNormal

        nlgh'' :: Natural # LGH
        nlgh'' = joinConjugatedHarmonium (fst $ split nlgh) $ toNatural stnd

    --print $ roundSD 2 <$> listCoordinates (snd $ splitConjugatedHarmonium nlgh)


    --putStrLn "Randomish:"
    --comparison nlgh . toNatural $ toMean nlgh

    --putStrLn "Standard Normal Prior:"
    --comparison nlgh'' . toNatural $ toMean nlgh''
    --print $ roundSD 2 <$> listCoordinates nlgh''
    --print $ roundSD 2 <$> listCoordinates (toNatural $ toMean nlgh'')

    --print $ roundSD 2 <$> listCoordinates nlgh'
    --print $ roundSD 2 <$> listCoordinates (toMean nlgh')
    --print $ roundSD 2 <$> listCoordinates (toNatural $ toSource nlgh)
    --print $ roundSD 2 <$> listCoordinates (toNatural . toSource . toNatural $ toSource nlgh)

    --putStrLn "Isotransform0 (Source - Natural) error:"
    --print . euclideanDistance snrmx . toSource $ toNatural snrmx

    --putStrLn "Isotransform0 (Source - Mean) error:"
    --print . euclideanDistance snrmx . toSource $ toMean snrmx

    --putStrLn "Isotransform (Source - Natural) error:"
    --print . euclideanDistance slgh . toSource $ toNatural slgh
    --print slgh
    --print . toSource $ toNatural slgh

    --putStrLn "Isotransform (Source - Mean) error:"
    --print . euclideanDistance slgh . toSource $ toMean slgh

    putStrLn "Isotransform (Natural - Mean) error:"
    print . euclideanDistance nlgh . toNatural $ toMean nlgh
    print . euclideanDistance hmog . toNatural $ toMean hmog

    xyzs <- realize $ sample nsmps hmog
    --let (xs,zs) = unzip xzs
    --let nprr = snd $ splitConjugatedHarmonium nlgh
    ----zs' <- realize $ sample nsmps nprr

    let hmog' = averageSufficientStatistic xyzs
    let nmog' = expectationMaximization (fst <$> xyzs) hmog

    --let (maff,mprr) = split mlgh
    --    (maff',mprr') = split mlgh'

    putStrLn "HMoG Sampling error:"
    print $ euclideanDistance (toMean hmog) hmog'

    putStrLn "New Isotransform (Natural - Mean) error:"
    print . euclideanDistance hmog $ toNatural hmog'
    print $ euclideanDistance hmog nmog'

    ----putStrLn "Natural Sampling error:"
    ----print $ euclideanDistance nlgh $ transition mlgh'
    ----putStrLn "Source Sampling error:"
    ----print $ euclideanDistance slgh $ transition mlgh'

    ----print . euclideanDistance mprr' $ averageSufficientStatistic zs'
    --let cvr0 :: S.Matrix 5 5 Double
    --    cvr0 =
    --        let top = S.horizontalConcat (toMatrix $ toTensor strngx) (toMatrix stnsxz)
    --            btm = S.horizontalConcat (S.transpose $ toMatrix stnsxz) (toMatrix $ toTensor strngz)
    --         in S.verticalConcat top btm

    --let cvr = fromTensor $ fromMatrix cvr0
    --    smvn :: Source # FullNormal 5
    --    smvn = join (fromTuple (1,1,1,-1,-1)) cvr

    --xzs0' <- realize $ sample nsmps smvn
    --let xs = fst <$> xzs
    --    xzs' = S.splitAt <$> xzs0'
    --    (xs',zs') = unzip xzs'
    --    xmrg0,xmrg',xmrg'' :: Mean # FullNormal 3
    --    xmrg0 = toMean snrmx
    --    xmrg' = averageSufficientStatistic xs
    --    xmrg'' = averageSufficientStatistic xs'
    --    xmrg = snd . split $ transposeHarmonium mlgh

    ----print $ euclideanDistance xmrg0 xmrg
    ----print . euclideanDistance (averageSufficientStatistic zs) . transition . snd $ splitConjugatedHarmonium nlgh
    ----print . euclideanDistance (averageSufficientStatistic zs') . transition . snd $ splitConjugatedHarmonium nlgh
    --print . euclideanDistance (toMean smvn) $ averageSufficientStatistic xzs0'
    --print . euclideanDistance (toMean smvn) . averageSufficientStatistic $ zipWith (S.++) xs zs
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


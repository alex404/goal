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

trngx :: Source # Symmetric (MVNMean 3) (MVNMean 3)
trngx = fromTuple (1,2,3,4,5,6)

snrmx :: Source # SymmetricNormal 3
snrmx = join 1 trngx

trngz :: Source # Symmetric (MVNMean 2) (MVNMean 2)
trngz = fromTuple (1,2,2)

snrmz :: Source # SymmetricNormal 2
snrmz = join (-1) trngz

stnsxz :: Source # Tensor (MVNMean 3) (MVNMean 2)
stnsxz = fromTuple (2,2,2,1,1,1)

slgh :: Source # SymmetricGaussianHarmonium 3 2
slgh = join (join snrmx stnsxz) snrmz

--diag :: Source # Diagonal (Euclidean 2) (Euclidean 2)
--diag = fromTuple (1,2)
--
--scl :: Source # Scale (Euclidean 2) (Euclidean 2)
--scl = singleton 2
--
--tns :: Source # Tensor (Euclidean 2) (Euclidean 3)
--tns = fromTuple (1,2,2,3,2,1)
--
--tns' :: Source # Tensor (Euclidean 2) (Euclidean 2)
--tns' = fromTuple (1,2,2,3)
--
--hmog ::
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
    --print . map round $ listCoordinates slgh
    --print . map round . listCoordinates . toSource $ toNatural slgh
    let (sfxz,sz) = split slgh
        (sx,svrxz) = split sfxz
        (smux,svrx) = split sx
        (smuz,svrz) = split sz
        (nvrx0,nvrxz0,nvrz0) = blockSymmetricMatrixInversion svrx svrxz svrz
        nmux = breakPoint $ nvrx0 >.> smux + nvrxz0 >.> smuz
        nmuz = breakPoint $ nvrz0 >.> smuz + transpose nvrxz0 >.> smux
        nvrx = naturalPrecisionToSymmetric . breakPoint . toTensor $ -0.5 .> nvrx0
        nvrxz = breakPoint $ -nvrxz0
        nvrz = naturalPrecisionToSymmetric . breakPoint . toTensor $ -0.5 .> nvrz0
        nx = join nmux nvrx
        nz = join nmuz nvrz
        nfxz = join nx nvrxz
        nlgh :: Natural # SymmetricGaussianHarmonium 3 2
        nlgh = join nfxz nz
        --(nfxz,nz) = split nlgh
        --(nx,nvrxz) = split nfxz
        --(nmux,nvrx) = split nx
        --(nmuz,nvrz) = split nz
        (svrx0',svrxz0',svrz0') = blockSymmetricMatrixInversion
            (naturalSymmetricToPrecision nvrx) (2 /> nvrxz)
            (fromTensor $ naturalSymmetricToPrecision nvrz)
        svrx' = breakPoint $ -0.5 .> svrx0'
        svrxz' = breakPoint $ -0.5 .> svrxz0'
        svrz' = breakPoint $ -0.5 .> svrz0'
        smux' = svrx' >.> nmux + svrxz' >.> nmuz
        smuz' = svrz' >.> nmuz + transpose svrxz' >.> nmux
        sx' = join smux' svrx'
        sz' = join smuz' svrz'
        sfxz' = join sx' svrxz'
        slgh0 :: Mean # SymmetricGaussianHarmonium 3 2
        slgh0 = join sfxz' sz'
        slgh' :: Source # SymmetricGaussianHarmonium 3 2
        slgh' = breakPoint slgh'

    print . euclideanDistance svrx' $ breakPoint svrx
    print . euclideanDistance svrz' $ breakPoint svrz
    print . euclideanDistance svrxz' $ breakPoint svrxz
    print . euclideanDistance smux' $ breakPoint smux
    print . euclideanDistance smuz' $ breakPoint smuz

    --let bigMatrix :: S.Matrix 5 5 Double
    --    bigMatrix =
    --        let top = S.horizontalConcat (toMatrix $ toTensor svrx) (toMatrix svrxz)
    --            btm = S.horizontalConcat (S.transpose $ toMatrix svrxz) (toMatrix $ toTensor svrz)
    --         in S.verticalConcat top btm


    ----print bigMatrix
    --putStrLn "Simple Inverse"
    --printBig $ S.inverse bigMatrix
    --putStrLn "Advanced Inverse"
    --printBlocks nvrx0 nvrxz0 nvrz0
    --print nvrz


{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork.Convolutional where
    {-
    ( -- * Neural Networks
      InterLayer
    , splitInterLayer
    , joinInterLayer
    ) where
    -}


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily.NeuralNetwork

import GHC.TypeLits.Extra


--- Convolutional ---

data Convolutional (nk :: Nat) (rd :: Nat) (st :: Nat) (nr :: Nat) (nc :: Nat) (nd :: Nat) om im

type KnownConvolutional nk rd st nr nc nd om im =
    (KnownNat nk, KnownNat rd, KnownNat st, KnownNat nr, KnownNat nc, KnownNat nd
    , Manifold om, Manifold im, Dimension im ~ 1, Dimension om ~ 1, 1 <= st)

instance (KnownConvolutional nk rd st nr nc nd om im) => Manifold (Convolutional nk rd st nr nc nd om im) where
    type Dimension (Convolutional nk rd st nr nc nd om im) = nk * (2*rd+1) * (2*rd+1) * nd

instance (KnownConvolutional nk rd st nr nc nd om im) => Map (Convolutional nk rd st nr nc nd om im) where
    type Domain (Convolutional nk rd st nr nc nd om im) = Replicated (nr * nc * nd) im
    type Codomain (Convolutional nk rd st nr nc nd om im) = Replicated (Div nr st * Div nc st * nk) om

instance (KnownConvolutional nk rd st nr nc nd om im)
  => Apply Mean Natural (Convolutional nk rd st nr nc nd om im) where
      (>.>) cnv img = Point . toVector $ matrixMatrixMultiply (toKernelMatrix cnv) (toWindowMatrix cnv img)

toKernelMatrix
    :: KnownConvolutional nk rd st nr nc nd om im
    => Point c (Convolutional nk rd st nr nc nd om im) x -> Matrix nk ((2*rd+1) * (2*rd+1) * nd) x
toKernelMatrix prms = Matrix $ coordinates prms

toWindowMatrix
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Point d (Replicated (nr * nc * nd) im) x
    -> Matrix ((2*rd+1) * (2*rd+1) * nd) (Div nr st * Div nc st) x
toWindowMatrix cnv inp =
    fromColumns $ uncurry (cutWindow cnv Proxy Proxy Proxy (padInput cnv inp)) <$> convSteps cnv Proxy Proxy


--- Pooling ---

data Pooling (pl :: Nat) (nr :: Nat) (nc :: Nat) (nd :: Nat) im

type KnownPooling pl nr nc nd im =
    (KnownNat pl, KnownNat nr, KnownNat nc, KnownNat nd, Manifold im, Dimension im ~ 1, 1 <= pl)

instance (KnownPooling pl nr nc nd im) => Manifold (Pooling pl nr nc nd im) where
    type Dimension (Pooling pl nr nc nd im) = 0

instance (KnownPooling pl nr nc nd im) => Map (Pooling pl nr nc nd im) where
    type Domain (Pooling pl nr nc nd im) = Replicated (nr * nc * nd) im
    type Codomain (Pooling pl nr nc nd im) = Replicated (Div nr pl * Div nc pl * nd) (MeanNormal (1/1))



--- Internal ---


-- Pooling --

depthSlices
    :: (KnownPooling pl nr nc nd im, Num x)
    => Point c (Pooling pl nr nc nd im) x
    -> Vector (nr*nc*nd) x
    -> Vector nd (Matrix nr nc x)
depthSlices _ = fmap Matrix . toColumns . fromRows . breakEveryV

{-
shrinkSlice
    :: (KnownPooling pl nr nc nd im, Num x)
    => Point c (Pooling pl nr nc nd im) x
    -> Proxy pl
    -> Matrix nr nc x
    -> Matrix (Div nr pl) (Div nc pl) x
shrinkSlice _ prxypl =
    let pl = natValInt prxypl
     in fmap Matrix . toColumns . fromRows . breakEveryV
-}

-- Convolutional --

convSteps
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Proxy st
    -> Proxy nc
    -> Vector (Div nr st * Div nc st) (Int,Int)
convSteps _ prxyst prxync =
    let st = natValInt prxyst
        nc = natValInt prxync
        divMod' idx =
            let (dv,rm) = divMod idx nc
             in (dv*st,rm)
     in generateV (\idx -> divMod' $ idx*st)

cutWindow
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Proxy rd
    -> Proxy nc
    -> Proxy nd
    -> Vector ((2*rd + nr) * (2*rd + nc) * nd) x
    -> Int
    -> Int
    -> Vector ((2*rd+1) * (2*rd+1) * nd) x
cutWindow _ prxyrd prxync prxynd v r c =
    let rd = natValInt prxyrd
        wn = 2*rd + 1
        nc = natValInt prxync
        nd = natValInt prxynd
        reIndex idx =
            let (ws,wd) = divMod idx nd
                (wr,wc) = divMod ws wn
             in (r + wr) * (nc+2*rd) * nd + ((c + wc) * nd) + wd
     in backpermuteV v $ generateV reIndex

padInput
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Point d (Replicated (nr * nc * nd) im) x
    -> Vector ((2*rd + nr) * (2*rd + nc) * nd) x
padInput cnv = flattenV . toVector . padLeft cnv . padDown cnv . padRight cnv . padUp cnv . breakInput cnv

breakInput
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Point d (Replicated (nr * nc * nd) im) x
    -> Matrix nr nc (Vector nd x)
breakInput _ p =
    Matrix . breakEveryV $ coordinates p

padUp
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Matrix nr nc (Vector nd x)
    -> Matrix (rd + nr) nc (Vector nd x)
padUp _ mtx =
    fromRows . joinV (replicateV . replicateV $ replicateV 0) $ toRows mtx

padRight
    :: (KnownConvolutional nk rd st nr nc nd om im, Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Matrix (rd + nr) nc (Vector nd x)
    -> Matrix (rd + nr) (rd + nc) (Vector nd x)
padRight _ mtx =
    fromColumns . joinV (replicateV . replicateV $ replicateV 0) $ toColumns mtx

padDown
    :: (KnownConvolutional nk rd st nr nc nd om im, KnownNat nr', Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Matrix nr' (rd + nc) (Vector nd x)
    -> Matrix (rd + nr') (rd + nc) (Vector nd x)
padDown _ mtx =
    fromRows . joinV (toRows mtx) . replicateV . replicateV $ replicateV 0

padLeft
    :: (KnownConvolutional nk rd st nr nc nd om im, KnownNat nc', Num x)
    => Point c (Convolutional nk rd st nr nc nd om im) x
    -> Matrix (2*rd + nr) nc' (Vector nd x)
    -> Matrix (2*rd + nr) (rd + nc') (Vector nd x)
padLeft _ mtx =
    fromColumns . joinV (toColumns mtx) . replicateV . replicateV $ replicateV 0

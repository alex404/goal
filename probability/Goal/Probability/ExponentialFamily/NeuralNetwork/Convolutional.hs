{-# LANGUAGE BangPatterns,UndecidableInstances #-}

-- | Multilayer perceptrons aid backpropagation.
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

import qualified Goal.Core.Vector.Storable as S
import GHC.TypeLits.Extra


--- Convolutional ---

data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) (ih :: Nat) (oh :: Nat) om im

type KnownConvolutional rd r c ih oh om im =
    (KnownNat rd, KnownNat r, KnownNat c, KnownNat ih, KnownNat oh
    , Manifold om, Manifold im, Dimension im ~ 1, Dimension om ~ 1)

instance (KnownConvolutional rd r c ih oh om im) => Manifold (Convolutional rd r c ih oh om im) where
    type Dimension (Convolutional rd r c ih oh om im) = (2*rd+1) * (2*rd+1) * ih * oh

instance (KnownConvolutional rd r c ih oh om im) => Map (Convolutional rd r c ih oh om im) where
    type Domain (Convolutional rd r c ih oh om im) = Replicated (r * c * ih) im
    type Codomain (Convolutional rd r c ih oh om im) = Replicated (r * c * oh) om

instance (KnownConvolutional rd r c ih oh om im)
  => Apply Mean Natural (Convolutional rd r c ih oh om im) where
    {-# INLINE (>.>) #-}
    (>.>) cnv img = Point . convolve Proxy cnv $ coordinates img

convolve
    :: (Num x, KnownConvolutional rd r c ih oh om im)
    => Proxy oc
    -> Point (Function Mean Natural) (Convolutional rd r c ih oh om im) x
    -> Vector (st * st * or * oc * nd) x
    -> Vector (or * oc * nk) x
{-# INLINE convolve #-}
convolve prxoc cnv img =
    let kmtx = toKernelMatrix cnv
        oc = natValInt prxoc
        generator idx =
            let (r,c) = divMod idx oc
                !foo = matrixVectorMultiply kmtx $! cutWindow cnv Proxy Proxy Proxy Proxy img r c
             in foo
     in flattenV $! generateV generator

toKernelMatrix
    :: KnownConvolutional rd r c ih oh om im
    => Point c (Convolutional rd r c ih oh om im) x -> Matrix nk ((2*rd+1) * (2*rd+1) * nd) x
{-# INLINE toKernelMatrix #-}
toKernelMatrix prms = Matrix $ coordinates prms

cutWindow
    :: (KnownConvolutional rd r c ih oh om im, Num x)
    => Point c (Convolutional rd r c ih oh om im) x
    -> Proxy rd
    -> Proxy st
    -> Proxy oc
    -> Proxy nd
    -> Vector (st * st * or * oc * nd) x
    -> Int
    -> Int
    -> Vector ((2*rd+1) * (2*rd+1) * nd) x
{-# INLINE cutWindow #-}
cutWindow _ prxyrd prxyst prxyoc prxynd img r c =
    let rd = natValInt prxyrd
        st = natValInt prxyst
        oc = natValInt prxyoc
        nd = natValInt prxynd
        ic = st * oc
        wn = 2*rd + 1
        reIndex idx =
            let (ws,wd) = divMod idx nd
                (wr,wc) = divMod ws wn
             in (r + wr) * ic * nd + (c + wc) * nd + wd
      in fromMaybe 0 . lookupV img <$> generateV reIndex

--toWindowMatrix
--    :: (KnownConvolutional rd r c ih oh om im, Num x)
--    => Point c (Convolutional rd r c ih oh om im) x
--    -> Point d (Replicated (or * oc * nd) im) x
--    -> Matrix ((2*rd+1) * (2*rd+1) * nd) (Div or st * Div oc st) x
--toWindowMatrix cnv inp =
--    fromColumns $ uocurry (cutWindow cnv Proxy Proxy Proxy (padInput cnv inp)) <$> convSteps cnv Proxy Proxy
--
--
----- Pooling ---
--
--data Pooling (pl :: Nat) (or :: Nat) (oc :: Nat) (nd :: Nat) im
--
--type KnownPooling pl or oc nd im =
--    (KnownNat pl, KnownNat or, KnownNat oc, KnownNat nd, Manifold im, Dimension im ~ 1, 1 <= pl)
--
--instance (KnownPooling pl or oc nd im) => Manifold (Pooling pl or oc nd im) where
--    type Dimension (Pooling pl or oc nd im) = 0
--
--instance (KnownPooling pl or oc nd im) => Map (Pooling pl or oc nd im) where
--    type Domain (Pooling pl or oc nd im) = Replicated (or * oc * nd) im
--    type Codomain (Pooling pl or oc nd im) = Replicated (Div or pl * Div oc pl * nd) (MeanNormal (1/1))
--
--
--
----- Internal ---
--
--
---- Pooling --
--
----depthSlices
----    :: (KnownPooling pl or oc nd im, Num x)
----    => Point c (Pooling pl or oc nd im) x
----    -> Vector (or*oc*nd) x
----    -> Vector nd (Matrix or oc x)
----depthSlices _ = fmap Matrix . toColumns . fromRows . breakEveryV
----
----shrinkSlice
----    :: (KnownPooling pl or oc id im, Num x)
----    => Point c (Pooling pl or oc id im) x
----    -> Proxy pl
----    -> Matrix or oc x
----    -> Matrix (Div or pl) (Div oc pl) x
----shrinkSlice _ prxypl =
----    let pl = natValInt prxypl
----     in fmap Matrix . toColumns . fromRows . breakEveryV
--
---- Convolutional --
--
--convSteps
--    :: (KnownConvolutional nk rd st or oc id om im, Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Proxy st
--    -> Proxy oc
--    -> Vector (Div or st * Div oc st) (Int,Int)
--convSteps _ prxyst prxyoc =
--    let st = natValInt prxyst
--        oc = natValInt prxyoc
--        divMod' idx =
--            let (dv,rm) = divMod idx oc
--             in (dv*st,rm)
--     in generateV (\idx -> divMod' $ idx*st)
--
--
--padInput
--    :: (KnownConvolutional nk rd st or oc id om im, Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Point d (Replicated (or * oc * id) im) x
--    -> Vector ((2*rd + or) * (2*rd + oc) * id) x
--padInput cnv = flattenV . toVector . padLeft cnv . padDown cnv . padRight cnv . padUp cnv . breakInput cnv
--
--breakInput
--    :: (KnownConvolutional nk rd st or oc id om im, Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Point d (Replicated (or * oc * id) im) x
--    -> Matrix or oc (Vector id x)
--breakInput _ p =
--    Matrix . breakEveryV $ coordinates p
--
--padUp
--    :: (KnownConvolutional nk rd st or oc id om im, Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Matrix or oc (Vector id x)
--    -> Matrix (rd + or) oc (Vector id x)
--padUp _ mtx =
--    fromRows . joinV (replicateV . replicateV $ replicateV 0) $ toRows mtx
--
--padRight
--    :: (KnownConvolutional nk rd st or oc id om im, Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Matrix (rd + or) oc (Vector id x)
--    -> Matrix (rd + or) (rd + oc) (Vector id x)
--padRight _ mtx =
--    fromColumns . joinV (replicateV . replicateV $ replicateV 0) $ toColumns mtx
--
--padDown
--    :: (KnownConvolutional nk rd st or oc id om im, KnownNat or', Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Matrix or' (rd + oc) (Vector id x)
--    -> Matrix (rd + or') (rd + oc) (Vector id x)
--padDown _ mtx =
--    fromRows . joinV (toRows mtx) . replicateV . replicateV $ replicateV 0
--
--padLeft
--    :: (KnownConvolutional nk rd st or oc id om im, KnownNat oc', Num x)
--    => Point c (Convolutional nk rd st or oc id om im) x
--    -> Matrix (2*rd + or) oc' (Vector id x)
--    -> Matrix (2*rd + or) (rd + oc') (Vector id x)
--padLeft _ mtx =
--    fromColumns . joinV (toColumns mtx) . replicateV . replicateV $ replicateV 0

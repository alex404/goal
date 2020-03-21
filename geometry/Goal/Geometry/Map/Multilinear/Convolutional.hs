{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE ConstraintKinds,TypeApplications,UndecidableInstances #-}

-- | Manifolds of 'Convolutional' operators.
module Goal.Geometry.Map.Multilinear.Convolutional
    ( -- * Convolutional Manifolds
      Convolutional
    , KnownConvolutional
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Linear
import Goal.Geometry.Map.Multilinear
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S


-- Convolutional Layers --

-- | A 'Manifold' of correlational/convolutional transformations, defined by the
-- number of kernels, their radius, the depth of the input, and its number of
-- rows and columns.
data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) :: Type -> Type -> Type

-- | A convenience type for ensuring that all the type-level Nats of a
-- 'Convolutional' 'Manifold's are 'KnownNat's.
type KnownConvolutional rd r c z x
  = ( KnownNat rd, KnownNat r, KnownNat c, 1 <= (r*c)
    , Dimension x ~ (Div (Dimension x) (r*c) * r*c)
    , Dimension z ~ (Div (Dimension z) (r*c) * r*c)
    , Manifold (Convolutional rd r c z x)
    , Manifold x, Manifold z
    , KnownNat (Div (Dimension x) (r*c))
    , KnownNat (Div (Dimension z) (r*c))
    )

inputToImage
    :: (KnownConvolutional rd r c z x)
    => Function a b # Convolutional rd r c z x
    -> a # x
    -> S.Matrix (Div (Dimension x) (r*c)) (r*c) Double
{-# INLINE inputToImage #-}
inputToImage _ (Point img) = G.Matrix img

outputToImage
    :: (KnownConvolutional rd r c z x)
    => Function a b # Convolutional rd r c z x
    -> Dual b # z
    -> S.Matrix (Div (Dimension z) (r*c)) (r*c) Double
{-# INLINE outputToImage #-}
outputToImage _ (Point img) = G.Matrix img

layerToKernels
    :: ( KnownConvolutional rd r c z x)
    => a # Convolutional rd r c z x
    -> S.Matrix (Div (Dimension z) (r*c)) (Div (Dimension x) (r*c) * (2*rd+1)*(2*rd+1)) Double
{-# INLINE layerToKernels #-}
layerToKernels (Point krns) = G.Matrix krns

convolveApply
    :: forall a b rd r c z x
    . KnownConvolutional rd r c z x
    => Function a b # Convolutional rd r c z x
    -> a # x
    -> b # z
{-# INLINE convolveApply #-}
convolveApply cnv imp =
    let img :: S.Matrix (Div (Dimension x) (r*c)) (r*c) Double
        img = inputToImage cnv imp
        krns :: S.Matrix (Div (Dimension z) (r*c)) (Div (Dimension x) (r*c) * (2*rd+1)*(2*rd+1)) Double
        krns = layerToKernels cnv
     in Point . G.toVector
         $ S.crossCorrelate2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolveTranspose
    :: forall a b rd r c z x
    . KnownConvolutional rd r c z x
    => Function a b # Convolutional rd r c z x
    -> Function (Dual b) (Dual a) # Convolutional rd r c x z
{-# INLINE convolveTranspose #-}
convolveTranspose cnv =
    let krns = layerToKernels cnv
        pnk = Proxy :: Proxy (Div (Dimension z) (r*c))
        pmd = Proxy :: Proxy (Div (Dimension x) (r*c))
        krn' :: S.Matrix (Div (Dimension x) (r*c)) (Div (Dimension z) (r*c)*(2*rd+1)*(2*rd+1)) Double
        krn' = S.kernelTranspose pnk pmd (Proxy @ rd) (Proxy @ rd) krns
     in Point $ G.toVector krn'

--convolveTransposeApply
--    :: forall a rd r c z x . KnownConvolutional rd r c z x
--    => Dual a # z
--    -> a #> Convolutional rd r c z x
--    -> a # x
--{-# INLINE convolveTransposeApply #-}
--convolveTransposeApply imp cnv =
--    let img = outputToImage cnv imp
--        krns = layerToKernels cnv
--     in Point . G.toVector
--         $ S.convolve2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolutionalOuterProduct
    :: forall a b rd r c z x
    . KnownConvolutional rd r c z x
      => b # z
      -> Dual a # x
      -> Function a b # Convolutional rd r c z x
{-# INLINE convolutionalOuterProduct #-}
convolutionalOuterProduct (Point oimg) (Point iimg) =
    let omtx = G.Matrix oimg
        imtx = G.Matrix iimg
     in Point . G.toVector $ S.kernelOuterProduct (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) omtx imtx

convolvePropagate
    :: forall a b rd r c z x . KnownConvolutional rd r c z x
      => [b #* z]
      -> [a # x]
      -> Function a b # Convolutional rd r c z x
      -> (Function a b #* Convolutional rd r c z x, [b # z])
{-# INLINE convolvePropagate #-}
convolvePropagate omps imps cnv =
    let prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
        foldfun (omp,imp) (k,dkrns) =
            let img = inputToImage cnv imp
                dimg = outputToImage cnv omp
                dkrns' = Point . G.toVector $ S.kernelOuterProduct prdkr prdkc pmr pmc dimg img
             in (k+1,dkrns' + dkrns)
     in (uncurry (/>) . foldr foldfun (0,0) $ zip omps imps, cnv >$> imps)


--- Instances ---


-- Convolutional Manifolds --

instance ( 1 <= (r*c), Manifold x, Manifold y, KnownNat r, KnownNat c, KnownNat rd
         , KnownNat (Div (Dimension x) (r*c)) , KnownNat (Div (Dimension y) (r*c)) )
  => Manifold (Convolutional rd r c y x) where
      type Dimension (Convolutional rd r c y x)
        = (Div (Dimension y) (r * c) * ((Div (Dimension x) (r * c) * ((2 * rd) + 1)) * ((2 * rd) + 1)))


instance KnownConvolutional rd r c z x => Map a b (Convolutional rd r c) z x where
      {-# INLINE (>.>) #-}
      (>.>) = convolveApply
      {-# INLINE (>$>) #-}
      (>$>) cnv = map (convolveApply cnv)

instance KnownConvolutional rd r c z x => Bilinear (Convolutional rd r c) z x where
    {-# INLINE (>.<) #-}
    (>.<) = convolutionalOuterProduct
    {-# INLINE (>$<) #-}
    (>$<) ps qs = sum $ zipWith convolutionalOuterProduct ps qs
    {-# INLINE transpose #-}
    transpose = convolveTranspose

instance KnownConvolutional rd r c z x => Propagate a b (Convolutional rd r c) z x where
    {-# INLINE propagate #-}
    propagate = convolvePropagate

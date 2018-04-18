{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      InterLayer
    , type (<*<)
    , type (:+:)
    , splitInterLayer
    , joinInterLayer
    -- ** Convolutional Layers
    , Convolutional
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

--- Multilayer ---


data InterLayer f g

type (f :+: g) = InterLayer f g
infixr 3 :+:

type (m <*< g) = m <* Codomain g :+: g
infixr 3 <*<

splitInterLayer :: (Manifold m, Manifold n) => Point c (InterLayer m n) -> (Point c m, Point c n)
{-# INLINE splitInterLayer #-}
splitInterLayer (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

joinInterLayer :: (Manifold m, Manifold n) => Point c m -> Point c n -> Point c (InterLayer m n)
{-# INLINE joinInterLayer #-}
joinInterLayer (Point xms) (Point xns) =
    Point $ xms S.++ xns

-- Convolutional Layers --

data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) (ih :: Nat) (oh :: Nat) om im

inputToImage
    :: KnownConvolutional rd r c ih oh om im
    => Mean ~> Natural # Convolutional rd r c ih oh om im
    -> Mean # Domain (Convolutional rd r c ih oh om im)
    -> S.Vector ih (S.Matrix r c Double)
{-# INLINE inputToImage #-}
inputToImage _ (Point img) =
    S.map G.Matrix $ S.breakEvery img

outputToImage
    :: KnownConvolutional rd r c ih oh om im
    => Mean ~> Natural # Convolutional rd r c ih oh om im
    -> Mean # Codomain (Convolutional rd r c ih oh om im)
    -> S.Vector oh (S.Matrix r c Double)
{-# INLINE outputToImage #-}
outputToImage _ (Point img) =
    S.map G.Matrix $ S.breakEvery img

layerToKernels
    :: KnownConvolutional rd r c ih oh om im
    => a # Convolutional rd r c ih oh om im
    -> S.Vector ih (S.Vector oh (S.Matrix (2*rd+1) (2*rd+1) Double))
{-# INLINE layerToKernels #-}
layerToKernels (Point img) =
    S.map (S.map G.Matrix) . S.breakEvery $ S.breakEvery img

convolvePropagate :: KnownConvolutional rd r c ih oh om im
          => Point Mean (Codomain (Convolutional rd r c ih oh om im))
          -> Point Mean (Domain (Convolutional rd r c ih oh om im))
          -> Mean ~> Natural # Convolutional rd r c ih oh om im
          -> ( Natural ~> Mean # Convolutional rd r c ih oh om im
             , Point Natural (Codomain (Convolutional rd r c ih oh om im)) )
convolvePropagate omp imp cnv =
    let img = inputToImage cnv imp
        img' = outputToImage cnv omp
     in undefined

-- Convolutional Manifolds --

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
    (>.>) cnv imp = let img = inputToImage cnv imp
                        krns = layerToKernels cnv
                     in Point . S.concatMap G.toVector $ S.crossCorrelate2d krns img

instance (KnownConvolutional rd r c ih oh om im)
  => Bilinear Mean Natural (Convolutional rd r c ih oh om im) where
    {-# INLINE (<.<) #-}
    (<.<) imp' cnv = let img' = outputToImage cnv imp'
                         krns = layerToKernels cnv
                      in Point . S.concatMap G.toVector $ S.convolve2d krns img'

--instance (KnownConvolutional rd r c ih oh om im) => Propagate Mean Natural (Convolutional rd r c ih oh om im) where
--    {-# INLINE propagate #-}
--    propagate dps qs pq =
--        let dpss = splitReplicated dps
--            qss = splitReplicated qs
--            foldfun dmtx dps' qs' = (dps' >.< qs') <+> dmtx
--            n = S.length dpss
--         in (fromIntegral n /> S.zipFold foldfun zero dpss qss, pq >$> qs)

-- General Neural Networks --

instance (Manifold f, Manifold g) => Manifold (InterLayer f g) where
    type Dimension (InterLayer f g) = Dimension f + Dimension g

instance (Map f, Map g, Codomain g ~ Domain f) => Map (InterLayer f g) where
    type Domain (InterLayer f g) = Domain g
    type Codomain (InterLayer f g) = Codomain f

instance (d ~ Dual c, Apply c d f, Apply c d g, Transition d c (Codomain g), Codomain g ~ Domain f)
  => Apply c d (InterLayer f g) where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitInterLayer fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitInterLayer fg
         in f >$> mapReplicatedPoint transition (g >$> xs)

instance {-# OVERLAPPABLE #-} (Domain f ~ Codomain g, Manifold g, Propagate Mean Natural g, Legendre Natural (Codomain g), Riemannian Natural (Domain f), Bilinear Mean Natural f, Propagate Mean Natural f)
  => Propagate Mean Natural (InterLayer (Affine f) g) where
      propagate dps qs fg =
          let (f,g) = splitInterLayer fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              dhs = dualIsomorphism . detachTangentVector . flat
                  . joinTangentPair hs . Point . coordinates $ dps <$< fmtx
           in (joinInterLayer df dg, phts)


--instance {-# OVERLAPPING #-} (Replicated k Bernoulli ~ Codomain g, Manifold g, Propagate Mean Natural g, Legendre Natural (Codomain g), Manifold m, KnownNat k)
--  => Propagate Mean Natural (InterLayer (Affine (Product m (Replicated k Bernoulli))) g) where
--      propagate dps qs fg =
--          let (f,g) = splitInterLayer fg
--              fmtx = snd $ splitAffine f
--              mhs = mapReplicatedPoint dualTransition hs
--              (df,phts) = propagate dps mhs f
--              (dg,hs) = propagate dhs qs g
--              thts = S.map (\x -> x * (1-x)) $ coordinates mhs
--              dhs = Point . S.zipWith (*) thts . coordinates $ dps <$< fmtx
--           in (joinInterLayer df dg, phts)

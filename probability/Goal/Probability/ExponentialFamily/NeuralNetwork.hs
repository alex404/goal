{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      NeuralNetworkLayer
    , HiddenNeuralNetwork
    , NeuralNetwork
    --, type (<*<)
    , splitNeuralNetworkLayer
    , joinNeuralNetworkLayer
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S

--- Multilayer ---


data NeuralNetworkLayer (f :: * -> * -> *) (g :: * -> * -> *) n m o

type family NeuralNetwork' (fs :: [* -> * -> *]) (hs :: [*]) where
    NeuralNetwork' '[f] '[] = f
    NeuralNetwork' (f : fs) (h : hs) = NeuralNetworkLayer f (NeuralNetwork' fs hs) h

type family NeuralNetwork (fs :: [* -> * -> *]) (hs :: [*]) where
    NeuralNetwork '[f] '[] = f
    NeuralNetwork (f : fs) (h : hs) = NeuralNetworkLayer f (NeuralNetwork' fs hs) h

--type (f :+: g) = NeuralNetworkLayer f g
--infixr 3 :+:

--type (m <*< g) n o = NeuralNetworkLayer (Affine Product) g n m o
--infixr 3 <*<

splitNeuralNetworkLayer
    :: (Manifold (f m n), Manifold (g n o))
    => c # NeuralNetworkLayer f g n m o
    -> (c # f m n, c # g n o)
{-# INLINE splitNeuralNetworkLayer #-}
splitNeuralNetworkLayer (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

joinNeuralNetworkLayer
    :: (Manifold m, Manifold n)
    => c # f m n
    -> c # g n o
    -> c # NeuralNetworkLayer f g n m o
{-# INLINE joinNeuralNetworkLayer #-}
joinNeuralNetworkLayer (Point xms) (Point xns) =
    Point $ xms S.++ xns

---- Convolutional Layers --
--
--data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) (ih :: Nat) (oh :: Nat) om im
--
--inputToImage
--    :: KnownConvolutional rd r c ih oh om im
--    => Mean ~> Natural # Convolutional rd r c ih oh om im
--    -> Mean # Domain (Convolutional rd r c ih oh om im)
--    -> S.Matrix ih (r*c) Double
--{-# INLINE inputToImage #-}
--inputToImage _ (Point img) = G.Matrix img
--
--outputToImage
--    :: KnownConvolutional rd r c ih oh om im
--    => Mean ~> Natural # Convolutional rd r c ih oh om im
--    -> Mean # Codomain (Convolutional rd r c ih oh om im)
--    -> S.Matrix oh (r*c) Double
--{-# INLINE outputToImage #-}
--outputToImage _ (Point img) = G.Matrix img
--
--layerToKernels
--    :: KnownConvolutional rd r c ih oh om im
--    => a # Convolutional rd r c ih oh om im
--    -> S.Matrix oh (ih*(2*rd+1)*(2*rd+1)) Double
--{-# INLINE layerToKernels #-}
--layerToKernels (Point krns) = G.Matrix krns
--
--convolveApply
--    :: forall rd r c ih oh om im
--    . KnownConvolutional rd r c ih oh om im
--    => Mean ~> Natural # Convolutional rd r c ih oh om im
--    -> Point Mean (Domain (Convolutional rd r c ih oh om im))
--    -> Point Natural (Codomain (Convolutional rd r c ih oh om im))
--{-# INLINE convolveApply #-}
--convolveApply cnv imp =
--    let img = inputToImage cnv imp
--        krns = layerToKernels cnv
--        prdkr = Proxy :: Proxy rd
--        prdkc = Proxy :: Proxy rd
--        pmr = Proxy :: Proxy r
--        pmc = Proxy :: Proxy c
--     in Point . G.toVector $ S.crossCorrelate2d prdkr prdkc pmr pmc krns img
--
--convolveTransposeApply
--    :: forall rd r c ih oh om im
--    . KnownConvolutional rd r c ih oh om im
--    => Point Mean (Codomain (Convolutional rd r c ih oh om im))
--    -> Mean ~> Natural # Convolutional rd r c ih oh om im
--    -> Point Natural (Domain (Convolutional rd r c ih oh om im))
--{-# INLINE convolveTransposeApply #-}
--convolveTransposeApply imp cnv =
--    let img = outputToImage cnv imp
--        krns = layerToKernels cnv
--        prdkr = Proxy :: Proxy rd
--        prdkc = Proxy :: Proxy rd
--        pmr = Proxy :: Proxy r
--        pmc = Proxy :: Proxy c
--     in Point . G.toVector $ S.convolve2d prdkr prdkc pmr pmc krns img
--
--convolvePropagate
--    :: forall rd r c ih oh om im k
--    . (KnownConvolutional rd r c ih oh om im, KnownNat k)
--    => Point Mean (Replicated k (Codomain (Convolutional rd r c ih oh om im)))
--    -> Point Mean (Replicated k (Domain (Convolutional rd r c ih oh om im)))
--    -> Mean ~> Natural # Convolutional rd r c ih oh om im
--    -> ( Natural ~> Mean # Convolutional rd r c ih oh om im
--       , Point Natural (Replicated k (Codomain (Convolutional rd r c ih oh om im))) )
--{-# INLINE convolvePropagate #-}
--convolvePropagate omps0 imps0 cnv =
--    let prdkr = Proxy :: Proxy rd
--        prdkc = Proxy :: Proxy rd
--        pmr = Proxy :: Proxy r
--        pmc = Proxy :: Proxy c
--        omps = splitReplicated omps0
--        imps = splitReplicated imps0
--        foldfun dkrns omp imp =
--            let img = inputToImage cnv imp
--                dimg = outputToImage cnv omp
--                dkrns' = Point . G.toVector $ S.kernelDifferential prdkr prdkc pmr pmc dimg img
--             in dkrns' <+> dkrns
--        n = S.length omps
--     in (fromIntegral n /> S.zipFold foldfun zero omps imps, cnv >$> imps0)
--
---- Convolutional Manifolds --
--
--type KnownConvolutional rd r c ih oh om im =
--    (KnownNat rd, KnownNat r, KnownNat c, KnownNat ih, KnownNat oh
--    , Manifold om, Manifold im, Dimension im ~ 1, Dimension om ~ 1)
--
--instance (KnownConvolutional rd r c ih oh om im) => Manifold (Convolutional rd r c ih oh om im) where
--    type Dimension (Convolutional rd r c ih oh om im) = (2*rd+1) * (2*rd+1) * ih * oh
--
--instance (KnownConvolutional rd r c ih oh om im) => Map (Convolutional rd r c ih oh om im) where
--    type Domain (Convolutional rd r c ih oh om im) = Replicated (r * c * ih) im
--    type Codomain (Convolutional rd r c ih oh om im) = Replicated (r * c * oh) om
--
--instance (KnownConvolutional rd r c ih oh om im)
--  => Apply Mean Natural (Convolutional rd r c ih oh om im) where
--    {-# INLINE (>.>) #-}
--    (>.>) = convolveApply
--
--instance (KnownConvolutional rd r c ih oh om im)
--  => Bilinear Mean Natural (Convolutional rd r c ih oh om im) where
--    {-# INLINE (<.<) #-}
--    (<.<) = convolveTransposeApply
--
--instance (KnownConvolutional rd r c ih oh om im) => Propagate Mean Natural (Convolutional rd r c ih oh om im) where
--    {-# INLINE propagate #-}
--    propagate = convolvePropagate
--
-- General Neural Networks --

instance (Manifold (f m n), Manifold (g n o)) => Manifold (NeuralNetworkLayer f g n m o) where
      type Dimension (NeuralNetworkLayer f g n m o) = Dimension (f m n) + Dimension (g n o)

instance (Map c d f m n, Map c d g n o, Transition d c n)
  => Map c d (NeuralNetworkLayer f g n) m o where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitNeuralNetworkLayer fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitNeuralNetworkLayer fg
         in f >$> mapReplicatedPoint transition (g >$> xs)

instance {-# OVERLAPPABLE #-}
    ( Propagate Mean Natural (Affine f) m n, Propagate Mean Natural g n o
    , Legendre Natural n, Riemannian Natural n, Bilinear f m n)
  => Propagate Mean Natural (NeuralNetworkLayer (Affine f) g n) m o where
      propagate dps qs fg =
          let (f,g) = splitNeuralNetworkLayer fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              dhs = dualIsomorphism . detachTangentVector . flat
                  . joinTangentPair hs . Point . coordinates $ dps <$< fmtx
           in (joinNeuralNetworkLayer df dg, phts)

--instance {-# OVERLAPPING #-}
--    ( KnownNat k, n ~ Replicated k Bernoulli, Propagate Mean Natural (Affine f) m n
--    , Propagate Mean Natural g n o, Bilinear f m n)
--    => Propagate Mean Natural (NeuralNetworkLayer (Affine f) g (Replicated k Bernoulli)) m o where
--      propagate dps qs fg =
--          let (f,g) = splitNeuralNetworkLayer fg
--              fmtx = snd $ splitAffine f
--              mhs = mapReplicatedPoint dualTransition hs
--              (df,phts) = propagate dps mhs f
--              (dg,hs) = propagate dhs qs g
--              thts = S.map (\x -> x * (1-x)) $ coordinates mhs
--              dhs = Point . S.zipWith (*) thts . coordinates $ dps <$< fmtx
--           in (joinNeuralNetworkLayer df dg, phts)

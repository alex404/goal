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

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S

--- Multilayer ---


data NeuralNetworkLayer (f :: * -> * -> *) (g :: * -> * -> *) n m o

type family HiddenNeuralNetwork (fs :: [* -> * -> *]) (hs :: [*]) where
    HiddenNeuralNetwork '[f] '[] = Affine f
    HiddenNeuralNetwork (f : fs) (h : hs) = NeuralNetworkLayer (Affine f) (HiddenNeuralNetwork fs hs) h

type NeuralNetwork fs ms = HiddenNeuralNetwork fs (Init (Tail ms)) (Head ms) (Last ms)

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

-- Convolutional Layers --

data Convolutional (nk :: Nat) (rd :: Nat) (d :: Nat) (r :: Nat) (c :: Nat) :: * -> * -> *

type KnownConvolutional nk rd d r c
  = (KnownNat nk, KnownNat rd, KnownNat d, KnownNat r, KnownNat c, 1 <= (r*c))

inputToImage
    :: (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d))
    => Mean ~> Natural # Convolutional nk rd d r c m n
    -> Mean # n
    -> S.Matrix d (r*c) Double
{-# INLINE inputToImage #-}
inputToImage _ (Point img) = G.Matrix img

outputToImage
    :: (KnownConvolutional nk rd d r c, Dimension m ~ (nk*r*c))
    => Mean ~> Natural # Convolutional nk rd d r c m n
    -> Mean # m
    -> S.Matrix nk (r*c) Double
{-# INLINE outputToImage #-}
outputToImage _ (Point img) = G.Matrix img

layerToKernels
    :: KnownConvolutional nk rd d r c
    => a # Convolutional nk rd d r c m n
    -> S.Matrix nk (d*(2*rd+1)*(2*rd+1)) Double
{-# INLINE layerToKernels #-}
layerToKernels (Point krns) = G.Matrix krns

convolveApply
    :: forall nk rd d r c m n
    . (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r))
    => Mean ~> Natural # Convolutional nk rd d r c m n
    -> Mean # n
    -> Natural # m
{-# INLINE convolveApply #-}
convolveApply cnv imp =
    let img = inputToImage cnv imp
        krns = layerToKernels cnv
        prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
     in Point . G.toVector $ S.crossCorrelate2d prdkr prdkc pmr pmc krns img

convolveTransposeApply
    :: forall nk rd d r c m n
    . (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r))
    => Mean # m
    -> Mean ~> Natural # Convolutional nk rd d r c m n
    -> Natural # n
{-# INLINE convolveTransposeApply #-}
convolveTransposeApply imp cnv =
    let img = outputToImage cnv imp
        krns = layerToKernels cnv
        prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
     in Point . G.toVector $ S.convolve2d prdkr prdkc pmr pmc krns img

convolvePropagate
    :: forall nk rd d r c m n k
    . (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r), KnownNat k)
    => Mean # Replicated k m
    -> Mean # Replicated k n
    -> Mean ~> Natural # Convolutional nk rd d r c m n
    -> ( Natural ~> Mean # Convolutional nk rd d r c m n, Natural # Replicated k m )
{-# INLINE convolvePropagate #-}
convolvePropagate omps0 imps0 cnv =
    let prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
        omps = splitReplicated omps0
        imps = splitReplicated imps0
        foldfun dkrns omp imp =
            let img = inputToImage cnv imp
                dimg = outputToImage cnv omp
                dkrns' = Point . G.toVector $ S.kernelOuterProduct prdkr prdkc pmr pmc dimg img
             in dkrns' <+> dkrns
        n = S.length omps
     in (fromIntegral n /> S.zipFold foldfun zero omps imps, cnv >$> imps0)


--- Instances ---


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
    ( Propagate Mean Natural f m n, Propagate Mean Natural g n o
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

instance {-# OVERLAPPING #-}
    ( KnownNat k, n ~ Replicated k Bernoulli, Propagate Mean Natural (Affine f) m n
    , Propagate Mean Natural g n o, Bilinear f m n)
    => Propagate Mean Natural (NeuralNetworkLayer (Affine f) g (Replicated k Bernoulli)) m o where
      propagate dps qs fg =
          let (f,g) = splitNeuralNetworkLayer fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              thts = S.map (\x -> x * (1-x)) $ coordinates mhs
              dhs = Point . S.zipWith (*) thts . coordinates $ dps <$< fmtx
           in (joinNeuralNetworkLayer df dg, phts)


-- Convolutional Manifolds --

instance ( KnownConvolutional nk rd d r c, Manifold m, Manifold n
         , Dimension n ~ (rd*r*c), Dimension m ~ (nk*r*c) )
  => Manifold (Convolutional nk rd d r c m n) where
      type Dimension (Convolutional nk rd d r c m n)
        = (2*rd+1) * (2*rd+1) * nk * d

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



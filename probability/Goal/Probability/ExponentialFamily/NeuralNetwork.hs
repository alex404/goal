{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      HiddenNeuralNetwork
    , NeuralNetwork
    , splitHiddenNeuralNetwork
    , joinHiddenNeuralNetwork
    , Convolutional
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


data HiddenNeuralNetwork (fs :: [* -> * -> *]) (hs :: [*]) m n

type NeuralNetwork fs ms = HiddenNeuralNetwork fs (Init (Tail ms)) (Head ms) (Last ms)


fromSingleLayerNetwork :: c # HiddenNeuralNetwork '[f] '[] m n -> c # f m n
{-# INLINE fromSingleLayerNetwork #-}
fromSingleLayerNetwork = Point . coordinates

toSingleLayerNetwork :: c # f m n -> c # HiddenNeuralNetwork '[f] '[] m n
{-# INLINE toSingleLayerNetwork #-}
toSingleLayerNetwork = Point . coordinates

splitHiddenNeuralNetwork
    :: (Manifold (f m h), Manifold (HiddenNeuralNetwork fs hs h n))
    => c # HiddenNeuralNetwork (f : fs) (h : hs) m n
    -> (c # f m h, c # HiddenNeuralNetwork fs hs h n)
{-# INLINE splitHiddenNeuralNetwork #-}
splitHiddenNeuralNetwork (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

joinHiddenNeuralNetwork
    :: (Manifold (f m h), Manifold (HiddenNeuralNetwork fs hs h n))
    => c # f m h
    -> c # HiddenNeuralNetwork fs hs h n
    -> c # HiddenNeuralNetwork (f : fs) (h : hs) m n
{-# INLINE joinHiddenNeuralNetwork #-}
joinHiddenNeuralNetwork (Point xms) (Point xns) =
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
    => chrt1 ~> chrt2 # Convolutional nk rd d r c m n
    -> Dual chrt2 # m
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
    :: forall chrt1 chrt2 nk rd d r c m n
    . (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r))
    => Dual chrt2 # m
    -> chrt1 ~> chrt2 # Convolutional nk rd d r c m n
    -> Dual chrt1 # n
{-# INLINE convolveTransposeApply #-}
convolveTransposeApply imp cnv =
    let img = outputToImage cnv imp
        krns = layerToKernels cnv
        prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
     in Point . G.toVector $ S.convolve2d prdkr prdkc pmr pmc krns img

convolutionalOuterProduct
    :: forall chrt1 chrt2 nk rd d r c m n
    . (KnownConvolutional nk rd d r c, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r))
      => Point chrt2 m
      -> Point chrt1 n
      -> Dual chrt1 ~> chrt2 # Convolutional nk rd d r c m n
{-# INLINE convolutionalOuterProduct #-}
convolutionalOuterProduct (Point oimg) (Point iimg) =
    let prdkr = Proxy :: Proxy rd
        prdkc = Proxy :: Proxy rd
        pmr = Proxy :: Proxy r
        pmc = Proxy :: Proxy c
        omtx :: S.Matrix nk (r*c) Double
        omtx = G.Matrix oimg
        imtx :: S.Matrix d (r*c) Double
        imtx = G.Matrix iimg
     in Point . G.toVector $ S.kernelOuterProduct prdkr prdkc pmr pmc omtx imtx

convolvePropagate
    :: forall nk rd d r c m n k
    . ( KnownConvolutional nk rd d r c, Manifold m, Manifold n, KnownNat k
      , Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r), Manifold (Convolutional nk rd d r c m n) )
      => Mean # Replicated k m
      -> Mean # Replicated k n
      -> Mean ~> Natural # Convolutional nk rd d r c m n
      -> ( Natural ~> Mean # Convolutional nk rd d r c m n, Point Natural (Replicated k m) )
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

instance Manifold (f m n) => Manifold (HiddenNeuralNetwork '[f] '[] m n) where
      type Dimension (HiddenNeuralNetwork '[f] '[] m n) = Dimension (f m n)

instance (Manifold (f m h), Manifold (HiddenNeuralNetwork fs hs h n))
  => Manifold (HiddenNeuralNetwork (f : fs) (h : hs) m n) where
      type Dimension (HiddenNeuralNetwork (f : fs) (h : hs) m n)
        = Dimension (f m h) + Dimension (HiddenNeuralNetwork fs hs h n)

instance (Map c d f m n) => Map c d (HiddenNeuralNetwork '[f] '[]) m n where
    {-# INLINE (>.>) #-}
    (>.>) f x = fromSingleLayerNetwork f >.> x
    {-# INLINE (>$>) #-}
    (>$>) f xs = fromSingleLayerNetwork f >$> xs

instance (Map c d f m h, Map c d (HiddenNeuralNetwork fs hs) h n, Transition d c h)
  => Map c d (HiddenNeuralNetwork (f : fs) (h : hs)) m n where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitHiddenNeuralNetwork fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitHiddenNeuralNetwork fg
         in f >$> mapReplicatedPoint transition (g >$> xs)

instance (Propagate Mean Natural f m n) => Propagate Mean Natural (HiddenNeuralNetwork '[f] '[]) m n where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df,ps) = propagate dps qs $ fromSingleLayerNetwork f
         in (toSingleLayerNetwork df,ps)

instance {-# OVERLAPPABLE #-}
    ( Propagate Mean Natural f m h, Propagate Mean Natural (HiddenNeuralNetwork fs hs) h n
    , Legendre Natural h, Riemannian Natural h, Bilinear f m h)
  => Propagate Mean Natural (HiddenNeuralNetwork (Affine f : fs) (h : hs)) m n where
      propagate dps qs fg =
          let (f,g) = splitHiddenNeuralNetwork fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              dhs = dualIsomorphism . detachTangentVector . flat
                  . joinTangentPair hs . Point . coordinates $ dps <$< fmtx
           in (joinHiddenNeuralNetwork df dg, phts)

instance {-# OVERLAPPING #-}
    ( KnownNat k, h ~ Replicated k Bernoulli, Propagate Mean Natural f m h
    , Propagate Mean Natural (HiddenNeuralNetwork fs hs) h n, Bilinear f m h)
  => Propagate Mean Natural (HiddenNeuralNetwork (Affine f : fs) (Replicated k Bernoulli : hs)) m n where
      propagate dps qs fg =
          let (f,g) = splitHiddenNeuralNetwork fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              thts = S.map (\x -> x * (1-x)) $ coordinates mhs
              dhs = Point . S.zipWith (*) thts . coordinates $ dps <$< fmtx
           in (joinHiddenNeuralNetwork df dg, phts)

-- Convolutional Manifolds --

instance ( KnownConvolutional nk rd d r c, Manifold m, Manifold n
         , Dimension n ~ (rd*r*c), Dimension m ~ (nk*r*c) )
  => Manifold (Convolutional nk rd d r c m n) where
      type Dimension (Convolutional nk rd d r c m n)
        = (2*rd+1) * (2*rd+1) * nk * d

instance ( KnownConvolutional nk rd d r c, Manifold (Convolutional nk rd d r c m n)
         , Manifold m, Manifold n, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r) )
  => Map Mean Natural (Convolutional nk rd d r c) m n where
      {-# INLINE (>.>) #-}
      (>.>) = convolveApply

instance ( KnownConvolutional nk rd d r c, Manifold (Convolutional nk rd d r c m n)
         , Manifold m, Manifold n, Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r) )
  => Bilinear (Convolutional nk rd d r c) m n where
    {-# INLINE (<.<) #-}
    (<.<) = convolveTransposeApply
    {-# INLINE (>.<) #-}
    (>.<) = convolutionalOuterProduct

instance ( KnownConvolutional nk rd d r c, Manifold m, Manifold n
         , Dimension n ~ (r*c*d), Dimension m ~ (nk*c*r), Manifold (Convolutional nk rd d r c m n) )
  => Propagate Mean Natural (Convolutional nk rd d r c) m n where
    {-# INLINE propagate #-}
    propagate = convolvePropagate




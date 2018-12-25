{-# LANGUAGE
    RankNTypes,
    KindSignatures,
    DataKinds,
    TypeOperators,
    ConstraintKinds,
    TypeFamilies,
    FlexibleContexts,
    FlexibleInstances,
    MultiParamTypeClasses,
    ScopedTypeVariables,
    TypeApplications,
    UndecidableInstances
    #-}

-- | Multilayer perceptrons and backpropagation. The core type is the
-- 'NeuralNetwork', which is defined by two type-lists. The first type list is a
-- collection of maps. The second type-list is a list of 'Manifold's, which
-- defines the size and activation function of each layer of the network.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      NeuralNetwork
    , HiddenNeuralNetwork
    , splitNeuralNetwork
    , joinNeuralNetwork
    -- ** Convolutional Layers
    , Convolutional
    , KnownConvolutional
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Generic as G
import qualified Goal.Core.Vector.Storable as S


--- Multilayer ---


-- | A 'NeuralNetwork' where the input and output manifolds are explicit arguments.
data HiddenNeuralNetwork (fs :: [Type -> Type -> Type]) (hs :: [Type]) m n

-- | A 'NeuralNetwork' defined by a type-list of transformations and a type-list of layers.
type NeuralNetwork fs ms = HiddenNeuralNetwork fs (Init (Tail ms)) (Head ms) (Last ms)

fromSingleLayerNetwork :: c # HiddenNeuralNetwork '[f] '[] m n -> c # f m n
{-# INLINE fromSingleLayerNetwork #-}
fromSingleLayerNetwork = breakPoint

toSingleLayerNetwork :: c # f m n -> c # HiddenNeuralNetwork '[f] '[] m n
{-# INLINE toSingleLayerNetwork #-}
toSingleLayerNetwork = breakPoint

-- | Seperates a 'NeuralNetwork' into the final layer and the rest of the network.
splitNeuralNetwork
    :: (Manifold (f m h), Manifold (HiddenNeuralNetwork fs hs h n))
    => c # HiddenNeuralNetwork (f : fs) (h : hs) m n
    -> (c # f m h, c # HiddenNeuralNetwork fs hs h n)
{-# INLINE splitNeuralNetwork #-}
splitNeuralNetwork (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

-- | Joins a layer onto a 'NeuralNetwork' into a deeper network.
joinNeuralNetwork
    :: (Manifold (f m h), Manifold (HiddenNeuralNetwork fs hs h n))
    => c # f m h
    -> c # HiddenNeuralNetwork fs hs h n
    -> c # HiddenNeuralNetwork (f : fs) (h : hs) m n
{-# INLINE joinNeuralNetwork #-}
joinNeuralNetwork (Point xms) (Point xns) =
    Point $ xms S.++ xns

-- Convolutional Layers --

-- | A 'Manifold' of correlational/convolutional transformations, defined by the
-- number of kernels, their radius, the depth of the input, and its number of
-- rows and columns.
data Convolutional (rd :: Nat) (r :: Nat) (c :: Nat) :: Type -> Type -> Type

-- | A convenience type for ensuring that all the type-level Nats of a
-- 'Convolutional' 'Manifold's are 'KnownNat's.
type KnownConvolutional rd r c m n
  = ( KnownNat rd, KnownNat r, KnownNat c, 1 <= (r*c)
    , Dimension m ~ (Div (Dimension m) (r*c) * r*c)
    , Dimension n ~ (Div (Dimension n) (r*c) * r*c)
    , Manifold (Convolutional rd r c m n)
    , Manifold m, Manifold n
    , KnownNat (Div (Dimension n) (r*c))
    , KnownNat (Div (Dimension m) (r*c))
    )

inputToImage
    :: (KnownConvolutional rd r c m n)
    => Mean #> Natural # Convolutional rd r c m n
    -> Mean # n
    -> S.Matrix (Div (Dimension n) (r*c)) (r*c) Double
{-# INLINE inputToImage #-}
inputToImage _ (Point img) = G.Matrix img

outputToImage
    :: (KnownConvolutional rd r c m n)
    => chrt1 #> chrt2 # Convolutional rd r c m n
    -> Dual chrt2 # m
    -> S.Matrix (Div (Dimension m) (r*c)) (r*c) Double
{-# INLINE outputToImage #-}
outputToImage _ (Point img) = G.Matrix img

layerToKernels
    :: ( KnownConvolutional rd r c m n)
    => a # Convolutional rd r c m n
    -> S.Matrix (Div (Dimension m) (r*c)) (Div (Dimension n) (r*c) * (2*rd+1)*(2*rd+1)) Double
{-# INLINE layerToKernels #-}
layerToKernels (Point krns) = G.Matrix krns

convolveApply
    :: forall rd r c m n
    . KnownConvolutional rd r c m n
    => Mean #> Natural # Convolutional rd r c m n
    -> Mean # n
    -> Natural # m
{-# INLINE convolveApply #-}
convolveApply cnv imp =
    let img :: S.Matrix (Div (Dimension n) (r*c)) (r*c) Double
        img = inputToImage cnv imp
        krns :: S.Matrix (Div (Dimension m) (r*c)) (Div (Dimension n) (r*c) * (2*rd+1)*(2*rd+1)) Double
        krns = layerToKernels cnv
     in Point . G.toVector
         $ S.crossCorrelate2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolveTranspose
    :: forall chrt1 chrt2 rd r c m n
    . KnownConvolutional rd r c m n
    => chrt1 #> chrt2 # Convolutional rd r c m n
    -> Dual chrt2 #> Dual chrt1 # Convolutional rd r c n m
{-# INLINE convolveTranspose #-}
convolveTranspose cnv =
    let krns = layerToKernels cnv
        pnk = Proxy :: Proxy (Div (Dimension m) (r*c))
        pmd = Proxy :: Proxy (Div (Dimension n) (r*c))
        krn' :: S.Matrix (Div (Dimension n) (r*c)) (Div (Dimension m) (r*c)*(2*rd+1)*(2*rd+1)) Double
        krn' = S.kernelTranspose pnk pmd (Proxy @ rd) (Proxy @ rd) krns
     in Point $ G.toVector krn'

convolveTransposeApply
    :: forall chrt1 chrt2 rd r c m n
    . KnownConvolutional rd r c m n
    => Dual chrt2 # m
    -> chrt1 #> chrt2 # Convolutional rd r c m n
    -> Dual chrt1 # n
{-# INLINE convolveTransposeApply #-}
convolveTransposeApply imp cnv =
    let img = outputToImage cnv imp
        krns = layerToKernels cnv
     in Point . G.toVector
         $ S.convolve2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolutionalOuterProduct
    :: forall chrt1 chrt2 rd r c m n
    . KnownConvolutional rd r c m n
      => Point chrt2 m
      -> Point chrt1 n
      -> Dual chrt1 #> chrt2 # Convolutional rd r c m n
{-# INLINE convolutionalOuterProduct #-}
convolutionalOuterProduct (Point oimg) (Point iimg) =
    let omtx :: S.Matrix (Div (Dimension m) (r*c)) (r*c) Double
        omtx = G.Matrix oimg
        imtx :: S.Matrix (Div (Dimension n) (r*c)) (r*c) Double
        imtx = G.Matrix iimg
     in Point . G.toVector $ S.kernelOuterProduct (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) omtx imtx

convolvePropagate
    :: forall rd r c m n . KnownConvolutional rd r c m n
      => [Mean # m]
      -> [Mean # n]
      -> Mean #> Natural # Convolutional rd r c m n
      -> (Natural #> Mean # Convolutional rd r c m n, [Natural # m])
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
             in (k+1,dkrns' <+> dkrns)
     in (uncurry (/>) . foldr foldfun (0,zero) $ zip omps imps, cnv >$> imps)


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
        let (f,g) = splitNeuralNetwork fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitNeuralNetwork fg
         in f >$> map transition (g >$> xs)

instance (Propagate Mean Natural f m n) => Propagate Mean Natural (HiddenNeuralNetwork '[f] '[]) m n where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df,ps) = propagate dps qs $ fromSingleLayerNetwork f
         in (toSingleLayerNetwork df,ps)

instance {-# OVERLAPPABLE #-}
    ( Propagate Mean Natural f m h, Propagate Mean Natural (HiddenNeuralNetwork fs hs) h n
    , Transition Natural Mean h, Legendre Natural h, Riemannian Natural h, Bilinear f m h)
  => Propagate Mean Natural (HiddenNeuralNetwork (Affine f : fs) (h : hs)) m n where
      propagate dps qs fg =
          let (f,g) = splitNeuralNetwork fg
              fmtx = snd $ splitAffine f
              mhs = dualTransition <$> hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              dhs0 = dps <$< fmtx
              dhs = do
                  (h,dh0) <- zip hs dhs0
                  return . dualIsomorphism . detachTangentVector . flat . joinTangentPair h $ breakPoint dh0
           in (joinNeuralNetwork df dg, phts)

-- Convolutional Manifolds --

instance ( 1 <= (r*c), Manifold m, Manifold n, KnownNat r, KnownNat c, KnownNat rd
         , KnownNat (Div (Dimension n) (r*c)) , KnownNat (Div (Dimension m) (r*c)) )
  => Manifold (Convolutional rd r c m n) where
      type Dimension (Convolutional rd r c m n)
        = (Div (Dimension m) (r * c) * ((Div (Dimension n) (r * c) * ((2 * rd) + 1)) * ((2 * rd) + 1)))


instance KnownConvolutional rd r c m n => Map Mean Natural (Convolutional rd r c) m n where
      {-# INLINE (>.>) #-}
      (>.>) = convolveApply
      {-# INLINE (>$>) #-}
      (>$>) cnv = map (convolveApply cnv)

instance KnownConvolutional rd r c m n => Bilinear (Convolutional rd r c) m n where
    {-# INLINE (<.<) #-}
    (<.<) = convolveTransposeApply
    {-# INLINE (<$<) #-}
    (<$<) xs f = map (`convolveTransposeApply` f) xs
    {-# INLINE (>.<) #-}
    (>.<) = convolutionalOuterProduct
    {-# INLINE transpose #-}
    transpose = convolveTranspose

instance KnownConvolutional rd r c m n => Propagate Mean Natural (Convolutional rd r c) m n where
    {-# INLINE propagate #-}
    propagate = convolvePropagate

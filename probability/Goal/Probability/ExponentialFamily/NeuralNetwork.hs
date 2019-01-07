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
data HiddenNeuralNetwork (fs :: [Type -> Type -> Type]) (ys :: [Type]) z x

-- | A 'NeuralNetwork' defined by a type-list of transformations and a type-list of layers.
type NeuralNetwork fs ys = HiddenNeuralNetwork fs (Init (Tail ys)) (Head ys) (Last ys)

froysingleLayerNetwork :: c # HiddenNeuralNetwork '[f] '[] z x -> c # f z x
{-# INLINE froysingleLayerNetwork #-}
froysingleLayerNetwork = breakPoint

toSingleLayerNetwork :: c # f z x -> c # HiddenNeuralNetwork '[f] '[] z x
{-# INLINE toSingleLayerNetwork #-}
toSingleLayerNetwork = breakPoint

-- | Seperates a 'NeuralNetwork' into the final layer and the rest of the network.
splitNeuralNetwork
    :: (Manifold (f z y), Manifold (HiddenNeuralNetwork fs ys y x))
    => c # HiddenNeuralNetwork (f : fs) (y : ys) z x
    -> (c # f z y, c # HiddenNeuralNetwork fs ys y x)
{-# INLINE splitNeuralNetwork #-}
splitNeuralNetwork (Point xs) =
    let (xys,xns) = S.splitAt xs
     in (Point xys, Point xns)

-- | Joins a layer onto a 'NeuralNetwork' into a deeper network.
joinNeuralNetwork
    :: (Manifold (f z y), Manifold (HiddenNeuralNetwork fs ys y x))
    => c # f z y
    -> c # HiddenNeuralNetwork fs ys y x
    -> c # HiddenNeuralNetwork (f : fs) (y : ys) z x
{-# INLINE joinNeuralNetwork #-}
joinNeuralNetwork (Point xys) (Point xns) =
    Point $ xys S.++ xns

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
    => Mean #> Natural # Convolutional rd r c z x
    -> Mean # x
    -> S.Matrix (Div (Dimension x) (r*c)) (r*c) Double
{-# INLINE inputToImage #-}
inputToImage _ (Point img) = G.Matrix img

outputToImage
    :: (KnownConvolutional rd r c z x)
    => chrt1 #> chrt2 # Convolutional rd r c z x
    -> Dual chrt2 # z
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
    :: forall rd r c z x
    . KnownConvolutional rd r c z x
    => Mean #> Natural # Convolutional rd r c z x
    -> Mean # x
    -> Natural # z
{-# INLINE convolveApply #-}
convolveApply cnv imp =
    let img :: S.Matrix (Div (Dimension x) (r*c)) (r*c) Double
        img = inputToImage cnv imp
        krns :: S.Matrix (Div (Dimension z) (r*c)) (Div (Dimension x) (r*c) * (2*rd+1)*(2*rd+1)) Double
        krns = layerToKernels cnv
     in Point . G.toVector
         $ S.crossCorrelate2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolveTranspose
    :: forall chrt1 chrt2 rd r c z x
    . KnownConvolutional rd r c z x
    => chrt1 #> chrt2 # Convolutional rd r c z x
    -> Dual chrt2 #> Dual chrt1 # Convolutional rd r c x z
{-# INLINE convolveTranspose #-}
convolveTranspose cnv =
    let krns = layerToKernels cnv
        pnk = Proxy :: Proxy (Div (Dimension z) (r*c))
        pmd = Proxy :: Proxy (Div (Dimension x) (r*c))
        krn' :: S.Matrix (Div (Dimension x) (r*c)) (Div (Dimension z) (r*c)*(2*rd+1)*(2*rd+1)) Double
        krn' = S.kernelTranspose pnk pmd (Proxy @ rd) (Proxy @ rd) krns
     in Point $ G.toVector krn'

convolveTransposeApply
    :: forall chrt1 chrt2 rd r c z x
    . KnownConvolutional rd r c z x
    => Dual chrt2 # z
    -> chrt1 #> chrt2 # Convolutional rd r c z x
    -> Dual chrt1 # x
{-# INLINE convolveTransposeApply #-}
convolveTransposeApply imp cnv =
    let img = outputToImage cnv imp
        krns = layerToKernels cnv
     in Point . G.toVector
         $ S.convolve2d (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) krns img

convolutionalOuterProduct
    :: forall chrt1 chrt2 rd r c z x
    . KnownConvolutional rd r c z x
      => chrt2 # z
      -> chrt1 # x
      -> Dual chrt1 #> chrt2 # Convolutional rd r c z x
{-# INLINE convolutionalOuterProduct #-}
convolutionalOuterProduct (Point oimg) (Point iimg) =
    let omtx :: S.Matrix (Div (Dimension z) (r*c)) (r*c) Double
        omtx = G.Matrix oimg
        imtx :: S.Matrix (Div (Dimension x) (r*c)) (r*c) Double
        imtx = G.Matrix iimg
     in Point . G.toVector $ S.kernelOuterProduct (Proxy @ rd) (Proxy @ rd) (Proxy @ r) (Proxy @ c) omtx imtx

convolvePropagate
    :: forall rd r c z x . KnownConvolutional rd r c z x
      => [Mean # z]
      -> [Mean # x]
      -> Mean #> Natural # Convolutional rd r c z x
      -> (Natural #> Mean # Convolutional rd r c z x, [Natural # z])
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

instance Manifold (f z x) => Manifold (HiddenNeuralNetwork '[f] '[] z x) where
      type Dimension (HiddenNeuralNetwork '[f] '[] z x) = Dimension (f z x)

instance (Manifold (f z y), Manifold (HiddenNeuralNetwork fs ys y x))
  => Manifold (HiddenNeuralNetwork (f : fs) (y : ys) z x) where
      type Dimension (HiddenNeuralNetwork (f : fs) (y : ys) z x)
        = Dimension (f z y) + Dimension (HiddenNeuralNetwork fs ys y x)

instance (Map c d f z x) => Map c d (HiddenNeuralNetwork '[f] '[]) z x where
    {-# INLINE (>.>) #-}
    (>.>) f x = froysingleLayerNetwork f >.> x
    {-# INLINE (>$>) #-}
    (>$>) f xs = froysingleLayerNetwork f >$> xs

instance (Map c d f z y, Map c d (HiddenNeuralNetwork fs ys) y x, Transition d c y)
  => Map c d (HiddenNeuralNetwork (f : fs) (y : ys)) z x where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitNeuralNetwork fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitNeuralNetwork fg
         in f >$> map transition (g >$> xs)

instance (Propagate Mean Natural f z x) => Propagate Mean Natural (HiddenNeuralNetwork '[f] '[]) z x where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df,ps) = propagate dps qs $ froysingleLayerNetwork f
         in (toSingleLayerNetwork df,ps)

instance {-# OVERLAPPABLE #-}
    ( Propagate Mean Natural f z y, Propagate Mean Natural (HiddenNeuralNetwork fs ys) y x
    , Transition Natural Mean y, Legendre Natural y, Riemannian Natural y, Bilinear f z y)
  => Propagate Mean Natural (HiddenNeuralNetwork (Affine f : fs) (y : ys)) z x where
      propagate dzs xs fg =
          let (f,g) = splitNeuralNetwork fg
              fmtx = snd $ splitAffine f
              mys = dualTransition <$> ys
              (df,zhts) = propagate dzs mys f
              (dg,ys) = propagate dys xs g
              dys0 = dzs <$< fmtx
              dys = do
                  (y,dy0) <- zip ys dys0
                  return . dualIsomorphism . detachTangentVector . flat . joinTangentPair y $ breakPoint dy0
           in (joinNeuralNetwork df dg, zhts)

-- Convolutional Manifolds --

instance ( 1 <= (r*c), Manifold x, Manifold y, KnownNat r, KnownNat c, KnownNat rd
         , KnownNat (Div (Dimension x) (r*c)) , KnownNat (Div (Dimension y) (r*c)) )
  => Manifold (Convolutional rd r c y x) where
      type Dimension (Convolutional rd r c y x)
        = (Div (Dimension y) (r * c) * ((Div (Dimension x) (r * c) * ((2 * rd) + 1)) * ((2 * rd) + 1)))


instance KnownConvolutional rd r c z x => Map Mean Natural (Convolutional rd r c) z x where
      {-# INLINE (>.>) #-}
      (>.>) = convolveApply
      {-# INLINE (>$>) #-}
      (>$>) cnv = map (convolveApply cnv)

instance KnownConvolutional rd r c z x => Bilinear (Convolutional rd r c) z x where
    {-# INLINE (<.<) #-}
    (<.<) = convolveTransposeApply
    {-# INLINE (<$<) #-}
    (<$<) xs f = map (`convolveTransposeApply` f) xs
    {-# INLINE (>.<) #-}
    (>.<) = convolutionalOuterProduct
    {-# INLINE transpose #-}
    transpose = convolveTranspose

instance KnownConvolutional rd r c z x => Propagate Mean Natural (Convolutional rd r c) z x where
    {-# INLINE propagate #-}
    propagate = convolvePropagate

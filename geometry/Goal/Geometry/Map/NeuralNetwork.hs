{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation. The core type is the
-- 'NeuralNetwork', which is defined by two type-lists. The first type list is a
-- collection of maps. The second type-list is a list of 'Manifold's, which
-- defines the size and activation function of each layer of the network.
module Goal.Geometry.Map.NeuralNetwork
    ( -- * Neural Networks
      NeuralNetwork
    , splitNeuralNetwork
    , joinNeuralNetwork
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Vector
import Goal.Geometry.Map.Linear
import Goal.Geometry.Differential

import qualified Goal.Core.Vector.Storable as S


--- Multilayer ---


-- | A multilayer, artificial neural network.
data NeuralNetwork (gys :: [(Type -> Type -> Type,Type)]) (f :: (Type -> Type -> Type)) z x


fromSingleLayerNetwork :: c # NeuralNetwork '[] f z x -> c # f z x
{-# INLINE fromSingleLayerNetwork #-}
fromSingleLayerNetwork = breakPoint

toSingleLayerNetwork :: c # f z x -> c # NeuralNetwork '[] f z x
{-# INLINE toSingleLayerNetwork #-}
toSingleLayerNetwork = breakPoint

-- | Seperates a 'NeuralNetwork' into the final layer and the rest of the network.
splitNeuralNetwork
    :: (Manifold (f z y), Manifold (NeuralNetwork gys g y x))
    => c # NeuralNetwork ('(g,y):gys) f z x
    -> (c # f z y, c # NeuralNetwork gys g y x)
{-# INLINE splitNeuralNetwork #-}
splitNeuralNetwork (Point xs) =
    let (xys,xns) = S.splitAt xs
     in (Point xys, Point xns)

-- | Joins a layer onto the end of a 'NeuralNetwork'.
joinNeuralNetwork
    :: (Manifold (f z y), Manifold (NeuralNetwork gys g y x))
    => c # f z y
    -> c # NeuralNetwork gys g y x
    -> c # NeuralNetwork ('(g,y):gys) f z x
{-# INLINE joinNeuralNetwork #-}
joinNeuralNetwork (Point xys) (Point xns) =
    Point $ xys S.++ xns


--- Instances ---


-- General Neural Networks --

instance Manifold (f z x) => Manifold (NeuralNetwork '[] f z x) where
      type Dimension (NeuralNetwork '[] f z x) = Dimension (f z x)

instance (Manifold (f z y), Manifold (NeuralNetwork gys g y x))
  => Manifold (NeuralNetwork ('(g,y) : gys) f z x) where
      type Dimension (NeuralNetwork ('(g,y) : gys) f z x)
        = Dimension (f z y) + Dimension (NeuralNetwork gys g y x)

instance Map c f z x => Map c (NeuralNetwork '[] f) z x where
    {-# INLINE (>.>) #-}
    (>.>) f x = fromSingleLayerNetwork f >.> x
    {-# INLINE (>$>) #-}
    (>$>) f xs = fromSingleLayerNetwork f >$> xs

instance (Map c f z y, Map c (NeuralNetwork gys g) y x, Transition c (Dual c) y)
  => Map c (NeuralNetwork ('(g,y) : gys) f) z x where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitNeuralNetwork fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitNeuralNetwork fg
         in f >$> map transition (g >$> xs)

instance (Propagate c f z x) => Propagate c (NeuralNetwork '[] f) z x where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df,ps) = propagate dps qs $ fromSingleLayerNetwork f
         in (toSingleLayerNetwork df,ps)

instance
    ( Propagate c f z y, Propagate c (NeuralNetwork gys g) y x, Map c f y z
    , Transition c (Dual c) y, Legendre y, Riemannian c y, Bilinear f z y)
  => Propagate c (NeuralNetwork ('(g,y) : gys) (Affine f z)) z x where
      {-# INLINE propagate #-}
      propagate dzs xs fg =
          let (f,g) = splitNeuralNetwork fg
              fmtx = snd $ split f
              mys = transition <$> ys
              (df,zhts) = propagate dzs mys f
              (dg,ys) = propagate dys xs g
              dys0 = dzs <$< fmtx
              dys = zipWith flat ys dys0
           in (joinNeuralNetwork df dg, zhts)

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation. The core type is the
-- 'NeuralNetwork', which is defined by two type-lists. The first type list is a
-- collection of maps. The second type-list is a list of 'Manifold's, which
-- defines the size and activation function of each layer of the network.
module Goal.Geometry.Map.NeuralNetwork
    ( -- * Neural Networks
      NeuralNetwork (NeuralNetwork)
    , type (<<*)
    ) where


--- Imports ---


-- Goal --

import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Vector
import Goal.Geometry.Map.Linear
import Goal.Geometry.Differential


--- Multilayer ---


-- | A multilayer, artificial neural network.
newtype NeuralNetwork f g y z x = NeuralNetwork (f z y, g y x)

deriving instance (Manifold (f z y), Manifold (g y x)) => Manifold (NeuralNetwork f g y z x)
deriving instance (Manifold (f z y), Manifold (g y x)) => Product (NeuralNetwork f g y z x)

type (<<*) f g = NeuralNetwork f g
infixr 6 <<*


--- Instances ---


-- General Neural Networks --

instance (Map c f z y, Map c g y x, Transition c (Dual c) y)
  => Map c (NeuralNetwork f g y) z x where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = split fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = split fg
         in f >$> map transition (g >$> xs)

instance
    ( Propagate c (Affine f z0) z y, Propagate c g y x, Transition c (Dual c) y
    , Legendre y, Riemannian c y, Bilinear f z0 y, Map c f y z0, Translation z z0 )
  => Propagate c (NeuralNetwork (Affine f z0) g y) z x where
      {-# INLINE propagate #-}
      propagate dzs xs fg =
          let (f,g) = split fg
              fmtx = snd $ split f
              mys = transition <$> ys
              (df,zhts) = propagate dzs mys f
              (dg,ys) = propagate dys xs g
              dys0 = (anchor <$> dzs) <$< fmtx
              dys = zipWith flat ys dys0
           in (join df dg, zhts)

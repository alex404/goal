{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons which instantiate backpropagation through laziness.
-- They are fast (for CPUs), but right now they're type level description is
-- very ugly, and needs to be simplified.
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

deriving instance (Manifold (f z y), Manifold (g y x))
  => Manifold (NeuralNetwork f g y z x)
deriving instance (Manifold (f z y), Manifold (g y x))
  => Product (NeuralNetwork f g y z x)

-- | A first attempt at simplifying ugly 'NeuralNetwork' types.
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

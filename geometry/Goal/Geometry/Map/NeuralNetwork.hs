{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons which instantiate backpropagation through laziness.
-- They are fast (for CPUs), but right now they're type level description is
-- very ugly, and needs to be simplified.
module Goal.Geometry.Map.NeuralNetwork
    ( -- * Neural Networks
      NeuralNetwork
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
data NeuralNetwork (gys :: [(Type -> Type -> Type,Type)])
    (f :: (Type -> Type -> Type)) z x


--- Instances ---


instance Manifold (f z x) => Manifold (NeuralNetwork '[] f z x) where
      type Dimension (NeuralNetwork '[] f z x) = Dimension (f z x)

instance (Manifold (f z y), Manifold (NeuralNetwork gys g y x))
  => Manifold (NeuralNetwork ('(g,y) : gys) f z x) where
      type Dimension (NeuralNetwork ('(g,y) : gys) f z x)
        = Dimension (f z y) + Dimension (NeuralNetwork gys g y x)


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

instance (Manifold (f z y), Manifold (NeuralNetwork gys g y x))
  => Product (NeuralNetwork ('(g,y) : gys) f z x) where
      type First (NeuralNetwork ('(g,y) : gys) f z x)
        = f z y
      type Second (NeuralNetwork ('(g,y) : gys) f z x)
        = NeuralNetwork gys g y x
      join = joinNeuralNetwork
      split = splitNeuralNetwork

instance (Map c f z y, Map c (NeuralNetwork gys g) y x, Transition c (Dual c) y)
  => Map c (NeuralNetwork ('(g,y) : gys) f) z x where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = split fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = split fg
         in f >$> map transition (g >$> xs)

instance Map c f z x => Map c (NeuralNetwork '[] f) z x where
    {-# INLINE (>.>) #-}
    (>.>) f x = fromSingleLayerNetwork f >.> x
    {-# INLINE (>$>) #-}
    (>$>) f xs = fromSingleLayerNetwork f >$> xs

instance (Propagate c f z x) => Propagate c (NeuralNetwork '[] f) z x where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df,ps) = propagate dps qs $ fromSingleLayerNetwork f
         in (toSingleLayerNetwork df,ps)

instance
    ( Propagate c f z y, Propagate c (NeuralNetwork gys g) y x, Map c f y z
    , Transition c (Dual c) y, Legendre y, Riemannian c y, Bilinear f z y)
  => Propagate c (NeuralNetwork ('(g,y) : gys) (Affine f)) z x where
      {-# INLINE propagate #-}
      propagate dzs xs fg =
          let (f,g) = splitNeuralNetwork fg
              fmtx = snd $ splitAffine f
              mys = transition <$> ys
              (df,zhts) = propagate dzs mys f
              (dg,ys) = propagate dys xs g
              dys0 = dzs <$< fmtx
              dys = zipWith flat ys dys0
           in (joinNeuralNetwork df dg, zhts)


---- General Neural Networks --
--
--instance (Map c f z y, Map c g y x, Transition c (Dual c) y)
--  => Map c (NeuralNetwork f g y) z x where
--    {-# INLINE (>.>) #-}
--    (>.>) fg x =
--        let (f,g) = split fg
--         in f >.> transition (g >.> x)
--    {-# INLINE (>$>) #-}
--    (>$>) fg xs =
--        let (f,g) = split fg
--         in f >$> map transition (g >$> xs)
--
--instance
--    ( Propagate c (Affine f z0) z y, Propagate c g y x, Transition c (Dual c) y
--    , Legendre y, Riemannian c y, Bilinear f z0 y, Map c f y z0, Translation z z0 )
--  => Propagate c (NeuralNetwork (Affine f z0) g y) z x where
--      {-# INLINE propagate #-}
--      propagate dzs xs fg =
--          let (f,g) = split fg
--              fmtx = snd $ split f
--              mys = transition <$> ys
--              (df,zhts) = propagate dzs mys f
--              (dg,ys) = propagate dys xs g
--              dys0 = (anchor <$> dzs) <$< fmtx
--              dys = zipWith flat ys dys0
--           in (join df dg, zhts)
--



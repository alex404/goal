{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

{- | Multilayer perceptrons which instantiate backpropagation through laziness.
Right now the structure is simplier than it could be, but it leads to nice
types. If anyone ever wants to use a DNN with super-Affine biases, the code
is willing.
-}
module Goal.Geometry.Map.NeuralNetwork (
    -- * Neural Networks
    NeuralNetwork,
) where

--- Imports ---

-- Goal --

import Goal.Core

import Goal.Geometry.Differential
import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Map.Linear
import Goal.Geometry.Vector

import Goal.Core.Vector.Storable qualified as S
import Goal.Core.Vector.Storable.Linear qualified as L

--- Multilayer ---

-- | A multilayer, artificial neural network.
data NeuralNetwork (s :: L.LinearRep) (tys :: [(L.LinearRep, Type)]) z x

--- Instances ---

instance (Manifold (Affine t z z x)) => Manifold (NeuralNetwork t '[] z x) where
    type Dimension (NeuralNetwork t '[] z x) = Dimension (Affine t z z x)

instance
    (Manifold (Affine s z z y), Manifold (NeuralNetwork t tys y x)) =>
    Manifold (NeuralNetwork s ('(t, y) : tys) z x)
    where
    type
        Dimension (NeuralNetwork s ('(t, y) : tys) z x) =
            Dimension (Affine s z z y) + Dimension (NeuralNetwork t tys y x)

fromSingleLayerNetwork :: c # NeuralNetwork t '[] z x -> c # Affine t z z x
{-# INLINE fromSingleLayerNetwork #-}
fromSingleLayerNetwork = breakManifold

toSingleLayerNetwork :: c # Affine t z z x -> c # NeuralNetwork t '[] z x
{-# INLINE toSingleLayerNetwork #-}
toSingleLayerNetwork = breakManifold

-- | Seperates a 'NeuralNetwork' into the final layer and the rest of the network.
splitNeuralNetwork ::
    (Manifold (Affine s z z y), Manifold (NeuralNetwork t tys y x)) =>
    c # NeuralNetwork s ('(t, y) : tys) z x ->
    (c # Affine s z z y, c # NeuralNetwork t tys y x)
{-# INLINE splitNeuralNetwork #-}
splitNeuralNetwork (Point xs) =
    let (xfs, xnets) = S.splitAt xs
     in (Point xfs, Point xnets)

-- | Joins a layer onto the end of a 'NeuralNetwork'.
joinNeuralNetwork ::
    (Manifold (Affine s z z y), Manifold (NeuralNetwork t tys y x)) =>
    c # Affine s z z y ->
    c # NeuralNetwork t tys y x ->
    c # NeuralNetwork s ('(t, y) : tys) z x
{-# INLINE joinNeuralNetwork #-}
joinNeuralNetwork (Point xfs) (Point xnets) =
    Point $ xfs S.++ xnets

instance
    (Manifold (Affine s z z y), Manifold (NeuralNetwork t tys y x)) =>
    Product (NeuralNetwork s ('(t, y) : tys) z x)
    where
    type
        First (NeuralNetwork s ('(t, y) : tys) z x) =
            Affine s z z y
    type
        Second (NeuralNetwork s ('(t, y) : tys) z x) =
            NeuralNetwork t tys y x
    join = joinNeuralNetwork
    split = splitNeuralNetwork

instance
    (Manifold (Affine s z z y), KnownLinear s z y, Map c (NeuralNetwork t tys) y x, Transition c (Dual c) y) =>
    Map c (NeuralNetwork s ('(t, y) : tys)) z x
    where
    {-# INLINE (>.>) #-}
    (>.>) fnet x =
        let (f, net) = split fnet
         in f >.> transition (net >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fnet xs =
        let (f, net) = split fnet
         in f >$> map transition (net >$> xs)

instance (Manifold (Affine t z z x), KnownLinear t z x) => Map c (NeuralNetwork t '[]) z x where
    {-# INLINE (>.>) #-}
    (>.>) f x = fromSingleLayerNetwork f >.> x
    {-# INLINE (>$>) #-}
    (>$>) f xs = fromSingleLayerNetwork f >$> xs

instance (Manifold (Affine t z z x), KnownLinear t z x) => Propagate c (NeuralNetwork t '[]) z x where
    {-# INLINE propagate #-}
    propagate dps qs f =
        let (df, ps) = propagate dps qs $ fromSingleLayerNetwork f
         in (toSingleLayerNetwork df, ps)

instance
    ( Manifold (Affine s z z y)
    , KnownLinear s z y
    , KnownLinear s y z
    , Propagate c (NeuralNetwork t tys) y x
    , Transition c (Dual c) y
    , Legendre y
    , Riemannian c y
    ) =>
    Propagate c (NeuralNetwork s ('(t, y) : tys)) z x
    where
    {-# INLINE propagate #-}
    propagate dzs xs fg =
        let (f, g) = split fg
            fmtx = snd $ split f
            mys = transition <$> ys
            (df, zhts) = propagate dzs mys f
            (dg, ys) = propagate dys xs g
            dys0 = dzs <$< fmtx
            dys = zipWith flat ys dys0
         in (join df dg, zhts)

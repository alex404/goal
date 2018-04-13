{-# LANGUAGE UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      InterLayer
    , type (<*<)
    , type (:+:)
    , splitInterLayer
    , joinInterLayer
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S

--- Multilayer ---


data InterLayer f g

type (f :+: g) = InterLayer f g
infixr 3 :+:

type (m <*< g) = m <* Codomain g :+: g
infixr 3 <*<

splitInterLayer :: (Manifold m, Manifold n) => Point c (InterLayer m n) -> (Point c m, Point c n)
{-# INLINE splitInterLayer #-}
splitInterLayer (Point xs) =
    let (xms,xns) = S.splitAt xs
     in (Point xms, Point xns)

joinInterLayer :: (Manifold m, Manifold n) => Point c m -> Point c n -> Point c (InterLayer m n)
{-# INLINE joinInterLayer #-}
joinInterLayer (Point xms) (Point xns) =
    Point $ xms S.++ xns

instance (Manifold f, Manifold g) => Manifold (InterLayer f g) where
    type Dimension (InterLayer f g) = Dimension f + Dimension g

instance (Map f, Map g, Codomain g ~ Domain f) => Map (InterLayer f g) where
    type Domain (InterLayer f g) = Domain g
    type Codomain (InterLayer f g) = Codomain f

instance (d ~ Dual c, Apply c d f, Apply c d g, Transition d c (Codomain g), Codomain g ~ Domain f)
  => Apply c d (InterLayer f g) where
    {-# INLINE (>.>) #-}
    (>.>) fg x =
        let (f,g) = splitInterLayer fg
         in f >.> transition (g >.> x)
    {-# INLINE (>$>) #-}
    (>$>) fg xs =
        let (f,g) = splitInterLayer fg
         in f >$> mapReplicatedPoint transition (g >$> xs)

--instance (n ~ Codomain g, Manifold g, Manifold m, Propagate Mean Natural g, Legendre Natural (Codomain g), Riemannian Natural n)
--  => Propagate Mean Natural (InterLayer (Affine (Product m n)) g) where
--      propagate dps qs fg =
--          let (f,g) = splitInterLayer fg
--              fmtx = snd $ splitAffine f
--              mhs = mapReplicatedPoint dualTransition hs
--              (df,phts) = propagate dps mhs f
--              (dg,hs) = propagate dhs qs g
--              dhs = dualIsomorphism . detachTangentVector . flat
--                  . joinTangentPair hs . Point . coordinates $ dps <$< fmtx
----              dhs = joinReplicated . S.map (dualIsomorphism . detachTangentVector . flat)
----                  . S.zipWith joinTangentPair (splitReplicated hs) . S.map (Point . coordinates) . splitReplicated $ dps <$< fmtx
--           in (joinInterLayer df dg, phts)

instance (Replicated k Bernoulli ~ Codomain g, Manifold g, Propagate Mean Natural g, Legendre Natural (Codomain g), Manifold m, KnownNat k)
  => Propagate Mean Natural (InterLayer (Affine (Product m (Replicated k Bernoulli))) g) where
      propagate dps qs fg =
          let (f,g) = splitInterLayer fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              thts = S.map (\x -> x * (1-x)) $ coordinates mhs
              dhs = Point . S.zipWith (*) thts . coordinates $ dps <$< fmtx
           in (joinInterLayer df dg, phts)

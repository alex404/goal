{-# LANGUAGE BangPatterns,UndecidableInstances #-}

-- | Multilayer perceptrons and backpropagation.
module Goal.Probability.ExponentialFamily.NeuralNetwork
    ( -- * Neural Networks
      InterLayer
    , type (<*<)
    , type (:+:)
    , splitInterLayer
    , joinInterLayer
--    , Propagate2 (forwardPropagate,backwardPropagate)
--    , stochasticConditionalCrossEntropyDifferential2
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Storable as S

--- Multilayer ---


---- | A function for computing the relative entropy, also known as the KL-divergence.
--stochasticConditionalCrossEntropyDifferential2
--    :: ( Propagate2 f k
--       , ExponentialFamily (Domain f)
--       , ClosedFormExponentialFamily (Codomain f)
--       , KnownNat k)
--    => S.Vector k (Sample (Domain f))
--    -> S.Vector k (Sample (Codomain f))
--    -> Mean ~> Natural # f
--    -> CotangentVector (Mean ~> Natural) f
--{-# INLINE stochasticConditionalCrossEntropyDifferential2 #-}
--stochasticConditionalCrossEntropyDifferential2 xs ys f =
--    let !(!yhts,!qs) = forwardPropagate (S.map sufficientStatistic xs) f
--        mys = S.zipWith differentiator ys $!! yhts
--     in primalIsomorphism . backwardPropagate f mys $!! qs
--        where differentiator y yht =
--                  dualIsomorphism $ stochasticCrossEntropyDifferential (S.singleton y) yht


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

--class (1 <= n, Apply Mean Natural f, KnownNat n, NFData (Propagation f n)) => Propagate2 f n where
--    type Propagation f n :: *
--    forwardPropagate
--              :: S.Vector n (Mean # Domain f)
--              -> Function Mean Natural # f
--              -> (S.Vector n (Natural # Codomain f), Propagation f n)
--    backwardPropagate
--              :: Function Mean Natural # f
--              -> S.Vector n (Mean # Codomain f)
--              -> Propagation f n
--              -> Function Natural Mean # f
--
--instance (Apply Mean Natural (Product m n), KnownNat k, 1 <= k) => Propagate2 (Product m n) k where
--    type Propagation (Product m n) k = S.Vector k (Point Mean n)
--    {-# INLINE forwardPropagate #-}
--    forwardPropagate qs pq = (pq >$> qs, qs)
--    {-# INLINE backwardPropagate #-}
--    backwardPropagate _ dps qs = averagePoint $ S.zipWith (>.<) dps qs
--
--instance (Apply Mean Natural (Affine f), Propagate2 f k) => Propagate2 (Affine f) k where
--    type Propagation (Affine f) k = Propagation f k
--    {-# INLINE forwardPropagate #-}
--    forwardPropagate qs pq =
--        let (p,pq') = splitAffine pq
--            (p',qs') = forwardPropagate qs pq'
--         in (S.map (p <+>) p', qs')
--    {-# INLINE backwardPropagate #-}
--    backwardPropagate pq dps qs =
--        let (_,pq') = splitAffine pq
--         in joinAffine (averagePoint dps) $ backwardPropagate pq' dps qs
--
--
--instance (n ~ Codomain g, Manifold g, Manifold m, Propagate2 g k, Legendre Natural (Codomain g), Riemannian Natural n, Legendre Mean n)
--  => Propagate2 (InterLayer (Affine (Product m n)) g) k where
--      type Propagation (InterLayer (Affine (Product m n)) g) k =
--          (Propagation (Affine (Product m n)) k, Propagation g k)
--      {-# INLINE forwardPropagate #-}
--      forwardPropagate !qs !fg =
--          let !(!f,!g) = splitInterLayer fg
--              !(!hs,!qs') = forwardPropagate qs g
--              !mhs = S.map dualTransition hs
--              !(!ps,!qs0) = forwardPropagate mhs f
--           in (ps, (qs0,qs'))
--      {-# INLINE backwardPropagate #-}
--      backwardPropagate !fg !dps (!mhs,!qs) =
--          let !(!f,!g) = splitInterLayer fg
--              !df = backwardPropagate f dps mhs
--              !fmtx = snd $ splitAffine f
--              !dhs = S.map (dualIsomorphism . detachTangentVector . flat)
--                  $ S.zipWith joinTangentPair (S.map dualTransition mhs) (S.map (Point . coordinates) $ dps <$< fmtx)
--           in joinInterLayer df $ backwardPropagate g dhs qs


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


instance (n ~ Codomain g, Manifold g, Manifold m, Propagate Mean Natural g, Legendre Natural (Codomain g), Riemannian Natural n)
  => Propagate Mean Natural (InterLayer (Affine (Product m n)) g) where
      {-# INLINABLE propagate #-}
      propagate dps qs fg =
          let (f,g) = splitInterLayer fg
              fmtx = snd $ splitAffine f
              mhs = mapReplicatedPoint dualTransition hs
              (df,phts) = propagate dps mhs f
              (dg,hs) = propagate dhs qs g
              dhs = dualIsomorphism . detachTangentVector . flat
                  . joinTangentPair hs . Point . coordinates $ dps <$< fmtx
           in (joinInterLayer df dg, phts)

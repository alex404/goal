{-# LANGUAGE UndecidableInstances #-}

{- | Tools for modelling the differential and Riemannian geometry of a
'Manifold'.
-}
module Goal.Geometry.Differential (
    -- * Riemannian Manifolds
    Riemannian (metric, flat, sharp),
    euclideanDistance,

    -- * Backpropagation
    Propagate (propagate),
    backpropagation,

    -- * Legendre Manifolds
    PotentialCoordinates,
    Legendre (potential),
    DuallyFlat (dualPotential),
    canonicalDivergence,
) where

--- Imports ---

--- Goal

import Goal.Core

import Goal.Core.Vector.Storable qualified as S

import Goal.Geometry.Manifold
import Goal.Geometry.Map
import Goal.Geometry.Map.Linear
import Goal.Geometry.Vector

--- Misc

import Data.Kind (Type)

{- | A class of 'Map's which can 'propagate' errors. That is, given an error
derivative on the output, the input which caused the output, and a
'Map' to derive, return the derivative of the error with respect to the
parameters of the 'Map', as well as the output of the 'Map'.
-}
class (Map c f y x) => Propagate c f y x where
    propagate ::
        -- | The error differential
        [c #* y] ->
        -- | A vector of inputs
        [c #* x] ->
        -- | The function to differentiate
        c # f y x ->
        -- | The derivative, and function output
        (c #* f y x, [c # y])

-- | Distance between two 'Point's based on the 'Euclidean' metric (l2 distance).
euclideanDistance ::
    (Manifold x) => c # x -> c # x -> Double
{-# INLINE euclideanDistance #-}
euclideanDistance xs ys = S.l2Norm (coordinates $ xs - ys)

{- | An implementation of backpropagation using the 'Propagate' class. The first
argument is a function which takes a generalized target output and function
output and returns an error. The second argument is a list of target outputs
and function inputs. The third argument is the parameteric function to be
optimized, and its differential is what is returned.
-}
backpropagation ::
    (Propagate c f x y) =>
    (a -> c # x -> c #* x) ->
    [(a, c #* y)] ->
    c # f x y ->
    c #* f x y
{-# INLINE backpropagation #-}
backpropagation grd ysxs f =
    let (yss, xs) = unzip ysxs
        (df, yhts) = propagate dys xs f
        dys = zipWith grd yss yhts
     in df

--- Riemannian Manifolds ---

{- | 'Riemannian' 'Manifold's are differentiable 'Manifold's associated with a
smoothly varying 'Tensor' known as the Riemannian 'metric'. 'flat' and
'sharp' correspond to applying this 'metric' to elements of the 'Primal' and
'Dual' spaces, respectively.
-}
class (Primal c, Manifold x) => Riemannian c x where
    metric :: c # x -> c #* Tensor x x
    flat :: c # x -> c # x -> c #* x
    {-# INLINE flat #-}
    flat p v = metric p >.> v
    sharp :: c # x -> c #* x -> c # x
    {-# INLINE sharp #-}
    sharp p v = inverse (metric p) >.> v

--- Dually Flat Manifolds ---

{- | Although convex analysis is usually developed seperately from differential
geometry, it arises naturally out of the theory of dually flat 'Manifold's (<https://books.google.com/books?hl=en&lr=&id=vc2FWSo7wLUC&oi=fnd&pg=PR7&dq=methods+of+information+geometry&ots=4HsxHD_5KY&sig=gURe0tA3IEO-z-Cht_2TNsjjOG8#v=onepage&q=methods%20of%20information%20geometry&f=false Amari and Nagaoka, 2000>).

A 'Manifold' is 'Legendre' if it is associated with a particular convex
function known as a 'potential'.
-}
class (Primal (PotentialCoordinates x), Manifold x) => Legendre x where
    potential :: PotentialCoordinates x # x -> Double

{- | The (natural) coordinates of the given 'Manifold', on which the 'potential'
is defined.
-}
type family PotentialCoordinates x :: Type

{- | A 'Manifold' is 'DuallyFlat' when we can describe the 'dualPotential', which
is the convex conjugate of 'potential'.
-}
class (Legendre x) => DuallyFlat x where
    dualPotential :: PotentialCoordinates x #* x -> Double

{- | Computes the 'canonicalDivergence' between two points. Note that relative
to the typical definition of the KL-Divergence/relative entropy, the
arguments of this function are flipped.
-}
canonicalDivergence ::
    (DuallyFlat x) => PotentialCoordinates x # x -> PotentialCoordinates x #* x -> Double
{-# INLINE canonicalDivergence #-}
canonicalDivergence pp dq = potential pp + dualPotential dq - (pp <.> dq)

--- Instances ---

-- Euclidean --

instance (KnownNat k) => Riemannian Cartesian (Euclidean k) where
    {-# INLINE metric #-}
    metric _ =
        let diag :: Cartesian # Diagonal (Euclidean k)
            diag = 1
         in toTensor diag
    {-# INLINE flat #-}
    flat _ = breakChart
    {-# INLINE sharp #-}
    sharp _ = breakChart

-- Replicated Riemannian Manifolds --

-- instance {-# OVERLAPPABLE #-} (Riemannian c x, KnownNat k) => Riemannian c (Replicated k x) where
--    metric = error "Do not call metric on a replicated manifold"
--    {-# INLINE flat #-}
--    flat = S.map flat
--    {-# INLINE sharp #-}
--    sharp = S.map sharp

-- Backprop --

instance (KnownLinear t y x) => Propagate c (Linear t) y x where
    {-# INLINE propagate #-}
    propagate dps qs pq = (dps >$< qs, pq >$> qs)

-- instance (Bilinear Tensor y x, Primal c) => Propagate c Tensor y x where
--    {-# INLINE propagate #-}
--    propagate dps qs pq =
--        let foldfun (dp,q) (k,dpq) = (k+1,(dp >.< q) + dpq)
--         in (uncurry (/>) . foldr foldfun (0,0) $ zip dps qs, pq >$> qs)

instance
    (LinearSubspace y y0, KnownLinear t y0 x) =>
    Propagate c (Affine t y0) y x
    where
    {-# INLINE propagate #-}
    propagate dxs ys fxy =
        let (x, fx0y) = split fxy
            dx0s = linearProjection <$> dxs
            (dfx0y, x0s) = propagate dx0s ys fx0y
         in (join (average dxs) dfx0y, (x >+>) <$> x0s)

-- Sums --

type instance PotentialCoordinates (x, y) = PotentialCoordinates x

instance
    (Legendre x, Legendre y, PotentialCoordinates x ~ PotentialCoordinates y) =>
    Legendre (x, y)
    where
    {-# INLINE potential #-}
    potential pmn =
        let (pm, pn) = split pmn
         in potential pm + potential pn

type instance PotentialCoordinates (Replicated k x) = PotentialCoordinates x

instance (Legendre x, KnownNat k) => Legendre (Replicated k x) where
    {-# INLINE potential #-}
    potential ps =
        S.sum $ mapReplicated potential ps

instance (DuallyFlat x, KnownNat k) => DuallyFlat (Replicated k x) where
    {-# INLINE dualPotential #-}
    dualPotential ps =
        S.sum $ mapReplicated dualPotential ps

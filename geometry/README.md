In this package we find all the basic types and classes which drive the
manifold/geometry based approach of Goal. In particular, points and manifolds,
dual and multilinear spaces, function spaces and multilayer neural networks, and
generic optimization routines such as gradient pursuit. What follows is very
brief introduction to how we define points on a manifold in Goal.

The fundamental class in Goal is the `Manifold`
```haskell
class KnownNat (Dimension x) => Manifold x where
    type Dimension x :: Nat
```
`Manifold`s have an associated type, which is the `Dimension` of the `Manifold`.
The `Dimension` of a `Manifold` tells us the size required of vector to
represent a 'Point's on the given `Manifold`. In turn a `Point` is defined as
```haskell
newtype Point c x =
    Point { coordinates :: S.Vector (Dimension x) Double }
```
At the value level, a `Point` is a wrapper around an `S.Vector`, which is a
storable vector from the
[vector-sized](https://hackage.haskell.org/package/vector-sized) package, with
size `Dimension x`. In general, numerical operations in Goal are defined in
terms of [vector-sized](https://hackage.haskell.org/package/vector-sized) and
[hmatrix](https://hackage.haskell.org/package/hmatrix), with specific functions
for applying operations in bulk. Although I make no promises, Goal should be
quite efficient, at least for a CPU-based numerical library.

To continue, a `Point` is defined at the type-level by a `Manifold` `x`, and the
mysterious phantom type `c`.  In Goal `c` is referred to as a coordinate system,
or more succinctly as a chart.  A coordinate system describes how the abstract
elements of a `Manifold` may be uniquely represented by a vector of numbers. In
Goal we usually refer to `Point`s with the following infix type synonym
```haskell
type (c # x) = Point c x
```
which we may read as a `Point` in `c` coordinates on the `x` `Manifold`. I chose
the `#` symbol because it is reminiscent of the grid of a coordinate system.

Finally, with the notion of a coordinate system in hand, we may definition
`transition` functions for re-representing `Point`s in alternative coordinate
systems
```haskell
class Transition c d x where
    transition :: c # x -> d # x
```

As an example, where we define `Euclidean` space
```haskell
data Euclidean (n :: Nat)

instance (KnownNat n) => Manifold (Euclidean n) where
    type Dimension (Euclidean n) = n
```
and two coordinate systems on Euclidean space with an appropriate transition function
```haskell
data Cartesian
data Polar

instance Transition Cartesian Polar (Euclidean 2) where
    {-# INLINE transition #-}
    transition p =
        let [x,y] = listCoordinates p
            r = sqrt $ (x*x) + (y*y)
            phi = atan2 y x
         in fromTuple (r,phi)
```
we may create a `Point` in `Cartesian` coordinates an easily convert it to `Polar` coordinates
```haskell
xcrt :: Cartesian # Euclidean 2
xcrt = fromTuple (1,2)

xplr :: Polar # Euclidean 2
xplr = transition xcrt
```

So what has this bought us? Why would we make use of not only one, but
essentially two phantom types for describing vectors? Intuitively, the
`Manifold` under investigation is what we care about. If, for example, we
consider a `Manifold` of probability distributions, it is the distributions
themselves we care about. But distributions are abstract things, and so we
represent them in various coordinate systems (e.g. mean and variance) to handle
them numerically.

The charts available for a given `Manifold` are thus different (but isomorphic)
representations of the same thing. In particular, many coordinate systems have a
dual coordinate system that describes function differentials, which is critical
for numerical optimization. In general, many optimization problems can be
greatly simplified by finding the right coordinate system, and many complex
optimization problems can be solved by sequence of coordinate transformations
and simple numerical operations. Numerically the resulting computation is not
trivial, but theoretically it becomes an intuitive thing.

For in-depth tutorials visit my
[blog](https://sacha-sokoloski.gitlab.io/website/pages/blog.html).


# Goal: Geometric OptimizAtion Libraries

Goal (Geometric OptimizAtion Libraries) is a collection of Haskell libraries for
numerical optimization and machine learning. Building on the development of
vectors with static sizes, Goal furnishes vectors with additional type-level
structure based on ideas from differential geometry. The goals of Goal are to
provide types and functions that are:
- **practical**, and provide an effective yet simple means of expressing and
  applying optimization algorithms.
- *intuitive*, such that types can be easily read and understood, and
  meaningfully describe their values and functions.
- *evocative*, such that the correspondence between Goal types and mathematical
  constructs facilitates the user with an understanding of mathematical
  optimization and geometry.

The fundamental class in Goal is the `Manifold`:
```haskell
class KnownNat (Dimension x) => Manifold x where
    type Dimension x :: Nat
```
`Manifold`s have an associated type, which is the `Dimension` of the `Manifold`.
The `Dimension` of a `Manifold` tells us the size required of by vector to
represent a 'Point's on the given `Manifold`. In turn a `Point` is defined as:
```haskell
newtype Point c x =
    Point { coordinates :: S.Vector (Dimension x) Double }
```
At the value level, a `Point` is a wrapper around an `S.Vector`, which is a
storable vector form the @sized-vector@ package, with size `Dimension x`. At the
type level, a `Point` is defined by a `Manifold` `x`, and the mysterious phantom
type `c`.

In Goal `c` is referred to as a chart. A `chart` is a particular coordinate some
on the given `Manifold`, which assigns numbers to the abstract `Point`s of
`Manifold`. In Goal we usually refer to `Point`s with the following infix type
synonym:
```haskell
type (c # x) = Point c x
```
which we may read as a `Point` in `c` coordinates on the `x` `Manifold`. I chose
the `#` symbol because it is reminiscent of a little coordinate system.

Finally, with the notion of a coordinate system in hand, we may definition
`transition` functions for re-representing `Point`s in alternative coordinate
systems:
```haskell
class Transition c d x where
    transition :: c # x -> d # x
```

So what has this bought us? Why would we make use of not only one, but
essentially two phantom types for describing vectors? The answer lies of course
in how they are used, and the effectiveness of their use. So let us break this
down.

Goal is consists of the following packages.

## Core

*goal-core* provides a small set of generic libraries and data structures, and
re-exports a number of other libraries useful for scientific computing.

## Geometry

*goal-geometry* provides the basic types and classes which drive the manifold/geometry
based approach of Goal. Points and manifolds, multilinear and dual spaces,
function spaces and multilayer neural networks, and generic optimization
routines are defined here.

## Probability

*goal-probability* provides tools for implementing and applying machine learning
algorithms. The core concept of goal-probability are statistical manifolds, i.e.
manifold of probability distributions, with a focus on exponential family
distributions. Various graphical models are also defined here, e.g. mixture
models and restricted Boltzmann machines, as well as algorithms for fitting them
e.g. expectation maximization and contrastive divergence.

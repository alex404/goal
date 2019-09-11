# Goal: Geometric OptimizAtion Libraries

Goal (Geometric OptimizAtion Libraries) is a collection of Haskell libraries for
numerical optimization and machine learning. Building on the development of
vectors with static sizes, Goal furnishes vectors with additional type-level
structure based on ideas from differential geometry.

The fundamental class in Goal are `Manifold`s:
```haskell
class KnownNat (Dimension x) => Manifold x where
    type Dimension x :: Nat
```

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

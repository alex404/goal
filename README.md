# Goal: Geometric OptimizAtion Libraries

Goal (Geometric OptimizAtion Libraries) is a collection of Haskell libraries for
numerical optimization and machine learning. Expanding on vectors with static
sizes, Goal furnishes vectors with additional type-level structure based on
ideas from differential geometry. The goals of Goal are to provide types,
classes, and functions that are:
- **practical**, and provide a safe and effective means of expressing and
  applying efficient optimization algorithms.
- **intuitive**, such that types and classes can be easily read and understood,
  and meaningfully describe their values and functions.
- **evocative**, such that the correspondence between Goal types and
  mathematical constructs promotes an appreciation and understanding of
  mathematical optimization and geometry.

For more detailed explanations of the Goal libraries, visit the individual
package pages. At my
[blog](https://sacha-sokoloski.gitlab.io/website/pages/blog.html) you may find
examples and tutorials for Goal.

The central packages of the Goal libraries are:

## Core

[goal-core](https://gitlab.com/sacha-sokoloski/goal/tree/master/core)
re-exports a number of other libraries, and provides a set of additional
utility functions useful for scientific computing. In particular,
implementations of HLists and Mealy Automata (Circuits), tools for working with
CSV files and gnuplot, and a module which combines
[vector-sized](https://hackage.haskell.org/package/vector-sized) vectors with
[hmatrix](https://hackage.haskell.org/package/hmatrix).

## Geometry

[goal-geometry](https://gitlab.com/sacha-sokoloski/goal/tree/master/geometry)
provides the basic types and classes which drive the manifold/geometry based
approach of Goal. Points and manifolds, dual and multilinear spaces, function
spaces and multilayer neural networks, and generic optimization routines are
defined here.

## Probability

[goal-probability](https://gitlab.com/sacha-sokoloski/goal/tree/master/probability)
provides tools for implementing and applying basic statistical models. The core
concept of goal-probability are statistical manifolds, i.e.  manifold of
probability distributions, with a focus on exponential family distributions.

## Graphical

[goal-graphical](https://gitlab.com/sacha-sokoloski/goal/tree/master/graphical)
provides tools for with dynamical and graphical models. Various graphical models
are defined here, e.g. mixture models and restricted Boltzmann machines,
dynamical models such as HMMs and Kalman filters, and in both cases algorithms
for fitting them e.g.  expectation maximization and contrastive divergence
minimization.

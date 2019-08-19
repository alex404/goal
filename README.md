# Goal: Geometric OptimizAtion Libraries

A collection of Haskell libraries for numerical optimization and machine learning inspired by
differential and information geometry.

## Core

*goal-core* provides a small set of generic libraries and data structures, and
re-exports a number of other libraries useful for scientific computing.

## Geometry

*goal-geometry* provides the basic types and classes which drive the manifold/geometry
based approach of Goal. Points and manifold, multilinear and dual spaces, and
function spaces and multilayer neural networks, and generic optimization
routines are defined here.

## Probability

*goal-probability* provides tools for implementing and applying machine learning
algorithms. The core concept of goal-probability are statistical manifolds, i.e.
manifold of probability distributions, with a focus on exponential family
distributions. Various graphical models are also defined here, e.g. mixture
models and restricted Boltzmann machines, as well as algorithms for fitting them
e.g. expectation maximization and contrastive divergence.

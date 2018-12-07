# Goal: Geometric OptimizAtion Libraries

A collection of Haskell libraries for scientific computing and numerical
optimization. To import each package entirely, import Goal.Core, Goal.Geometry,
Goal.Probability, or Goal.Simulation.

## Core

*goal-core* provides a basic set of imports which which are generally useful for
scientific computing, as well as plotting functionality.

## Geometry

*goal-geometry* provides the basic types and classes which drive the manifold/geometry
based approach of this library. Function spaces and Multilinear algebra tools are
also defined here.

## Probability

*goal-probability* provides the tools to work with probability distributions and random
number generation. A key component of the approach taken here is the focus on
exponential family distributions; harmoniums (e.g. restricted Boltzmann
machines) and multi-layer perceptions are provided by this library as they are
implemented in terms of exponential families.

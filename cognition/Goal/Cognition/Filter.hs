{-# LANGUAGE Arrows #-}

-- | This set of modules represents the part of goal-cognition concerned with
-- solving the filtering problem. Different approaches to parametric filtering
-- are wrapped up in chains and flows over triples, where the triples are given
-- by the latent state, the observation, and the belief parameters.
module Goal.Cognition.Filter
    ( -- * Simulating filters
      parametricFilter
    , filterChain
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Simulation


--- Learning ---

-- | 'parametricFilter' returns a Bayesian filter (a 'Mealy' from observations to beliefs)
-- given a prediction function and an update function. The predictions and
-- beliefs are assumed to be part of the same parametric family of
-- distributions.
parametricFilter
    :: (Point c m x -> Point c m x) -- ^ Prediction function
    -> (z -> Point c m x -> Point c m x) -- ^ Update function
    -> Point c m x -- ^ Prior
    -> Circuit z (Point c m x) -- ^ Filter
parametricFilter pf uf prd0 = accumulateFunction prd0 $ \y prd ->
    let blf = uf y prd
        prd' = pf blf
     in (blf,prd')

-- | Constructs a 'Chain' over the latent states, observations, and belief
-- parameters, given a 'Chain' over the latent states, a 'Mealy' which
-- represents the emission distribution, and a 'Mealy' which computes updated
-- beliefs.
filterChain
    :: Chain x -- ^ Latent state chain
    -> Circuit x z -- ^ Emission distribution
    -> Circuit z tht -- ^ Filter
    -> Chain (x, z, tht) -- ^ coupled chain
filterChain xdnmcs nmdl zdnmcs = proc () -> do
    x <- xdnmcs -< ()
    n <- nmdl -< x
    z <- zdnmcs -< n
    returnA -< (x, n, z)

{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

-- | 'Statistical' models where the observations depend on known conditions.
module Goal.Probability.Conditional (
    SampleMap,

    -- ** Markov Kernels
    (>.>*),
    (>$>*),
    (*<.<),
    (*<$<),

    -- ** Conditional Distributions
    conditionalLogLikelihood,
    conditionalLogLikelihoodDifferential,
    conditionalDataMap,
    kFoldMap,
    kFoldMap',
    -- , mapToConditionalData
    mapConditionalLogLikelihood,
    mapConditionalLogLikelihoodDifferential,
    parMapConditionalLogLikelihood,
    parMapConditionalLogLikelihoodDifferential,
) where

--- Imports  ---

--- Goal

import Goal.Core
import Goal.Geometry

import Goal.Probability.ExponentialFamily
import Goal.Probability.Statistical

--- Misc

import Data.List qualified as L
import Data.Map.Strict qualified as M

import Control.Parallel.Strategies (parMap, rdeepseq)

--- Generic ---

-- | Evalutes the given conditional distribution at a 'SamplePoint'.
(>.>*) ::
    (Map Natural f y x, ExponentialFamily x) =>
    Natural # f y x ->
    SamplePoint x ->
    Natural # y
(>.>*) p x = p >.> sufficientStatistic x

-- | Mapped application of conditional distributions on a 'Sample'.
(>$>*) ::
    (Map Natural f y x, ExponentialFamily x) =>
    Natural # f y x ->
    Sample x ->
    [Natural # y]
(>$>*) p xs = p >$> (sufficientStatistic <$> xs)

infix 8 >.>*
infix 8 >$>*

-- | Applies the transpose of a 'Bilinear' 'Map' to a 'SamplePoint'.
(*<.<) ::
    (ExponentialFamily x, KnownLinear t x y, KnownLinear t y x) =>
    SamplePoint x ->
    Natural # Linear t x y ->
    Natural # y
(*<.<) x p = sufficientStatistic x <.< p

-- | Mapped transpose application on a 'Sample'.
(*<$<) ::
    (ExponentialFamily x, KnownLinear t x y, KnownLinear t y x) =>
    Sample x ->
    Natural # Linear t x y ->
    [Natural # y]
(*<$<) xs p = (sufficientStatistic <$> xs) <$< p

infix 8 *<.<
infix 8 *<$<

{- | A synonym for Maps from Inputs to Outputs that matches the confusing,
backwards style of Goal.
-}
type SampleMap z x = M.Map (SamplePoint x) (Sample z)

dependantLogLikelihood ::
    (LogLikelihood Natural y s, Map Natural f y x) =>
    [([s], Mean # x)] ->
    Natural # f y x ->
    Double
dependantLogLikelihood ysxs chrm =
    let (yss, xs) = unzip ysxs
     in average . zipWith logLikelihood yss $ chrm >$> xs

dependantLogLikelihoodDifferential ::
    (LogLikelihood Natural y s, Propagate Natural f y x) =>
    [([s], Mean # x)] ->
    Natural # f y x ->
    Mean # f y x
dependantLogLikelihoodDifferential ysxs chrm =
    let (yss, xs) = unzip ysxs
        (df, yhts) = propagate mys xs chrm
        mys = zipWith logLikelihoodDifferential yss yhts
     in df

dependantLogLikelihoodPar ::
    (LogLikelihood Natural y s, Map Natural f y x) =>
    [([s], Mean # x)] ->
    Natural # f y x ->
    Double
dependantLogLikelihoodPar ysxs chrm =
    let (yss, xs) = unzip ysxs
     in average . parMap rdeepseq (uncurry logLikelihood) . zip yss $ chrm >$> xs

dependantLogLikelihoodDifferentialPar ::
    (LogLikelihood Natural y s, Propagate Natural f y x) =>
    [([s], Mean # x)] ->
    Natural # f y x ->
    Mean # f y x
dependantLogLikelihoodDifferentialPar ysxs chrm =
    let (yss, xs) = unzip ysxs
        (df, yhts) = propagate mys xs chrm
        mys = parMap rdeepseq (uncurry logLikelihoodDifferential) $ zip yss yhts
     in df

{- | Turns a list of input/output pairs into a Map, by collecting into lists the
different outputs to each particular input.
-}
conditionalDataMap ::
    (Ord x) =>
    -- | Output/Input Pairs
    [(t, x)] ->
    -- | Input Output map
    M.Map x [t]
conditionalDataMap = L.foldl' folder M.empty
  where
    folder mp (t, x) =
        let ts = M.lookup x mp
            ts' = maybe [t] (t :) ts
         in M.insert x ts' mp

-- M.fromListWith (++) [(x, [y]) | (y, x) <- yxs]

{- | Partition a conditional dataset into k > 1 (training,validation) pairs,
where each dataset condition is partitioned to match its size.
-}
kFoldMap ::
    (Ord x) => Int -> M.Map x [y] -> [(M.Map x [y], M.Map x [y])]
kFoldMap k ixzmp =
    let ixzmps = kFold k <$> ixzmp
        ixs = M.keys ixzmp
        tvzss = M.elems ixzmps
        tvxzmps = M.fromList . zip ixs <$> L.transpose tvzss
     in zip (fmap fst <$> tvxzmps) (fmap snd <$> tvxzmps)

{- | Partition a conditional dataset into k > 2 (training,test,validation) triplets,
where each dataset condition is partitioned to match its size.
-}
kFoldMap' ::
    (Ord x) => Int -> M.Map x [y] -> [(M.Map x [y], M.Map x [y], M.Map x [y])]
kFoldMap' k ixzmp =
    let ixzmps = kFold' k <$> ixzmp
        ixs = M.keys ixzmp
        tvzss = M.elems ixzmps
        tvxzmps = M.fromList . zip ixs <$> L.transpose tvzss
     in zip3
            (fmap (\(x, _, _) -> x) <$> tvxzmps)
            (fmap (\(_, x, _) -> x) <$> tvxzmps)
            (fmap (\(_, _, x) -> x) <$> tvxzmps)

-- mapToConditionalData :: M.Map x [y] -> [(y,x)]
-- mapToConditionalData mp =
--    let (xs,zss) = unzip $ M.toAscList mp
--     in concat $ zipWith (\x zs -> zip zs $ repeat x) xs zss

-- | The conditional 'logLikelihood' for a conditional distribution.
conditionalLogLikelihood ::
    (ExponentialFamily x, Map Natural f y x, LogLikelihood Natural y t) =>
    -- | Output/Input Pairs
    [(t, SamplePoint x)] ->
    -- | Function
    Natural # f y x ->
    -- | conditional cross entropy estimate
    Double
conditionalLogLikelihood yxs f =
    let ysxs = [([y], sufficientStatistic x) | (y, x) <- yxs]
     in dependantLogLikelihood ysxs f

-- | The conditional 'logLikelihoodDifferential' for a conditional distribution.
conditionalLogLikelihoodDifferential ::
    (ExponentialFamily x, LogLikelihood Natural y t, Propagate Natural f y x) =>
    -- | Output/Input Pairs
    [(t, SamplePoint x)] ->
    -- | Function
    Natural # f y x ->
    -- | Differential
    Mean # f y x
conditionalLogLikelihoodDifferential yxs f =
    let ysxs = [([y], sufficientStatistic x) | (y, x) <- yxs]
     in dependantLogLikelihoodDifferential ysxs f

{- | The conditional 'logLikelihood' for a conditional distribution, where
redundant conditions/inputs are combined. This can dramatically increase performance when
the number of distinct conditions/inputs is small.
-}
mapConditionalLogLikelihood ::
    (ExponentialFamily x, Map Natural f y x, LogLikelihood Natural y t) =>
    -- | Output/Input Pairs
    M.Map (SamplePoint x) [t] ->
    -- | Function
    Natural # f y x ->
    -- | conditional cross entropy estimate
    Double
mapConditionalLogLikelihood xtsmp =
    dependantLogLikelihood [(ts, sufficientStatistic x) | (x, ts) <- M.toList xtsmp]

{- | The conditional 'logLikelihoodDifferential', where redundant conditions are
combined. This can dramatically increase performance when the number of
distinct conditions is small.
-}
mapConditionalLogLikelihoodDifferential ::
    ( ExponentialFamily x
    , LogLikelihood Natural y t
    , Propagate Natural f y x
    , Ord (SamplePoint x)
    ) =>
    -- | Output/Input Pairs
    M.Map (SamplePoint x) [t] ->
    -- | Function
    Natural # f y x ->
    -- | Differential
    Mean # f y x
mapConditionalLogLikelihoodDifferential xtsmp =
    dependantLogLikelihoodDifferential [(ts, sufficientStatistic x) | (x, ts) <- M.toList xtsmp]

{- | The conditional 'logLikelihood' for a conditional distribution, where
redundant conditions/inputs are combined. This can dramatically increase performance when
the number of distinct conditions/inputs is small.
-}
parMapConditionalLogLikelihood ::
    (ExponentialFamily x, Map Natural f y x, LogLikelihood Natural y t) =>
    -- | Output/Input Pairs
    M.Map (SamplePoint x) [t] ->
    -- | Function
    Natural # f y x ->
    -- | conditional cross entropy estimate
    Double
parMapConditionalLogLikelihood xtsmp =
    dependantLogLikelihoodPar [(ts, sufficientStatistic x) | (x, ts) <- M.toList xtsmp]

{- | The conditional 'logLikelihoodDifferential', where redundant conditions are
combined. This can dramatically increase performance when the number of
distinct conditions is small.
-}
parMapConditionalLogLikelihoodDifferential ::
    ( ExponentialFamily x
    , LogLikelihood Natural y t
    , Propagate Natural f y x
    , Ord (SamplePoint x)
    ) =>
    -- | Output/Input Pairs
    M.Map (SamplePoint x) [t] ->
    -- | Function
    Natural # f y x ->
    -- | Differential
    Mean # f y x
parMapConditionalLogLikelihoodDifferential xtsmp =
    dependantLogLikelihoodDifferentialPar [(ts, sufficientStatistic x) | (x, ts) <- M.toList xtsmp]

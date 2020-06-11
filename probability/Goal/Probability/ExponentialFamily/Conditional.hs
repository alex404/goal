{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances,Arrows #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Probability.ExponentialFamily.Conditional
    ( SampleMap
    -- ** Conditional Distributions
    , conditionalLogLikelihood
    , conditionalLogLikelihoodDifferential
    , conditionalDataMap
    , kFoldMap
    , mapToConditionalData
    , mapConditionalLogLikelihood
    , mapConditionalLogLikelihoodDifferential
     -- * Conditional Harmoniums
    , ConditionalBias
    , ConditionalBiases
    , ConditionalDeepHarmonium
    , ConditionalMixture
    -- ** Construction
    , joinConditionalDeepHarmonium
    , splitConditionalDeepHarmonium
    -- ** Evaluation
    , conditionalHarmoniumConjugationDifferential
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily

import qualified Goal.Core.Vector.Storable as S
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Inference

import qualified Data.Map as M
import qualified Data.List as L

--import Control.Parallel.Strategies


--- Generic ---


type SampleMap z x = M.Map (SamplePoint x) (Sample z)


dependantLogLikelihood
    :: (LogLikelihood d y s, Map Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Double
{-# INLINE dependantLogLikelihood #-}
dependantLogLikelihood ysxs chrm =
    let (yss,xs) = unzip ysxs
     in average . zipWith logLikelihood yss $ chrm >$> xs

dependantLogLikelihoodDifferential
    :: (LogLikelihood d y s, Propagate Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Function Mean d #* f y x
{-# INLINE dependantLogLikelihoodDifferential #-}
dependantLogLikelihoodDifferential ysxs chrm =
    let (yss,xs) = unzip ysxs
        (df,yhts) = propagate mys xs chrm
        mys = zipWith logLikelihoodDifferential yss yhts
     in df

--dependantLogLikelihood
--    :: (LogLikelihood d y s, Map Mean d f y x)
--    => [([s], Mean # x)] -> Function Mean d # f y x -> Double
--{-# INLINE dependantLogLikelihood #-}
--dependantLogLikelihood ysxs chrm =
--    let (yss,xs) = unzip ysxs
--     in average . parMap rdeepseq (uncurry logLikelihood) . zip yss $ chrm >$> xs
--
--dependantLogLikelihoodDifferential
--    :: (LogLikelihood d y s, Propagate Mean d f y x)
--    => [([s], Mean # x)] -> Function Mean d # f y x -> Function Mean d #* f y x
--{-# INLINE dependantLogLikelihoodDifferential #-}
--dependantLogLikelihoodDifferential ysxs chrm =
--    let (yss,xs) = unzip ysxs
--        (df,yhts) = propagate mys xs chrm
--        mys = parMap rdeepseq (uncurry logLikelihoodDifferential) $ zip yss yhts
--     in df

conditionalDataMap
    :: Ord x
    => [(t, x)] -- ^ Output/Input Pairs
    -> M.Map x [t] -- ^ Input Output map
{-# INLINE conditionalDataMap #-}
conditionalDataMap yxs =
    M.fromListWith (++) [(x, [y]) | (y, x) <- yxs]

kFoldMap
    :: Ord x => Int -> M.Map x [y] -> [(M.Map x [y], M.Map x [y])]
kFoldMap k ixzmp =
    let ixzmps = kFold k <$> ixzmp
        ixs = M.keys ixzmp
        tvzss = M.elems ixzmps
        tvxzmps = M.fromList . zip ixs <$> L.transpose tvzss
     in zip (fmap fst <$> tvxzmps) (fmap snd <$> tvxzmps)

mapToConditionalData :: M.Map x [y] -> [(y,x)]
mapToConditionalData mp =
    let (xs,zss) = unzip $ M.toAscList mp
     in concat $ zipWith (\x zs -> zip zs $ repeat x) xs zss


-- | The conditional 'logLikelihood' for a conditional distribution.
conditionalLogLikelihood
    :: (ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t)
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE conditionalLogLikelihood #-}
conditionalLogLikelihood yxs f =
    let ysxs = [ ([y],sufficientStatistic x) | (y,x) <- yxs ]
     in dependantLogLikelihood ysxs f

-- | The conditional 'logLikelihoodDifferential' for a conditional distribution.
conditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t, Propagate Mean Natural f y x )
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
{-# INLINE conditionalLogLikelihoodDifferential #-}
conditionalLogLikelihoodDifferential yxs f =
    let ysxs = [ ([y],sufficientStatistic x) | (y,x) <- yxs ]
     in dependantLogLikelihoodDifferential ysxs f

-- | The conditional 'logLikelihood' for a conditional distribution, where
-- redundant conditions/inputs are combined. This can dramatically increase performance when
-- the number of distinct conditions/inputs is small.
mapConditionalLogLikelihood
    :: ( ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t )
    => M.Map (SamplePoint x) [t] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
{-# INLINE mapConditionalLogLikelihood #-}
mapConditionalLogLikelihood xtsmp f =
     dependantLogLikelihood [ (ts, sufficientStatistic x) | (x,ts) <- M.toList xtsmp] f

-- | The conditional 'logLikelihoodDifferential', where redundant conditions are
-- combined. This can dramatically increase performance when the number of
-- distinct conditions is small.
mapConditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t
       , Propagate Mean Natural f y x, Ord (SamplePoint x) )
    => M.Map (SamplePoint x) [t] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
{-# INLINE mapConditionalLogLikelihoodDifferential #-}
mapConditionalLogLikelihoodDifferential xtsmp f =
     dependantLogLikelihoodDifferential [ (ts, sufficientStatistic x) | (x,ts) <- M.toList xtsmp] f

-- | An approximate differntial for conjugating a harmonium likelihood.
conditionalHarmoniumConjugationDifferential
    :: ( Propagate Mean Natural f y z
       , LegendreExponentialFamily (Harmonium y g x)
       , LegendreExponentialFamily x, ExponentialFamily y, ExponentialFamily z )
    => Double -- ^ Conjugation shift
    -> Natural # z -- ^ Conjugation parameters
    -> Sample z -- ^ Sample points
    -> Natural #> ConditionalHarmonium f y g x z
    -> Mean #> ConditionalHarmonium f y g x z
{-# INLINE conditionalHarmoniumConjugationDifferential #-}
conditionalHarmoniumConjugationDifferential rho0 rprms xsmps chrm =
    let rcts = conjugationCurve rho0 rprms xsmps
        mhrms = transition <$> nhrms
        ptns = potential <$> nhrms
        dhrms = [ (ptn - rct) .> mhrm | (rct,mhrm,ptn) <- zip3 rcts mhrms ptns ]
        (dchrm,nhrms) = propagate dhrms (sufficientStatistic <$> xsmps) chrm
     in dchrm



--- Types ---


-- | A generic conditional model.
data ConditionalBias (f :: Type -> Type -> Type) z x

-- | Another generic conditional model.
data ConditionalBiases (f :: Type -> Type -> Type) z x

-- | A conditional 'DeepHarmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalDeepHarmonium f y (gxs :: [(Type -> Type -> Type,Type)])
  = ConditionalBias f (DeepHarmonium y gxs)

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium f y g x = ConditionalDeepHarmonium f y ('[ '(g,x)])

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture f y k = ConditionalHarmonium f y Tensor (Categorical k) -- ^ Function

-- | Splits a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalDeepHarmonium
    :: (Manifold (f y z), Manifold (DeepHarmonium y gxs))
    => c #> ConditionalDeepHarmonium f y gxs z -- ^ Conditional Harmonium
    -> (c # DeepHarmonium y gxs, c #> f y z) -- ^ Matrix function and upper part
{-# INLINE splitConditionalDeepHarmonium #-}
splitConditionalDeepHarmonium dhrm =
    let (dhrmcs,fcs) = S.splitAt $ coordinates dhrm
     in (Point dhrmcs,Point fcs)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalDeepHarmonium
    :: ( Manifold (f y z), Manifold (DeepHarmonium y gxs) )
    => c # DeepHarmonium y gxs
    -> c #> f y z
    -> c #> ConditionalDeepHarmonium f y gxs z -- ^ Conditional Harmonium
{-# INLINE joinConditionalDeepHarmonium #-}
joinConditionalDeepHarmonium (Point dcs) (Point fcs) = Point $ dcs S.++ fcs


--- Instances ---


instance (Map Mean Natural f y z, Manifold (DeepHarmonium y gxs))
    => Manifold (ConditionalDeepHarmonium f y gxs z) where
        type Dimension (ConditionalDeepHarmonium f y gxs z)
          = Dimension (DeepHarmonium y gxs) + Dimension (f y z)

instance ( Map Mean Natural f y z, Manifold (DeepHarmonium y gxs) )
     => Map Mean Natural (ConditionalBias f) (DeepHarmonium y gxs) z where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm q =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in biasBottom (pq >.> q) dhrm
    {-# INLINE (>$>) #-}
    (>$>) pdhrm qs =
        let (dhrm,pq) = splitConditionalDeepHarmonium pdhrm
         in flip biasBottom dhrm <$> (pq >$> qs)

instance (Propagate Mean Natural f y z, Manifold (DeepHarmonium y gxs))
  => Propagate Mean Natural (ConditionalBias f) (DeepHarmonium y gxs) z where
        {-# INLINE propagate #-}
        propagate dhrms dzs chrm =
            let dys = getBottomBias <$> dhrms
                (hrm,f) = splitConditionalDeepHarmonium chrm
                (df,hrmhts) = propagate dys dzs f
             in (joinConditionalDeepHarmonium (average dhrms) df, flip biasBottom hrm <$> hrmhts)


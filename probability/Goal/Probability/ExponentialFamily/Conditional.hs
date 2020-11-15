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
    , kFoldMap'
    --, mapToConditionalData
    , mapConditionalLogLikelihood
    , mapConditionalLogLikelihoodDifferential
    , parMapConditionalLogLikelihood
    , parMapConditionalLogLikelihoodDifferential
     -- * Conditional Harmoniums
    , ConditionalBias
    , ConditionalBiases
    , ConditionalDeepHarmonium
    , ConditionalMixture
    -- ** Construction
    , joinConditionalDeepHarmonium
    , joinConditionalHarmonium
    , splitConditionalDeepHarmonium
    , splitConditionalHarmonium
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

import qualified Data.Map.Strict as M
import qualified Data.List as L

import Control.Parallel.Strategies


--- Generic ---


type SampleMap z x = M.Map (SamplePoint x) (Sample z)


dependantLogLikelihood
    :: (LogLikelihood d y s, Map Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Double
dependantLogLikelihood ysxs chrm =
    let (yss,xs) = unzip ysxs
     in average . zipWith logLikelihood yss $ chrm >$> xs

dependantLogLikelihoodDifferential
    :: (LogLikelihood d y s, Propagate Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Function Mean d #* f y x
dependantLogLikelihoodDifferential ysxs chrm =
    let (yss,xs) = unzip ysxs
        (df,yhts) = propagate mys xs chrm
        mys = zipWith logLikelihoodDifferential yss yhts
     in df

dependantLogLikelihoodPar
    :: (LogLikelihood d y s, Map Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Double
dependantLogLikelihoodPar ysxs chrm =
    let (yss,xs) = unzip ysxs
     in average . parMap rdeepseq (uncurry logLikelihood) . zip yss $ chrm >$> xs

dependantLogLikelihoodDifferentialPar
    :: (LogLikelihood d y s, Propagate Mean d f y x)
    => [([s], Mean # x)] -> Function Mean d # f y x -> Function Mean d #* f y x
dependantLogLikelihoodDifferentialPar ysxs chrm =
    let (yss,xs) = unzip ysxs
        (df,yhts) = propagate mys xs chrm
        mys = parMap rdeepseq (uncurry logLikelihoodDifferential) $ zip yss yhts
     in df

conditionalDataMap
    :: Ord x
    => [(t, x)] -- ^ Output/Input Pairs
    -> M.Map x [t] -- ^ Input Output map
conditionalDataMap yxs = foldl' folder M.empty yxs
    where folder mp (t,x) =
            let ts = M.lookup x mp
                ts' = maybe [t] (t:) ts
             in M.insert x ts' mp
    --M.fromListWith (++) [(x, [y]) | (y, x) <- yxs]

-- | Partition a conditional dataset into k > 1 (training,validation) pairs,
-- where each dataset condition is partitioned to match its size.
kFoldMap
    :: Ord x => Int -> M.Map x [y] -> [(M.Map x [y], M.Map x [y])]
kFoldMap k ixzmp =
    let ixzmps = kFold k <$> ixzmp
        ixs = M.keys ixzmp
        tvzss = M.elems ixzmps
        tvxzmps = M.fromList . zip ixs <$> L.transpose tvzss
     in zip (fmap fst <$> tvxzmps) (fmap snd <$> tvxzmps)

kFoldMap'
    :: Ord x => Int -> M.Map x [y] -> [(M.Map x [y], M.Map x [y], M.Map x [y])]
kFoldMap' k ixzmp =
    let ixzmps = kFold' k <$> ixzmp
        ixs = M.keys ixzmp
        tvzss = M.elems ixzmps
        tvxzmps = M.fromList . zip ixs <$> L.transpose tvzss
     in zip3 (fmap (\(x,_,_) -> x) <$> tvxzmps)
             (fmap (\(_,x,_) -> x) <$> tvxzmps)
             (fmap (\(_,_,x) -> x) <$> tvxzmps)

--mapToConditionalData :: M.Map x [y] -> [(y,x)]
--mapToConditionalData mp =
--    let (xs,zss) = unzip $ M.toAscList mp
--     in concat $ zipWith (\x zs -> zip zs $ repeat x) xs zss


-- | The conditional 'logLikelihood' for a conditional distribution.
conditionalLogLikelihood
    :: (ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t)
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
conditionalLogLikelihood yxs f =
    let ysxs = [ ([y],sufficientStatistic x) | (y,x) <- yxs ]
     in dependantLogLikelihood ysxs f

-- | The conditional 'logLikelihoodDifferential' for a conditional distribution.
conditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t, Propagate Mean Natural f y x )
    => [(t, SamplePoint x)] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
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
mapConditionalLogLikelihoodDifferential xtsmp f =
     dependantLogLikelihoodDifferential [ (ts, sufficientStatistic x) | (x,ts) <- M.toList xtsmp] f

-- | The conditional 'logLikelihood' for a conditional distribution, where
-- redundant conditions/inputs are combined. This can dramatically increase performance when
-- the number of distinct conditions/inputs is small.
parMapConditionalLogLikelihood
    :: ( ExponentialFamily x, Map Mean Natural f y x, LogLikelihood Natural y t )
    => M.Map (SamplePoint x) [t] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Double -- ^ conditional cross entropy estimate
parMapConditionalLogLikelihood xtsmp f =
     dependantLogLikelihoodPar [ (ts, sufficientStatistic x) | (x,ts) <- M.toList xtsmp] f

-- | The conditional 'logLikelihoodDifferential', where redundant conditions are
-- combined. This can dramatically increase performance when the number of
-- distinct conditions is small.
parMapConditionalLogLikelihoodDifferential
    :: ( ExponentialFamily x, LogLikelihood Natural y t
       , Propagate Mean Natural f y x, Ord (SamplePoint x) )
    => M.Map (SamplePoint x) [t] -- ^ Output/Input Pairs
    -> Natural #> f y x -- ^ Function
    -> Natural #*> f y x -- ^ Differential
parMapConditionalLogLikelihoodDifferential xtsmp f =
     dependantLogLikelihoodDifferentialPar [ (ts, sufficientStatistic x) | (x,ts) <- M.toList xtsmp] f

-- | An approximate differntial for conjugating a harmonium likelihood.
conditionalHarmoniumConjugationDifferential
    :: ( Propagate Mean Natural f y z, Manifold (g y x)
       , LegendreExponentialFamily (Harmonium y g x)
       , LegendreExponentialFamily x, ExponentialFamily y, ExponentialFamily z )
    => Double -- ^ Conjugation shift
    -> Natural # z -- ^ Conjugation parameters
    -> Sample z -- ^ Sample points
    -> Natural #> ConditionalHarmonium f y g x z
    -> Mean #> ConditionalHarmonium f y g x z
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
    :: ( Map Mean Natural f y z, Manifold (g y x), Manifold (DeepHarmonium x gxs) )
    => c #> ConditionalDeepHarmonium f y ('(g,x) : gxs) z -- ^ Conditional Harmonium
    -> (c #> f y z, c #> g y x, c # DeepHarmonium x gxs) -- ^ Matrix function and upper part
splitConditionalDeepHarmonium dhrm =
    let (fcs,gdhrmcs) = S.splitAt $ coordinates dhrm
        (gcs,dhrmcs) = S.splitAt gdhrmcs
     in (Point fcs,Point gcs,Point dhrmcs)

splitConditionalHarmonium
    :: ( Map Mean Natural f y z, Manifold (g y x), Manifold x )
    => c #> ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
    -> (c #> f y z, c #> g y x, c # x) -- ^ Matrix function and upper part
splitConditionalHarmonium dhrm =
    let (fyz,gyx,nx0) = splitConditionalDeepHarmonium dhrm
     in (fyz,gyx,fromOneHarmonium nx0)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalDeepHarmonium
    :: ( Map Mean Natural f y z, Manifold (g y x), Manifold (DeepHarmonium x gxs) )
    => c #> f y z
    -> c #> g y x
    -> c # DeepHarmonium x gxs -- ^ Matrix function and upper part
    -> c #> ConditionalDeepHarmonium f y ('(g,x) : gxs) z -- ^ Conditional Harmonium
joinConditionalDeepHarmonium (Point fcs) (Point gcs) (Point dhrmcs) =
    Point $ fcs S.++ gcs S.++ dhrmcs

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalHarmonium
    :: ( Map Mean Natural f y z, Manifold (g y x), Manifold x )
    => c #> f y z
    -> c #> g y x
    -> c # x -- ^ Matrix function and upper part
    -> c #> ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
joinConditionalHarmonium fyz gyx x =
    joinConditionalDeepHarmonium fyz gyx $ toOneHarmonium x


--- Instances ---


instance ( Map Mean Natural f y z, Manifold (g y x), Manifold (DeepHarmonium x gxs) )
  => Manifold (ConditionalDeepHarmonium f y ('(g,x) : gxs) z) where
      type Dimension (ConditionalDeepHarmonium f y ('(g,x) : gxs) z)
        = Dimension (f y z) + Dimension (g y x) + Dimension (DeepHarmonium x gxs)

instance ( Map Mean Natural f y z, Manifold (g y x), Manifold x
         , Manifold (DeepHarmonium x gxs) )
     => Map Mean Natural (ConditionalBias f) (DeepHarmonium y ('(g,x) : gxs)) z where
    (>.>) chrm mz =
        let (fyz,gyx,dhrm) = splitConditionalDeepHarmonium chrm
            affyx = joinAffine (fyz >.> mz) gyx
         in joinBottomHarmonium affyx dhrm
    (>$>) chrm mz =
        let (fyz,gyx,dhrm) = splitConditionalDeepHarmonium chrm
            affyxs = flip joinAffine gyx <$> (fyz >$> mz)
         in flip joinBottomHarmonium dhrm <$> affyxs

instance ( Propagate Mean Natural f y z, Manifold (g y x), Manifold x
         , Manifold (DeepHarmonium x gxs) )
     => Propagate Mean Natural (ConditionalBias f) (DeepHarmonium y ('(g,x) : gxs)) z where
        propagate mdhrms mzs chrm =
            let (mys,mgyxs,mdhrms') = unzip3 $ splitDeepHarmonium <$> mdhrms
                (nfyz,ngyx,ndhrm) = splitConditionalDeepHarmonium chrm
                (mfyz,nyhts) = propagate mys mzs nfyz
             in ( joinConditionalDeepHarmonium mfyz (average mgyxs) (average mdhrms')
                , [joinDeepHarmonium nyht ngyx ndhrm | nyht <- nyhts] )


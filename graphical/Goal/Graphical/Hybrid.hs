{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Hybrid
    (
     -- * Conditional Harmoniums
      ConditionalBias
    , ConditionalBiases
    , ConditionalDeepHarmonium
    , ConditionalHarmonium
    , ConditionalMixture
    -- ** Construction
    , joinConditionalDeepHarmonium
    , joinConditionalHarmonium
    , splitConditionalDeepHarmonium
    , splitConditionalHarmonium
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import Goal.Graphical.Generative.Harmonium


--- Generic ---



--- Types ---


-- | A generic conditional model.
data ConditionalBias (f :: Type -> Type -> Type) z x

-- | Another generic conditional model.
data ConditionalBiases (f :: Type -> Type -> Type) z x

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium g f z x = ConditionalBias g (Harmonium f z x)

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture f y k = ConditionalHarmonium f y Tensor (Categorical k) -- ^ Function

-- | Splits a conditional 'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalHarmonium
    :: ( Map Natural f y z, Manifold (g y x), Manifold (Harmonium x gxs) )
    => c # ConditionalDeepHarmonium f y ('(g,x) : gxs) z -- ^ Conditional Harmonium
    -> (c # f y z, c # g y x, c # DeepHarmonium x gxs) -- ^ Matrix function and upper part
splitConditionalDeepHarmonium dhrm =
    let (fcs,gdhrmcs) = S.splitAt $ coordinates dhrm
        (gcs,dhrmcs) = S.splitAt gdhrmcs
     in (Point fcs,Point gcs,Point dhrmcs)

splitConditionalHarmonium
    :: ( Map Natural f y z, Manifold (g y x), Manifold x )
    => c # ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
    -> (c # f y z, c # g y x, c # x) -- ^ Matrix function and upper part
splitConditionalHarmonium dhrm =
    let (fyz,gyx,nx0) = splitConditionalDeepHarmonium dhrm
     in (fyz,gyx,fromOneHarmonium nx0)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalDeepHarmonium
    :: ( Map Natural f y z, Manifold (g y x), Manifold (DeepHarmonium x gxs) )
    => c # f y z
    -> c # g y x
    -> c # DeepHarmonium x gxs -- ^ Matrix function and upper part
    -> c # ConditionalDeepHarmonium f y ('(g,x) : gxs) z -- ^ Conditional Harmonium
joinConditionalDeepHarmonium (Point fcs) (Point gcs) (Point dhrmcs) =
    Point $ fcs S.++ gcs S.++ dhrmcs

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalHarmonium
    :: ( Map Natural f y z, Manifold (g y x), Manifold x )
    => c # f y z
    -> c # g y x
    -> c # x -- ^ Matrix function and upper part
    -> c # ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
joinConditionalHarmonium fyz gyx x =
    joinConditionalDeepHarmonium fyz gyx $ toOneHarmonium x


--- Instances ---


instance ( Map Natural f y z, Manifold (g y x), Manifold (DeepHarmonium x gxs) )
  => Manifold (ConditionalDeepHarmonium f y ('(g,x) : gxs) z) where
      type Dimension (ConditionalDeepHarmonium f y ('(g,x) : gxs) z)
        = Dimension (f y z) + Dimension (g y x) + Dimension (DeepHarmonium x gxs)

instance ( Map Natural f y z, Manifold (g y x), Manifold x
         , Manifold (DeepHarmonium x gxs) )
     => Map Natural (ConditionalBias f) (DeepHarmonium y ('(g,x) : gxs)) z where
    (>.>) chrm mz =
        let (fyz,gyx,dhrm) = splitConditionalDeepHarmonium chrm
            affyx = joinAffine (fyz >.> mz) gyx
         in joinBottomHarmonium affyx dhrm
    (>$>) chrm mz =
        let (fyz,gyx,dhrm) = splitConditionalDeepHarmonium chrm
            affyxs = flip joinAffine gyx <$> (fyz >$> mz)
         in flip joinBottomHarmonium dhrm <$> affyxs

instance ( Propagate Natural f y z, Manifold (g y x), Manifold x
         , Manifold (DeepHarmonium x gxs) )
     => Propagate Natural (ConditionalBias f) (DeepHarmonium y ('(g,x) : gxs)) z where
        propagate mdhrms mzs chrm =
            let (mys,mgyxs,mdhrms') = unzip3 $ splitDeepHarmonium <$> mdhrms
                (nfyz,ngyx,ndhrm) = splitConditionalDeepHarmonium chrm
                (mfyz,nyhts) = propagate mys mzs nfyz
             in ( joinConditionalDeepHarmonium mfyz (average mgyxs) (average mdhrms')
                , [joinDeepHarmonium nyht ngyx ndhrm | nyht <- nyhts] )


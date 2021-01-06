{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Hybrid
    (
     -- * Conditional Harmoniums
      ConditionalHarmonium
    , ConditionalMixture
    -- ** Construction
    , joinConditionalHarmonium
    , splitConditionalHarmonium
    ) where


--- Imports  ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Goal.Graphical.Generative.Harmonium

import qualified Goal.Core.Vector.Storable as S


--- Generic ---



--- Types ---


-- | A generic conditional model.
data ConditionalBias (f :: Type -> Type -> Type) z y

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium g f z x = ConditionalBias g (Harmonium f z x)

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture f z k = ConditionalHarmonium f Tensor z (Categorical k) -- ^ Function

-- | Splits a conditional 'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalHarmonium
    :: ( Map Natural f z x, Manifold (g z y), Manifold x )
    => c # ConditionalHarmonium g f z x y -- ^ Conditional Harmonium
    -> (c # g z y, c # f z x, c # x) -- ^ Matrix function and upper part
splitConditionalHarmonium dhrm =
    let (cg,cs') = S.splitAt $ coordinates dhrm
        (cf,cx) = S.splitAt cs'
     in (Point cg,Point cf,Point cx)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalHarmonium
    :: ( Map Natural f z x, Manifold (g z y), Manifold x )
    => c # g z y
    -> c # f z x
    -> c # x
    -> c # ConditionalHarmonium g f z x y -- ^ Conditional Harmonium
joinConditionalHarmonium gzy fzx x =
    Point $ coordinates gzy S.++ coordinates fzx S.++ coordinates x


--- Instances ---


instance ( Map Natural g z y, Manifold (f z x), Manifold x )
  => Manifold (ConditionalHarmonium g f z x y) where
      type Dimension (ConditionalHarmonium g f z x y)
        = Dimension (g z y) + Dimension (f z x) + Dimension x

instance ( Map Natural g z y, Map Natural f z x, Manifold (f z x), Manifold x )
     => Map Natural (ConditionalBias g) (Harmonium f z x) y where
    (>.>) chrm my =
        let (gzy,fzx,nx) = splitConditionalHarmonium chrm
            affzx = joinAffine (gzy >.> my) fzx
         in joinBottomHarmonium affzx nx
    (>$>) chrm mys =
        let (gzy,fzx,nx) = splitConditionalHarmonium chrm
            affzxs = flip joinAffine fzx <$> (gzy >$> mys)
         in flip joinBottomHarmonium nx <$> affzxs

instance ( Propagate Natural g z y, Map Natural f z x, Manifold (f z x), Manifold x )
     => Propagate Natural (ConditionalBias g) (Harmonium f z x) y where
        propagate mhrms mys chrm =
            let (mzs,mfzxs,mxs) = unzip3 $ splitHarmonium <$> mhrms
                (ngzy,nfzx,nx) = splitConditionalHarmonium chrm
                (mgzx,nzhts) = propagate mzs mys ngzy
             in ( joinConditionalHarmonium mgzx (average mfzxs) (average mxs)
                , [joinHarmonium nzht nfzx nx | nzht <- nzhts] )

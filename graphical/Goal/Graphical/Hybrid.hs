{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Graphical.Hybrid
    (
     -- * Conditional Harmoniums
      ConditionalHarmonium
    , ConditionalMixture
    , LatentProcess
    -- ** Construction
    , joinConditionalHarmonium
    , splitConditionalHarmonium
    , joinLatentProcess
    , splitLatentProcess
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

splitLatentProcess
    :: ( Manifold (f x x), Manifold (g z x), Manifold x, Manifold z )
    => c # LatentProcess f g z x
    -> (c # Affine f x x, c # Affine g z x, c # x)
splitLatentProcess ltnt =
    let (cf,cs') = S.splitAt $ coordinates ltnt
        (cg,cx) = S.splitAt cs'
     in (Point cf,Point cg,Point cx)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinLatentProcess
    :: ( Manifold (f x x), Manifold (g z x), Manifold x, Manifold z )
    => c # Affine f x x
    -> c # Affine g z x
    -> c # x
    -> c # LatentProcess f g z x -- ^ Conditional Harmonium
joinLatentProcess cf cg cx =
    Point $ coordinates cf S.++ coordinates cg S.++ coordinates cx



--- Latent Process ---

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
data LatentProcess0 (f :: Type -> Type -> Type) (g :: Type -> Type -> Type) z x

type LatentProcess f g z x = LatentProcess0 f g [z] [x]


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

instance ( Manifold (f x x), Manifold (g z x), Manifold x, Manifold z )
  => Manifold (LatentProcess f g z x) where
      type Dimension (LatentProcess f g z x)
        = Dimension (Affine f x x) + Dimension (Affine g z x) + Dimension x

instance Manifold (LatentProcess f g z x) => Statistical (LatentProcess f g z x) where
    type SamplePoint (LatentProcess f g z x) = [SamplePoint (z,x)]



{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | 'Statistical' models where the observable biases depend on additional inputs.
module Goal.Probability.ExponentialFamily.Harmonium.Conditional
    ( -- * Types
      ConditionalBias
    , ConditionalBiases
    , ConditionalDeepHarmonium
    , ConditionalHarmonium
    , ConditionalMixture
    , ConditionalMixture'
    -- * Construction
    , joinConditionalHarmonium
    , splitConditionalHarmonium
    , joinConditionalDeepHarmonium
    , splitConditionalDeepHarmonium
    -- * Analysis
    , conditionalHarmoniumExpectationStep
    , conditionalHarmoniumExpectationStep'
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
type ConditionalHarmonium f y g x = ConditionalBiases f (Harmonium y (g :: Type -> Type -> Type) x)

-- | A conditional 'Harmonium', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalHarmonium' f y g x = ConditionalDeepHarmonium f y ('[ '(g,x)])

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture f y k = ConditionalHarmonium f y Tensor (Categorical k) -- ^ Function

-- | A conditional 'Mixture', where the observable biases of the
-- 'Harmonium' model depend on additional variables.
type ConditionalMixture' f y k = ConditionalHarmonium' f y Tensor (Categorical k) -- ^ Function

-- | Splits a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' into the
-- unbiased harmonium and the function which models the dependence.
splitConditionalHarmonium
    :: Manifold (Harmonium y g x)
    => c #> ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
    -> (c # Harmonium y g x, c #> f (y,x) z) -- ^ Matrix function and upper part
{-# INLINE splitConditionalHarmonium #-}
splitConditionalHarmonium chrm =
    let (hrmcs,fcs) = S.splitAt $ coordinates chrm
     in (Point hrmcs, Point fcs)

-- | Creates a conditional 'DeepHarmonium'/'Harmonium'/'Mixture' given an
-- unbiased harmonium and a function which models the dependence.
joinConditionalHarmonium
    :: Manifold (Harmonium y g x)
    => c # Harmonium y g x
    -> c #> f (y,x) z
    -> c #> ConditionalHarmonium f y g x z -- ^ Conditional Harmonium
{-# INLINE joinConditionalHarmonium #-}
joinConditionalHarmonium (Point hrmcs) (Point fcs) =
    Point $ hrmcs S.++ fcs

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

-- | Empirical expectations of a conditional harmonium.
conditionalHarmoniumExpectationStep
    :: ( ExponentialFamily y, ExponentialFamily z, Bilinear g y x
       , Map Mean Natural g x y , Map Mean Natural f (y,x) z
       , Manifold (f (y,x) z), LegendreExponentialFamily x )
    => Sample (y,z) -- ^ Model Samples
    -> Natural #> ConditionalHarmonium f y g x z -- ^ Harmonium
    -> [Mean # Harmonium y g x] -- ^ Harmonium expected sufficient statistics
{-# INLINE conditionalHarmoniumExpectationStep #-}
conditionalHarmoniumExpectationStep yzs chrm =
    let (ys,zs) = unzip yzs
        (hrm,fyxz) = splitConditionalHarmonium chrm
        (_,nyx,nx) = splitHarmonium hrm
        nxs0 = snd . splitPair <$> fyxz >$>* zs
        mys = sufficientStatistic <$> ys
        mxs = transition . (+ nx) <$> zipWith (+) nxs0 (mys <$< nyx)
     in zipWith3 joinHarmonium mys (zipWith (>.<) mys mxs) mxs

-- | Empirical expectations of a conditional harmonium.
conditionalHarmoniumExpectationStep'
    :: ( ExponentialFamily y, Bilinear g y x, Map Mean Natural g x y
       , Manifold (f y z), LegendreExponentialFamily x )
    => Sample y -- ^ Model Samples
    -> Natural #> ConditionalHarmonium' f y g x z -- ^ Harmonium
    -> [Mean # Harmonium y g x] -- ^ Harmonium expected sufficient statistics
{-# INLINE conditionalHarmoniumExpectationStep' #-}
conditionalHarmoniumExpectationStep' ys chrm =
    let (hrm,_) = splitConditionalDeepHarmonium chrm
        (_,nyx,nx) = splitHarmonium hrm
        mys = sufficientStatistic <$> ys
        mxs = transition . (+ nx) <$> mys <$< nyx
     in zipWith3 joinHarmonium mys (zipWith (>.<) mys mxs) mxs


--- Instances ---


instance (Map Mean Natural f (y,x) z, Manifold (Harmonium y g x))
    => Manifold (ConditionalHarmonium f y g x z) where
        type Dimension (ConditionalHarmonium f y g x z)
          = Dimension (Harmonium y g x) + Dimension (f (y,x) z)

instance ( Map Mean Natural f (y,x) z, Manifold (g y x)
         , Manifold (Harmonium y g x), Manifold y, Manifold x )
     => Map Mean Natural (ConditionalBiases f) (Harmonium y g x) z where
    {-# INLINE (>.>) #-}
    (>.>) pdhrm mzs =
        let (hrm,fyxz) = splitConditionalHarmonium pdhrm
            (ny,nyx,nx) = splitHarmonium hrm
            (ny',nx') = splitPair $ fyxz >.> mzs
         in joinHarmonium (ny + ny') nyx (nx + nx')
    {-# INLINE (>$>) #-}
    (>$>) pdhrm mzs =
        let (hrm,fyxz) = splitConditionalHarmonium pdhrm
            (ny,nyx,nx) = splitHarmonium hrm
            nyxs = fyxz >$> mzs
         in [ joinHarmonium (ny + ny') nyx (nx + nx') | (ny',nx') <- splitPair <$> nyxs ]

instance ( Propagate Mean Natural f (y,x) z, Manifold (Harmonium y g x)
         , Manifold y, Manifold x, Manifold (g y x) )
  => Propagate Mean Natural (ConditionalBiases f) (Harmonium y g x) z where
        {-# INLINE propagate #-}
        propagate dhrms mzs chrm =
            let (dys,_,dxs) = unzip3 $ splitHarmonium <$> dhrms
                (hrm,f) = splitConditionalHarmonium chrm
                (ny,nyx,nx) = splitHarmonium hrm
                (df,nyxs) = propagate (zipWith joinPair dys dxs) mzs f
             in ( joinConditionalHarmonium (average dhrms) df
                , [ joinHarmonium (ny + ny') nyx (nx + nx') | (ny',nx') <- splitPair <$> nyxs ] )

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


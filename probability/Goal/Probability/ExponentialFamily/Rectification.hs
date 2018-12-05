{-# LANGUAGE
    DataKinds,
    TypeOperators,
    TypeFamilies,
    FlexibleContexts,
    ScopedTypeVariables,
    RankNTypes
    #-}
-- | Population codes and exponential families.
module Goal.Probability.ExponentialFamily.Rectification
    (
    -- * Rectification
      regressRectificationParameters
    , rectificationCurve
    , mixtureLikelihoodRectificationParameters
    -- , rectifiedBayesRule
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

import Goal.Probability.Statistical
import Goal.Probability.Distributions
import Goal.Probability.ExponentialFamily


--- Population Encoders ---


-- | Computes the rectification parameters of a likelihood defined by a categorical latent variable.
mixtureLikelihoodRectificationParameters
    :: (KnownNat k, Enum e, Legendre Natural z)
    => Mean #> Natural # z <* Categorical e k -- ^ Categorical likelihood
    -> (Double, Natural # Categorical e k) -- ^ Rectification parameters
{-# INLINE mixtureLikelihoodRectificationParameters #-}
mixtureLikelihoodRectificationParameters aff =
    let (nz,nzx) = splitAffine aff
        rho0 = potential nz
        rprms = S.map (\nzxi -> subtract rho0 . potential $ nz <+> Point nzxi) $ S.toColumns (toMatrix nzx)
     in (rho0, Point rprms)

-- Population Code Rectification


-- | Computes the rectification curve given a set of rectification parameters,
-- at the given set of points.
rectificationCurve
    :: ExponentialFamily m
    => Double -- ^ Rectification shift
    -> Natural # m -- ^ Rectification parameters
    -> Sample m -- ^ Samples points
    -> [Double] -- ^ Rectification curve at sample points
{-# INLINE rectificationCurve #-}
rectificationCurve rho0 rprms mus = (\x -> rprms <.> sufficientStatistic x + rho0) <$> mus

-- Linear Least Squares

-- | Returns the rectification parameters which best satisfy the rectification
-- equation for the given population code.
regressRectificationParameters
    :: (Map Mean Natural f z x, ExponentialFamily x, Legendre Natural z)
    => Mean #> Natural # f z x -- ^ PPC
    -> Sample x -- ^ Sample points
    -> (Double, Natural # x) -- ^ Approximate rectification parameters
{-# INLINE regressRectificationParameters #-}
regressRectificationParameters lkl mus =
    let dpnds = potential <$> lkl >$>* mus
        indpnds = independentVariables0 lkl mus
        (rho0,rprms) = S.splitAt $ S.linearLeastSquares indpnds dpnds
     in (S.head rho0, Point rprms)

---- | The posterior distribution given a prior and likelihood, where the
---- likelihood is rectified.
--rectifiedBayesRule
--    :: ( Map Mean Natural f z x )
--      => Natural # n -- ^ Rectification Parameters
--      -> Mean #> Natural # f z x -- ^ Likelihood
--      -> SamplePoint z -- ^ Observation
--      -> Natural # x -- ^ Prior
--      -> Natural # x -- ^ Updated prior
--{-# INLINE rectifiedBayesRule #-}
--rectifiedBayesRule rprms lkl x dhrm =
--    let dhrm' = joinBottomHarmonium lkl $ biasBottom ((-1) .> rprms) dhrm
--     in dhrm' lkl >.>* x
--

--- Internal ---

independentVariables0
    :: forall f x z
    . ExponentialFamily x
    => Mean #> Natural # f z x
    -> Sample x
    -> [S.Vector (Dimension x + 1) Double]
independentVariables0 _ mus =
    let sss :: [Mean # x]
        sss = sufficientStatistic <$> mus
     in (S.singleton 1 S.++) . coordinates <$> sss



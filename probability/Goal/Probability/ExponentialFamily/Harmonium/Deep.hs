{-# LANGUAGE UndecidableInstances #-}

-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium.Deep where
--    ( -- * Harmoniums
--      DeepHarmonium (DeepHarmonium)
--    -- ** Structure Manipulation
--    , splitDeepHarmonium
--    , joinDeepHarmonium
--    -- ** Conditional Distributions
--    , conditionalLatentHarmonium
--    -- ** Statistics
--    , deepHarmoniumGradientCalculator
--    , sampleStronglyRectifiedDeepHarmonium
--    , strongRectifierDifferentials
--    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
--import Goal.Probability.ExponentialFamily.Harmonium.Rectified

import qualified Goal.Core.Vector.Storable as S


--- Types ---

-- | An exponential family defined using a product of two other exponential
-- families. The first argument represents the so-called observable variables, and
-- the second argument the so-called latent variables.
data DeepHarmonium m n o

deepHarmoniumBaseMeasure
    :: forall m n o. (ExponentialFamily m, ExponentialFamily n, ExponentialFamily o)
    => Proxy (DeepHarmonium m n o)
    -> Sample (DeepHarmonium m n o)
    -> Double
deepHarmoniumBaseMeasure _ xyz =
    let (x,yz) = S.splitAt xyz
        (y,z) = S.splitAt yz
     in baseMeasure (Proxy :: Proxy m) x * baseMeasure (Proxy :: Proxy n) y * baseMeasure (Proxy :: Proxy o) z


-- Datatype manipulation --

-- | Splits a 'DeepHarmonium' into its components parts of a pair of biases and a 'Product'.
splitDeepHarmonium
    :: (Manifold m, Manifold n, Manifold o)
    => Point c (DeepHarmonium m n o) -- ^ The 'Harmonium'
    -> ( Point c m, Point c n, Point c o
       , Point (Function (Dual c) c) (Product m n), Point (Function (Dual c) c) (Product n o))
splitDeepHarmonium (Point cs) =
    let (mcs,cs') = S.splitAt cs
        (ncs,cs'') = S.splitAt cs'
        (ocs,cs''') = S.splitAt cs''
        (mncs,nocs) = S.splitAt cs'''
     in (Point mcs, Point ncs, Point ocs, Point mncs, Point nocs)

-- | Assembles a 'Harmonium' out of the component parameters.
joinDeepHarmonium
    :: (Manifold m, Manifold n, Manifold o)
    => Point c m
    -> Point c n
    -> Point c o
    -> Point (Function (Dual c) c) (Product m n)
    -> Point (Function (Dual c) c) (Product n o)
    -> Point c (DeepHarmonium m n o) -- ^ The 'Harmonium'
joinDeepHarmonium (Point mcs) (Point ncs) (Point ocs) (Point mncs) (Point nocs) =
    Point $ mcs S.++ ncs S.++ ocs S.++ mncs S.++ nocs

-- | Returns the conditional distribution of the latent variables given the
-- sufficient statistics of the observable state.
--conditionalLatentHarmonium
--    :: (Manifold m, Manifold n, Manifold o, KnownNat k)
--    => Point Natural (DeepHarmonium m n o)
--    -> S.Vector k (Point Mean o)
--    -> S.Vector k (Point Natural (m <*> n))
--conditionalLatentHarmonium dhrm mos =
--    let (nm,nn,_,nmn,nno) = splitDeepHarmonium dhrm
--     in S.map (\nn' -> joinHarmonium nm (nn' <+> nn) nmn) . splitReplicated $ nno >$> mos

--sampleStronglyRectifiedDeepHarmonium
--    :: ( ExponentialFamily x, SourceGenerative Natural x
--       , ExponentialFamily y, SourceGenerative Natural y, SourceGenerative Natural z)
--    => Int
--    -> Natural :#: x
--    -> Natural :#: y
--    -> Natural :#: DeepHarmonium x y z
--    -> RandST s [(Sample x, Sample y, Sample z)]
--sampleStronglyRectifiedDeepHarmonium n rx ry dhrm = do
--    xs <- replicateM n $ standardGenerate rx
--    let (_,_,nz,nxy,nyz) = splitDeepHarmonium dhrm
--    ys <- mapM standardGenerate $ joinAffine ry (matrixTranspose nxy) >$>* xs
--    zs <- mapM standardGenerate $ joinAffine nz (matrixTranspose nyz) >$>* ys
--    return $ zip3 xs ys zs
--
--strongRectifierDifferentials
--    :: ( ClosedFormExponentialFamily x, ClosedFormExponentialFamily y, ClosedFormExponentialFamily z
--       , SourceGenerative Natural x, SourceGenerative Natural y )
--    => Int
--    -> Natural :#: x
--    -> Natural :#: y
--    -> Natural :#: DeepHarmonium x y z
--    -> RandST s (Differentials :#: Tangent Natural x,Differentials :#: Tangent Natural y)
--strongRectifierDifferentials n rx ry dhrm = do
--    let (nx,ny,nz,nxy,nyz) = splitDeepHarmonium dhrm
--        hrm1 = joinHarmonium ny nz nyz
--        hrm2 = joinHarmonium nx ry nxy
--    drx <- rectifierDifferentials n rx hrm2
--    dry <- rectifierDifferentials n ry hrm1
--    return (drx,dry)
--
--deepHarmoniumGradientCalculator
--    :: [Mixture :#: x]
--    -> [Mixture :#: x]
--    -> [Mixture :#: y]
--    -> [Mixture :#: y]
--    -> [Mixture :#: z]
--    -> [Mixture :#: z]
--    -> Natural :#: DeepHarmonium x y z
--    -> Differentials :#: Tangent Natural (DeepHarmonium x y z)
--deepHarmoniumGradientCalculator mxs mxs' mys mys' mzs mzs' dhrm =
--    let dx = averagePoint $ zipWith (<->) mxs mxs'
--        dy = averagePoint $ zipWith (<->) mys mys'
--        dz = averagePoint $ zipWith (<->) mzs mzs'
--        dxy = averagePoint [ (mx >.< my) <-> (mx' >.< my') | (mx,mx',my,my') <- zip4 mxs mxs' mys mys' ]
--        dyz = averagePoint [ (my >.< mz) <-> (my' >.< mz') | (my,my',mz,mz') <- zip4 mys mys' mzs mzs' ]
--     in fromCoordinates (Tangent dhrm)
--        $ coordinates dx C.++ coordinates dy C.++ coordinates dz C.++ coordinates dxy C.++ coordinates dyz
--
--
----- Instances ---
--
--
---- Harmoniums --

instance (Manifold m, Manifold n, Manifold o) => Manifold (DeepHarmonium m n o) where
    type Dimension (DeepHarmonium m n o) = Dimension m + Dimension n + Dimension o + Dimension m * Dimension n + Dimension n * Dimension o

instance (Statistical m, Statistical n, Statistical o) => Statistical (DeepHarmonium m n o) where
    type SampleDimension (DeepHarmonium m n o) = SampleDimension m + SampleDimension n + SampleDimension o

instance (ExponentialFamily m, ExponentialFamily n, ExponentialFamily o)
  => ExponentialFamily (DeepHarmonium m n o) where
      sufficientStatistic xmno =
          let (xm,xno) = S.splitAt xmno
              (xn,xo) = S.splitAt xno
              sm = sufficientStatistic xm
              sn = sufficientStatistic xn
              so = sufficientStatistic xo
           in joinDeepHarmonium sm sn so (sm >.< sn) (sn >.< so)
      baseMeasure = deepHarmoniumBaseMeasure

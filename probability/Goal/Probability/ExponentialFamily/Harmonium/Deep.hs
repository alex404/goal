-- | Exponential Family Harmoniums. Gibbs sampling is defined in 'goal-simulation'.
module Goal.Probability.ExponentialFamily.Harmonium.Deep
    ( -- * Harmoniums
      DeepHarmonium (DeepHarmonium)
    -- ** Structure Manipulation
    , splitDeepHarmonium
    , joinDeepHarmonium
    -- ** Conditional Distributions
    , conditionalLatentHarmonium
    -- ** Statistics
    , deepHarmoniumGradientCalculator
    , sampleStronglyRectifiedDeepHarmonium
    , strongRectifierDifferentials
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.ExponentialFamily.Harmonium.Rectified

import System.Random.MWC.Monad hiding (save)


--- Types ---

-- | An exponential family defined using a product of two other exponential
-- families. The first argument represents the so-called observable variables, and
-- the second argument the so-called latent variables.
data DeepHarmonium x y z = DeepHarmonium x y z deriving (Eq, Read, Show)

-- Datatype manipulation --

-- | Splits a 'DeepHarmonium' into its components parts of a pair of biases and a 'Tensor'.
splitDeepHarmonium
    :: (Manifold x, Manifold y, Manifold z)
    => (Natural :#: DeepHarmonium x y z) -- ^ The 'Harmonium'
    -> (Natural :#: x, Natural :#: y, Natural :#: z, NaturalFunction :#: Tensor x y, NaturalFunction :#: Tensor y z) -- ^ The component parameters
splitDeepHarmonium dhrm =
    let (DeepHarmonium x y z) = manifold dhrm
        (xcs,css') = C.splitAt (dimension x) $ coordinates dhrm
        (ycs,css'') = C.splitAt (dimension y) css'
        (zcs,css''') = C.splitAt (dimension z) css''
        (xycs,yzcs) = C.splitAt (dimension $ Tensor x y) css'''
     in ( fromCoordinates x xcs, fromCoordinates y ycs, fromCoordinates z zcs
        , fromCoordinates (Tensor x y) xycs, fromCoordinates (Tensor y z) yzcs )

-- | Assembles a 'Harmonium' out of the component parameters.
joinDeepHarmonium
    :: (Natural :#: x) -- ^ Latent biases
    -> (Natural :#: y) -- ^ Observable biases
    -> (Natural :#: z) -- ^ Observable biases
    -> (NaturalFunction :#: Tensor x y) -- ^ Interaction matrix
    -> (NaturalFunction :#: Tensor y z) -- ^ Interaction matrix
    -> (Natural :#: DeepHarmonium x y z) -- ^ The resulting 'Harmonium'
joinDeepHarmonium nx ny nz nxy nxz =
    fromCoordinates (DeepHarmonium (manifold nx) (manifold ny) (manifold nz))
        $ coordinates nx C.++ coordinates ny C.++ coordinates nz C.++ coordinates nxy C.++ coordinates nxz

-- | Returns the conditional distribution of the latent variables given the
-- sufficient statistics of the observable state.
conditionalLatentHarmonium
    :: (Manifold x, Manifold y, Manifold z)
    => Natural :#: DeepHarmonium x y z -> [Mixture :#: z] -> [Natural :#: Harmonium x y]
conditionalLatentHarmonium dhrm mzs =
    let (nx,ny,_,nxy,nyz) = splitDeepHarmonium dhrm
     in [ joinHarmonium nx (ny <+> nyz >.> mz) nxy | mz <- mzs ]

sampleStronglyRectifiedDeepHarmonium
    :: ( ExponentialFamily x, SourceGenerative Natural x
       , ExponentialFamily y, SourceGenerative Natural y, SourceGenerative Natural z)
    => Int
    -> Natural :#: x
    -> Natural :#: y
    -> Natural :#: DeepHarmonium x y z
    -> RandST s [(Sample x, Sample y, Sample z)]
sampleStronglyRectifiedDeepHarmonium n rx ry dhrm = do
    xs <- replicateM n $ standardGenerate rx
    let (_,_,nz,nxy,nyz) = splitDeepHarmonium dhrm
    ys <- mapM standardGenerate $ joinAffine ry (matrixTranspose nxy) >$>* xs
    zs <- mapM standardGenerate $ joinAffine nz (matrixTranspose nyz) >$>* ys
    return $ zip3 xs ys zs

strongRectifierDifferentials
    :: ( ClosedFormExponentialFamily x, ClosedFormExponentialFamily y, ClosedFormExponentialFamily z
       , SourceGenerative Natural x, SourceGenerative Natural y )
    => Int
    -> Natural :#: x
    -> Natural :#: y
    -> Natural :#: DeepHarmonium x y z
    -> RandST s (Differentials :#: Tangent Natural x,Differentials :#: Tangent Natural y)
strongRectifierDifferentials n rx ry dhrm = do
    let (nx,ny,nz,nxy,nyz) = splitDeepHarmonium dhrm
        hrm1 = joinHarmonium ny nz nyz
        hrm2 = joinHarmonium nx ry nxy
    drx <- rectifierDifferentials n rx hrm2
    dry <- rectifierDifferentials n ry hrm1
    return (drx,dry)

deepHarmoniumGradientCalculator
    :: [Mixture :#: x]
    -> [Mixture :#: x]
    -> [Mixture :#: y]
    -> [Mixture :#: y]
    -> [Mixture :#: z]
    -> [Mixture :#: z]
    -> Natural :#: DeepHarmonium x y z
    -> Differentials :#: Tangent Natural (DeepHarmonium x y z)
deepHarmoniumGradientCalculator mxs mxs' mys mys' mzs mzs' dhrm =
    let dx = averagePoint $ zipWith (<->) mxs mxs'
        dy = averagePoint $ zipWith (<->) mys mys'
        dz = averagePoint $ zipWith (<->) mzs mzs'
        dxy = averagePoint [ (mx >.< my) <-> (mx' >.< my') | (mx,mx',my,my') <- zip4 mxs mxs' mys mys' ]
        dyz = averagePoint [ (my >.< mz) <-> (my' >.< mz') | (my,my',mz,mz') <- zip4 mys mys' mzs mzs' ]
     in fromCoordinates (Tangent dhrm)
        $ coordinates dx C.++ coordinates dy C.++ coordinates dz C.++ coordinates dxy C.++ coordinates dyz


--- Instances ---


-- Harmoniums --

instance (Manifold x, Manifold y, Manifold z) => Manifold (DeepHarmonium x y z) where
    dimension (DeepHarmonium x y z) = dimension x + dimension y + dimension z + dimension x * dimension y + dimension y * dimension z

instance (Chart Natural x, Chart Natural y, Chart Natural z) => Chart Natural (DeepHarmonium x y z) where

{-
instance (Statistical l, Statistical o) => Statistical (Harmonium l o) where
    type SampleSpace (Harmonium l o) = (SampleSpace l, SampleSpace o)
    sampleSpace (Harmonium l o) = (sampleSpace l, sampleSpace o)

instance (ExponentialFamily l, ExponentialFamily o) => ExponentialFamily (Harmonium l o) where
    sufficientStatistic (Harmonium l o) (xl, xo) =
        let slcs = coordinates $ sufficientStatistic l xl
            socs = coordinates $ sufficientStatistic o xo
            smtxcs = C.concat [C.map (slc*) socs | slc <- C.toList slcs]
         in fromCoordinates (Harmonium l o) $ C.concat [slcs,socs,smtxcs]
    baseMeasure (Harmonium l o) (xl, xo) = baseMeasure l xl * baseMeasure o xo
    -}

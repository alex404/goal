module Goal.Probability.ExponentialFamily.Harmonium.Rectification where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.Distributions

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G


--- Functions ---


categoricalHarmoniumRectificationParameters
    :: (Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily z)
    => Point Natural (Categorical e k <*> z)
    -> (Double, Point Natural (Categorical e k))
categoricalHarmoniumRectificationParameters hrm =
    let (_,nz,nxz) = splitHarmonium hrm
        rho0 = potential nz
        rprms = S.map (\nxzi -> subtract rho0 . potential $ nz <+> Point nxzi) $ S.toRows (toMatrix nxz)
     in (rho0, Point rprms)

sampleRectifiedHarmonium
    :: ( Bilinear Mean Natural f, ExponentialFamily (Codomain f)
       , Generative Natural (Codomain f), Generative Natural (Domain f)
       , KnownNat k )
      => Point Natural (Codomain f)
      -> Point Natural (Harmonium f)
      -> Random s (B.Vector k (Sample (Harmonium f)))
sampleRectifiedHarmonium rprms hrm = do
    let (lx,_,_) = splitHarmonium hrm
    xs <- B.replicateM $ sample (lx <+> rprms)
    zs <- B.mapM sample . G.convert $ conditionalObservableDistributions hrm xs
    return $ B.zip xs zs

sampleCategoricalHarmonium
    :: ( Enum e, KnownNat k
       , KnownNat n, 1 <= n, ClosedFormExponentialFamily o, Generative Natural o )
      => Point Natural (Categorical e n <*> o)
    -> Random s (B.Vector k (Sample (Categorical e n <*> o)))
sampleCategoricalHarmonium hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    sampleRectifiedHarmonium rx hrm

estimateRectifiedHarmoniumDifferentials
    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f), ExponentialFamily (Codomain f)
       , ExponentialFamily (Domain f), Generative Natural (Codomain f)
       , Generative Natural (Domain f), KnownNat k )
    => B.Vector k (Sample (Domain f))
    -> Point Natural (Codomain f) -- ^ Rectification Parameters
    -> Point Natural (Harmonium f)
    -> Random s (CotangentVector Natural (Harmonium f))
estimateRectifiedHarmoniumDifferentials pzs rprms hrm = do
    pxs <- B.mapM sample . G.convert $ conditionalLatentDistributions hrm pzs
    let pxzs = B.zip pxs pzs
    qxzs <- sampleRectifiedHarmonium rprms hrm
    return $ estimateStochasticCrossEntropyDifferential pxzs qxzs

estimateCategoricalHarmoniumDifferentials
    :: ( Enum e, KnownNat k, 1 <= k
       , ClosedFormExponentialFamily o, Generative Natural o, KnownNat n )
    => B.Vector n (Sample o)
    -> Point Natural (Categorical e k <*> o)
    -> Random s (CotangentVector Natural (Categorical e k <*> o))
estimateCategoricalHarmoniumDifferentials pos hrm = do
    pls <- B.mapM sample . G.convert $ conditionalLatentDistributions hrm pos
    let plos = B.zip (G.convert pls) (G.convert pos)
    qlos <- sampleCategoricalHarmonium hrm
    return $ estimateStochasticCrossEntropyDifferential plos qlos

rectifiedHarmoniumNegativeLogLikelihood
    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f)
       , ClosedFormExponentialFamily (Codomain f), ClosedFormExponentialFamily (Domain f) )
      => (Double, Point Natural (Codomain f)) -- ^ Rectification Parameters
      -> Point Natural (Harmonium f)
      -> Sample (Domain f)
      -> Double
rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
    let (nl,no,nlo) = splitHarmonium hrm
     in negate $ sufficientStatistic ox <.> no + potential (nl <+> nlo >.>* ox) - potential (nl <+> rprms) - rho0

categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o )
    => Point Natural (Categorical e k <*> o)
    -> Sample o
    -> Double
categoricalHarmoniumNegativeLogLikelihood hrm =
    rectifiedHarmoniumNegativeLogLikelihood (categoricalHarmoniumRectificationParameters hrm) hrm


{-
rectifiedHarmoniumRelativeEntropyDifferentials
    :: (ExponentialFamily x, SourceGenerative Natural x, ClosedFormExponentialFamily z, SourceGenerative Natural z)
    => [Sample z]
    -> Natural :#: x
    -> Natural :#: Harmonium x z
    -> RandST s (Differentials :#: Tangent Natural (Harmonium x z))
rectifiedHarmoniumRelativeEntropyDifferentials zs' rx hrm = do

    xzs <- sampleRectifiedHarmonium (length zs') rx hrm
    let Harmonium x z = manifold hrm
    let (xs,zs) = unzip xzs
    let mzs' = sufficientStatistic z <$> zs'

    xs' <- mapM standardGenerate $ conditionalLatentDistribution hrm >$> mzs'

    let mxs = sufficientStatistic x <$> xs
        mxs' = sufficientStatistic x <$> xs'
        mzs = sufficientStatistic z <$> zs

    return $ harmoniumGradientCalculator mxs mxs' mzs mzs' hrm

rectifierDifferentials
    :: (ClosedFormExponentialFamily x, ClosedFormExponentialFamily z, SourceGenerative Natural x)
    => Int
    -> Natural :#: x
    -> Natural :#: Harmonium x z
    -> RandST s (Differentials :#: Tangent Natural x)
rectifierDifferentials n rx hrm = do
    let (Harmonium xm _) = manifold hrm
        (nx,nz,imtx) = splitHarmonium hrm
        ssx = sufficientStatistic xm
        covariate x = ssx x <.> (rx <-> nx) - potential (nz <+> matrixTranspose imtx >.> ssx x)
    xs <- replicateM n $ standardGenerate rx
    let cv0 = averagePoint [covariate x .> ssx x | x <- xs ]
        cv1 = average [covariate x | x <- xs ] .> averagePoint [ ssx x | x <- xs ]
        cv = cv0 <-> cv1
    return . fromCoordinates (Tangent rx) $ coordinates cv


-- | Bayes' rule given a rectified harmonium generative model.
rectifiedBayesRule
    :: (ExponentialFamily x, ExponentialFamily z)
    => Function Mixture Natural :#: Affine z x -- ^ Likelihood
    -> Natural :#: x
    -> Sample z -- ^ Observation
    -> Natural :#: x -- ^ Prior
    -> Natural :#: x -- ^ Posterior
rectifiedBayesRule lklhd rprms z p0 =
    let mtx = matrixTranspose . snd $ splitAffine lklhd
     in mtx >.>* z <+> p0 <+> rprms
     -}


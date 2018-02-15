module Goal.Probability.ExponentialFamily.Harmonium.Rectification where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry

import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.ExponentialFamily.Harmonium
import Goal.Probability.Distributions


--- Functions ---


categoricalHarmoniumRectificationParameters
    :: (Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily z, RealFloat x)
    => Point Natural (Categorical e k <*> z) x
    -> (x,Point Natural (Categorical e k) x)
categoricalHarmoniumRectificationParameters hrm =
    let (_,nz,nxz) = splitHarmonium hrm
        rho0 = potential nz
        rprms = (\nxzi -> subtract rho0 . potential $ nz <+> Point nxzi) <$> toRows (toMatrix nxz)
     in (rho0, Point rprms)

sampleRectifiedHarmonium
    :: ( Bilinear Mean Natural f, ExponentialFamily (Codomain f)
       , Generative Natural (Codomain f), Generative Natural (Domain f)
       , KnownNat k, RealFloat x)
    => Point Natural (Codomain f) x
    -> Point Natural (Harmonium f) x
    -> Random s (Vector k (Sample (Harmonium f)))
sampleRectifiedHarmonium rprms hrm = do
    let (lx,_,_) = splitHarmonium hrm
    xs <- replicateMV $ generate (lx <+> rprms)
    zs <- mapM generate $ conditionalObservableDistributions hrm xs
    return $ zipV xs zs

sampleCategoricalHarmonium
    :: ( Enum e, KnownNat n, KnownNat k, 1 <= k, ClosedFormExponentialFamily o, RealFloat x, Generative Natural o )
    => Point Natural (Categorical e k <*> o) x
    -> Random s (Vector n (Sample (Categorical e k <*> o)))
sampleCategoricalHarmonium hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
    sampleRectifiedHarmonium rx hrm

estimateRectifiedHarmoniumDifferentials
    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f), ExponentialFamily (Codomain f), ExponentialFamily (Domain f)
       , Generative Natural (Codomain f), Generative Natural (Domain f), RealFloat x, KnownNat k )
    => Vector k (Sample (Domain f))
    -> Point Natural (Codomain f) x -- ^ Rectification Parameters
    -> Point Natural (Harmonium f) x
    -> Random s (CotangentVector Natural (Harmonium f) x)
estimateRectifiedHarmoniumDifferentials pzs rprms hrm = do
    pxs <- mapM generate $ conditionalLatentDistributions hrm pzs
    let pxzs = zipV pxs pzs
    qxzs <- sampleRectifiedHarmonium rprms hrm
    return $ estimateStochasticCrossEntropyDifferential pxzs qxzs

estimateCategoricalHarmoniumDifferentials
    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o, Generative Natural o, KnownNat n, RealFloat x )
    => Vector n (Sample o)
    -> Point Natural (Categorical e k <*> o) x
    -> Random s (CotangentVector Natural (Categorical e k <*> o) x)
estimateCategoricalHarmoniumDifferentials pos hrm = do
    pls <- mapM generate $ conditionalLatentDistributions hrm pos
    let plos = zipV pls pos
    qlos <- sampleCategoricalHarmonium hrm
    return $ estimateStochasticCrossEntropyDifferential plos qlos

rectifiedHarmoniumNegativeLogLikelihood
    :: ( Bilinear Mean Natural f, ExponentialFamily (Harmonium f)
       , ClosedFormExponentialFamily (Codomain f), ClosedFormExponentialFamily (Domain f), RealFloat x )
    => (x, Point Natural (Codomain f) x) -- ^ Rectification Parameters
    -> Point Natural (Harmonium f) x
    -> Sample (Domain f)
    -> x
rectifiedHarmoniumNegativeLogLikelihood (rho0,rprms) hrm ox =
    let (nl,no,nlo) = splitHarmonium hrm
     in negate $ sufficientStatistic ox <.> no + potential (nl <+> nlo >.>* ox) - potential (nl <+> rprms) - rho0

categoricalHarmoniumNegativeLogLikelihood
    :: ( Enum e, KnownNat k, 1 <= k, ClosedFormExponentialFamily o, RealFloat x )
    => Point Natural (Categorical e k <*> o) x
    -> Sample o
    -> x
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


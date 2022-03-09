{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances,TypeApplications #-}

-- | Various instances of statistical manifolds, with a focus on exponential
-- families. In the documentation we use \(X\) to indicate a random variable
-- with the distribution being documented.
module Goal.Probability.Distributions.Gaussian
    ( -- * Univariate
      Normal
    , NormalMean
    , NormalVariance
    -- * Multivariate
    , MVNMean
    , MVNCovariance
    , MultivariateNormal
    , IsotropicNormal
    , multivariateNormalCorrelations
    , bivariateNormalConfidenceEllipse
    , splitMultivariateNormal
    , splitMeanMultivariateNormal
    , splitNaturalMultivariateNormal
    , joinMultivariateNormal
    , joinMeanMultivariateNormal
    , joinNaturalMultivariateNormal
    , isotropicNormalToFull
    , fullNormalToIsotropic
    -- * Linear Models
    , SimpleLinearModel
    , LinearModel
    , FactorAnalysis
    , PrincipleComponentAnalysis
    ) where

-- Package --

import Goal.Core
import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified System.Random.MWC.Distributions as R

-- Normal Distribution --

-- | The Mean of a normal distribution. When used as a distribution itself, it
-- is a Normal distribution with unit variance.
data NormalMean

-- | The variance of a normal distribution.
data NormalVariance

-- | The 'Manifold' of 'Normal' distributions. The 'Source' coordinates are the
-- mean and the variance.
type Normal = LocationShape NormalMean NormalVariance

-- | The Mean of a normal distribution. When used as a distribution itself, it
-- is a Normal distribution with unit variance.
data MVNMean (n :: Nat)

-- | The variance of a normal distribution.
data MVNCovariance (n :: Nat)

-- | Linear models are linear functions with additive Guassian noise.
type LinearModel n k = Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)

-- | Linear models are linear functions with additive Guassian noise.
type SimpleLinearModel = Affine Tensor NormalMean Normal NormalMean

type FactorAnalysis n k = Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)
type PrincipleComponentAnalysis n k = Affine Tensor (MVNMean n) (IsotropicNormal n) (MVNMean k)


-- Multivariate Normal --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The 'Source'
-- coordinates are the (vector) mean and the covariance matrix. For the
-- coordinates of a multivariate normal distribution, the elements of the mean
-- come first, and then the elements of the covariance matrix in row major
-- order.
--
-- Note that we only store the lower triangular elements of the covariance
-- matrix, to better reflect the true dimension of a MultivariateNormal
-- Manifold. In short, be careful when using 'join' and 'split' to access the
-- values of the Covariance matrix, and consider using the specific instances
-- for MVNs.
type MultivariateNormal (n :: Nat) = LocationShape (MVNMean n) (MVNCovariance n)

-- | Split a MultivariateNormal into its Means and Covariance matrix.
splitMultivariateNormal
    :: KnownNat n
    => Source # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitMultivariateNormal mvn =
    let (mu,cvr) = split mvn
     in (coordinates mu, S.fromLowerTriangular $ coordinates cvr)

-- | Join a covariance matrix into a MultivariateNormal.
joinMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Source # MultivariateNormal n
joinMultivariateNormal mus sgma =
    join (Point mus) (Point $ S.lowerTriangular sgma)

-- | Split a MultivariateNormal into its Means and Covariance matrix.
splitMeanMultivariateNormal
    :: KnownNat n
    => Mean # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitMeanMultivariateNormal mvn =
    let (mu,cvr) = split mvn
     in (coordinates mu, S.fromLowerTriangular $ coordinates cvr)

-- | Join a covariance matrix into a MultivariateNormal.
joinMeanMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Mean # MultivariateNormal n
joinMeanMultivariateNormal mus sgma =
    join (Point mus) (Point $ S.lowerTriangular sgma)

-- | Split a MultivariateNormal into the precision weighted means and (-0.5*)
-- Precision matrix. Note that this performs an easy to miss computation for
-- converting the natural parameters in our reduced representation of MVNs into
-- the full precision matrix.
splitNaturalMultivariateNormal
    :: KnownNat n
    => Natural # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitNaturalMultivariateNormal np =
    let (nmu,cvrs) = split np
        nmu0 = coordinates nmu
        nsgma0' = (/2) . S.fromLowerTriangular $ coordinates cvrs
        nsgma0 = nsgma0' + S.diagonalMatrix (S.takeDiagonal nsgma0')
     in (nmu0, nsgma0)

-- | Joins a MultivariateNormal out of the precision weighted means and (-0.5)
-- Precision matrix. Note that this performs an easy to miss computation for
-- converting the full precision Matrix into the reduced, EF representation we use here.
joinNaturalMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Natural # MultivariateNormal n
joinNaturalMultivariateNormal nmu0 nsgma0 =
    let nmu = Point nmu0
        diag = S.diagonalMatrix $ S.takeDiagonal nsgma0
     in join nmu . Point . S.lowerTriangular $ 2*nsgma0 - diag

-- | Confidence elipses for bivariate normal distributions.
bivariateNormalConfidenceEllipse
    :: Int
    -> Double
    -> Source # MultivariateNormal 2
    -> [(Double,Double)]
bivariateNormalConfidenceEllipse nstps prcnt nrm =
    let (mu,cvr) = splitMultivariateNormal nrm
        chl = S.withMatrix (S.scale prcnt) $ S.unsafeCholesky cvr
        xs = range 0 (2*pi) nstps
        sxs = [ S.fromTuple (cos x, sin x) | x <- xs ]
     in S.toPair . (mu +) <$> S.matrixMap chl sxs

-- | Computes the correlation matrix of a 'MultivariateNormal' distribution.
multivariateNormalCorrelations
    :: KnownNat k
    => Source # MultivariateNormal k
    -> S.Matrix k k Double
multivariateNormalCorrelations mnrm =
    let cvrs = snd $ splitMultivariateNormal mnrm
        sds = S.map sqrt $ S.takeDiagonal cvrs
        sdmtx = S.outerProduct sds sds
     in G.Matrix $ S.zipWith (/) (G.toVector cvrs) (G.toVector sdmtx)

multivariateNormalLogBaseMeasure
    :: forall n . (KnownNat n)
    => Proxy (MultivariateNormal n)
    -> S.Vector n Double
    -> Double
multivariateNormalLogBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n/2 * log (2*pi)

isotropicNormalLogBaseMeasure
    :: forall n . (KnownNat n)
    => Proxy (IsotropicNormal n)
    -> S.Vector n Double
    -> Double
isotropicNormalLogBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n/2 * log (2*pi)

diagonalNormalLogBaseMeasure
    :: forall n . (KnownNat n)
    => Proxy (DiagonalNormal n)
    -> S.Vector n Double
    -> Double
diagonalNormalLogBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n/2 * log (2*pi)

mvnMeanLogBaseMeasure
    :: forall n . (KnownNat n)
    => Proxy (MVNMean n)
    -> S.Vector n Double
    -> Double
mvnMeanLogBaseMeasure _ x =
    let n = natValInt (Proxy :: Proxy n)
     in -fromIntegral n/2 * log pi - S.dotProduct x x / 2

-- | samples a multivariateNormal by way of a covariance matrix i.e. by taking
-- the square root.
sampleMultivariateNormal
    :: KnownNat n
    => Int
    -> Source # MultivariateNormal n
    -> Random [S.Vector n Double]
sampleMultivariateNormal n p = do
    let (mus,sgma) = splitMultivariateNormal p
    nrms <- replicateM n . S.replicateM $ Random (R.normal 0 1)
    let rtsgma = S.matrixRoot sgma
    return $ (mus +) <$> S.matrixMap rtsgma nrms

isotropicNormalToFull
    :: KnownNat n
    => Natural # IsotropicNormal n
    -> Natural # MultivariateNormal n
isotropicNormalToFull iso =
    let (mus,sgma0) = split iso
        sgma = realToFrac . S.head $ coordinates sgma0
     in joinNaturalMultivariateNormal (coordinates mus) $ sgma * S.matrixIdentity

fullNormalToIsotropic
    :: KnownNat n
    => Mean # MultivariateNormal n
    -> Mean # IsotropicNormal n
fullNormalToIsotropic iso =
    let (mus,sgma0) = splitMeanMultivariateNormal iso
        sgma = S.sum $ S.takeDiagonal sgma0
     in join (Point mus) $ singleton sgma

-- Restricted MVNs --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The 'Source'
-- coordinates are the (vector) mean and the covariance matrix. For the
-- coordinates of a multivariate normal distribution, the elements of the mean
-- come first, and then the elements of the covariance matrix in row major
-- order.
type IsotropicNormal (n :: Nat) = LocationShape (MVNMean n) NormalVariance

type DiagonalNormal (n :: Nat) = LocationShape (MVNMean n) (Replicated n NormalVariance)


--- Internal ---


--- Instances ---


-- NormalMean Distribution --

instance Manifold NormalMean where
    type Dimension NormalMean = 1

instance Statistical NormalMean where
    type SamplePoint NormalMean = Double

instance ExponentialFamily NormalMean where
    sufficientStatistic x = singleton x
    logBaseMeasure _ x = -square x/2 - sqrt (2*pi)

type instance PotentialCoordinates NormalMean = Natural

instance Transition Mean Natural NormalMean where
    transition = breakPoint

instance Transition Mean Source NormalMean where
    transition = breakPoint

instance Transition Source Natural NormalMean where
    transition = breakPoint

instance Transition Source Mean NormalMean where
    transition = breakPoint

instance Transition Natural Mean NormalMean where
    transition = breakPoint

instance Transition Natural Source NormalMean where
    transition = breakPoint

instance Legendre NormalMean where
    potential (Point cs) =
        let tht = S.head cs
         in square tht / 2

instance LogLikelihood Natural NormalMean Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- Normal Shape --


instance Manifold NormalVariance where
    type Dimension NormalVariance = 1


-- Normal Distribution --

instance ExponentialFamily Normal where
    sufficientStatistic x =
         Point . S.doubleton x $ x**2
    logBaseMeasure _ _ = -1/2 * log (2 * pi)

type instance PotentialCoordinates Normal = Natural

instance Legendre Normal where
    potential (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in -(square tht0 / (4*tht1)) - 0.5 * log(-2*tht1)

instance Transition Natural Mean Normal where
    transition p =
        let (tht0,tht1) = S.toPair $ coordinates p
            dv = tht0/tht1
         in Point $ S.doubleton (-0.5*dv) (0.25 * square dv - 0.5/tht1)

instance DuallyFlat Normal where
    dualPotential (Point cs) =
        let (eta0,eta1) = S.toPair cs
         in -0.5 * log(eta1 - square eta0) - 1/2

instance Transition Mean Natural Normal where
    transition p =
        let (eta0,eta1) = S.toPair $ coordinates p
            dff = eta1 - square eta0
         in Point $ S.doubleton (eta0 / dff) (-0.5 / dff)

instance Riemannian Natural Normal where
    metric p =
        let (tht0,tht1) = S.toPair $ coordinates p
            d00 = -1/(2*tht1)
            d01 = tht0/(2*square tht1)
            d11 = 0.5*(1/square tht1 - square tht0 / (tht1^(3 :: Int)))
         in Point $ S.doubleton d00 d01 S.++ S.doubleton d01 d11

instance Riemannian Mean Normal where
    metric p =
        let (eta0,eta1) = S.toPair $ coordinates p
            eta02 = square eta0
            dff2 = square $ eta1 - eta02
            d00 = (dff2 + 2 * eta02) / dff2
            d01 = -eta0 / dff2
            d11 = 0.5 / dff2
         in Point $ S.doubleton d00 d01 S.++ S.doubleton d01 d11

-- instance Riemannian Source Normal where
--     metric p =
--         let (_,vr) = S.toPair $ coordinates p
--          in Point $ S.doubleton (recip vr) 0 S.++ S.doubleton 0 (recip $ 2*square vr)

instance Transition Source Mean Normal where
    transition (Point cs) =
        let (mu,vr) = S.toPair cs
         in Point . S.doubleton mu $ vr + square mu

instance Transition Mean Source Normal where
    transition (Point cs) =
        let (eta0,eta1) = S.toPair cs
         in Point . S.doubleton eta0 $ eta1 - square eta0

instance Transition Source Natural Normal where
    transition (Point cs) =
        let (mu,vr) = S.toPair cs
         in Point $ S.doubleton (mu / vr) (negate . recip $ 2 * vr)

instance Transition Natural Source Normal where
    transition (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in Point $ S.doubleton (-0.5 * tht0 / tht1) (negate . recip $ 2 * tht1)

instance (Transition c Source Normal) => Generative c Normal where
    samplePoint p =
        let (Point cs) = toSource p
            (mu,vr) = S.toPair cs
         in Random $ R.normal mu (sqrt vr)

instance AbsolutelyContinuous Source Normal where
    densities (Point cs) xs = do
        let (mu,vr) = S.toPair cs
        x <- xs
        return $ recip (sqrt $ vr*2*pi) * exp (negate $ (x - mu) ** 2 / (2*vr))

instance AbsolutelyContinuous Mean Normal where
    densities = densities . toSource

instance AbsolutelyContinuous Natural Normal where
    logDensities = exponentialFamilyLogDensities

instance Transition Mean c Normal => MaximumLikelihood c Normal where
    mle = transition . averageSufficientStatistic

instance LogLikelihood Natural Normal Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- MVNMean --

instance KnownNat n => Manifold (MVNMean n) where
    type Dimension (MVNMean n) = n

instance (KnownNat n) => Statistical (MVNMean n) where
    type SamplePoint (MVNMean n) = S.Vector n Double

instance KnownNat n => ExponentialFamily (MVNMean n) where
    sufficientStatistic x = Point x
    logBaseMeasure = mvnMeanLogBaseMeasure

type instance PotentialCoordinates (MVNMean n) = Natural

-- MVNCovariance --

instance (KnownNat n, KnownNat (Triangular n)) => Manifold (MVNCovariance n) where
    type Dimension (MVNCovariance n) = Triangular n

-- Isotropic Normal

instance KnownNat n => ExponentialFamily (IsotropicNormal n) where
    sufficientStatistic xs = Point $ xs S.++ S.singleton (S.dotProduct xs xs)
    averageSufficientStatistic xs = Point $ average xs S.++ S.singleton ( average $ zipWith S.dotProduct xs xs )
    logBaseMeasure = isotropicNormalLogBaseMeasure

instance KnownNat n => Legendre (IsotropicNormal n) where
    potential p =
        let (nmu0,nsgma0) = split p
            nmu = coordinates nmu0
            [nsgma] = listCoordinates nsgma0
            n = fromIntegral $ natVal (Proxy @n)
            isgma = realToFrac $ recip nsgma
         in -0.25 * isgma * S.dotProduct nmu nmu - (n/2) * (log . negate $ 2 * nsgma)

instance KnownNat n => Transition Source Natural (IsotropicNormal n) where
    transition p =
        let (mu,vr0) = split p
            vr = head $ listCoordinates vr0
            nmu = breakPoint $ vr /> mu
         in join nmu . singleton $ (-0.5) * recip vr

instance KnownNat n => Transition Natural Source (IsotropicNormal n) where
    transition p =
        let (nmu,nvr0) = split p
            nvr = head $ listCoordinates nvr0
            vr = (-0.5) * recip nvr
            mu = breakPoint $ vr .> nmu
         in join mu $ singleton vr

instance KnownNat n => Transition Natural Mean (IsotropicNormal n) where
    transition p =
        let (nmu,nvr0) = split p
            nvr = S.head $ coordinates nvr0
            n = fromIntegral $ natVal (Proxy @n)
            mmu = breakPoint $ (-2 * nvr) /> nmu
            mmu0 = coordinates mmu
            mvr = singleton $ S.dotProduct mmu0 mmu0 - n/(2*nvr)
         in join mmu mvr

instance KnownNat n => Transition Mean Natural (IsotropicNormal n) where
    transition p =
        let (mmu,mvr0) = split p
            mvr = coordinates mvr0 `S.unsafeIndex` 0
            mmu0 = coordinates mmu
            n = fromIntegral $ natValInt (Proxy @n)
            nvr = recip . negate $ 2*(mvr - S.dotProduct mmu0 mmu0)/n
            nmu = breakPoint $ (-2 * nvr) .> mmu
         in join nmu $ singleton nvr


instance KnownNat n => LogLikelihood Natural (IsotropicNormal n) (S.Vector n Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance (KnownNat n, Transition Mean c (IsotropicNormal n))
  => MaximumLikelihood c (IsotropicNormal n) where
    mle = transition . averageSufficientStatistic

instance KnownNat n => AbsolutelyContinuous Natural (IsotropicNormal n) where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n, Transition c Source (IsotropicNormal n))
  => Generative c (IsotropicNormal n) where
      samplePoint p = do
          let (mus,vr) = split $ toSource p
              sd = sqrt . head $ listCoordinates vr
          S.mapM (\mu -> Random (R.normal mu sd)) $ coordinates mus

instance (KnownNat n)
  => AbsolutelyContinuous Source (IsotropicNormal n) where
      densities mvn xs = do
          let (mu,sgma0) = split mvn
              sgma = S.head $ coordinates sgma0
              n = fromIntegral $ natValInt (Proxy @n)
              scl = (2*pi*sgma)**(-n/2)
              isgma = recip sgma
          x <- xs
          let dff = x - coordinates mu
              expval = realToFrac isgma * S.dotProduct dff dff
          return $ scl * exp (-expval / 2)


-- Diagonal Normal --

instance KnownNat n => ExponentialFamily (DiagonalNormal n) where
    sufficientStatistic xs = Point $ xs S.++ (xs * xs)
    averageSufficientStatistic xs = Point $ average xs S.++ average (zipWith (*) xs xs)
    logBaseMeasure = diagonalNormalLogBaseMeasure

instance KnownNat n => Legendre (DiagonalNormal n) where
    potential p =
        let (nmu0,prcs0) = split p
            nmu = coordinates nmu0
            prcs = coordinates prcs0
            iprcs = recip prcs
         in -0.25 * S.dotProduct nmu (iprcs * nmu) -0.5 * (log . negate $ 2 * S.product prcs)

instance KnownNat n => Transition Source Natural (DiagonalNormal n) where
    transition p =
        let (mu0,vrs0) = split p
            mu = coordinates mu0
            vrs = coordinates vrs0
            nmu = Point $ mu / vrs
            prcs = Point $ -0.5 * recip vrs
         in join nmu prcs

instance KnownNat n => Transition Natural Source (DiagonalNormal n) where
    transition p =
        let (nmu0,prcs0) = split p
            nmus = coordinates nmu0
            prcs = coordinates prcs0
            vrs = Point $ -0.5 * recip prcs
            mus = Point $ -0.5 * nmus / prcs
         in join mus vrs

instance KnownNat n => Transition Natural Mean (DiagonalNormal n) where
    transition p =
        let (nmus0,prcs0) = split p
            nmus = coordinates nmus0
            prcs = coordinates prcs0
            mmus = Point $ -0.5 * nmus / prcs
            mmus0 = coordinates mmus
            mvrs = Point $ square mmus0 - 0.5 * recip prcs
         in join mmus mvrs

instance KnownNat n => Transition Mean Natural (DiagonalNormal n) where
    transition p =
        let (mmus0,mvrs0) = split p
            mvrs = coordinates mvrs0
            mmus = coordinates mmus0
            prcs0 = recip . negate $ 2 * (mvrs - square mmus)
            prcs = Point prcs0
            nmus = Point . negate $ 2 * mmus * prcs0
         in join nmus prcs


instance KnownNat n => LogLikelihood Natural (DiagonalNormal n) (S.Vector n Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance (KnownNat n, Transition Mean c (DiagonalNormal n))
  => MaximumLikelihood c (DiagonalNormal n) where
    mle = transition . averageSufficientStatistic

instance KnownNat n => AbsolutelyContinuous Natural (DiagonalNormal n) where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n, Transition c Source (DiagonalNormal n))
  => Generative c (DiagonalNormal n) where
      samplePoint p = do
          let (mus,vrs) = split $ toSource p
              sds = sqrt $ coordinates vrs
          S.zipWithM (\mu sd -> Random (R.normal mu sd)) (coordinates mus) sds

-- Multivariate Normal --

instance (KnownNat n, KnownNat (Triangular n))
  => AbsolutelyContinuous Source (MultivariateNormal n) where
      densities mvn xs = do
          let (mu,sgma) = splitMultivariateNormal mvn
              n = fromIntegral $ natValInt (Proxy @n)
              scl = (2*pi)**(-n/2) * S.determinant sgma**(-1/2)
              isgma = S.pseudoInverse sgma
          x <- xs
          let dff = x - mu
              expval = S.dotProduct dff $ S.matrixVectorMultiply isgma dff
          return $ scl * exp (-expval / 2)

instance (KnownNat n, KnownNat (Triangular n), Transition c Source (MultivariateNormal n))
  => Generative c (MultivariateNormal n) where
    sample n = sampleMultivariateNormal n . toSource

instance KnownNat n => Transition Source Natural (MultivariateNormal n) where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
            invsgma = S.pseudoInverse sgma
         in joinNaturalMultivariateNormal (S.matrixVectorMultiply invsgma mu) $ (-0.5) * invsgma

instance KnownNat n => Transition Natural Source (MultivariateNormal n) where
    transition p =
        let (nmu,nsgma) = splitNaturalMultivariateNormal p
            insgma = (-0.5) * S.pseudoInverse nsgma
         in joinMultivariateNormal (S.matrixVectorMultiply insgma nmu) insgma

instance KnownNat n => LogLikelihood Natural (MultivariateNormal n) (S.Vector n Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


instance (KnownNat n, KnownNat (Triangular n)) => ExponentialFamily (MultivariateNormal n) where
    sufficientStatistic xs = Point $ xs S.++ S.lowerTriangular (S.outerProduct xs xs)
    averageSufficientStatistic xs = Point $ average xs S.++ S.lowerTriangular ( S.averageOuterProduct $ zip xs xs )
    logBaseMeasure = multivariateNormalLogBaseMeasure

type instance PotentialCoordinates (MultivariateNormal n) = Natural

instance (KnownNat n, KnownNat (Triangular n)) => Legendre (MultivariateNormal n) where
    potential p =
        let (nmu,nsgma) = splitNaturalMultivariateNormal p
            insgma = S.pseudoInverse nsgma
         in -0.25 * S.dotProduct nmu (S.matrixVectorMultiply insgma nmu)
             -0.5 * (log . S.determinant . negate $ 2 * nsgma)

instance (KnownNat n, KnownNat (Triangular n)) => Transition Natural Mean (MultivariateNormal n) where
    transition = toMean . toSource

instance (KnownNat n, KnownNat (Triangular n)) => DuallyFlat (MultivariateNormal n) where
    dualPotential p =
        let sgma = snd . splitMultivariateNormal $ toSource p
            n = natValInt (Proxy @n)
            lndet = fromIntegral n*log (2*pi*exp 1) + log (S.determinant sgma)
         in -0.5 * lndet

instance (KnownNat n, KnownNat (Triangular n)) => Transition Mean Natural (MultivariateNormal n) where
    transition = toNatural . toSource

instance KnownNat n => Transition Source Mean (MultivariateNormal n) where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
         in joinMeanMultivariateNormal mu $ sgma + S.outerProduct mu mu

instance KnownNat n => Transition Mean Source (MultivariateNormal n) where
    transition p =
        let (mu,scnds) = splitMeanMultivariateNormal p
         in joinMultivariateNormal mu $ scnds - S.outerProduct mu mu

instance (KnownNat n, KnownNat (Triangular n)) => AbsolutelyContinuous Natural (MultivariateNormal n) where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n, Transition Mean c (MultivariateNormal n))
  => MaximumLikelihood c (MultivariateNormal n) where
    mle = transition . averageSufficientStatistic

--instance KnownNat n => MaximumLikelihood Source (MultivariateNormal n) where
--    mle _ xss =
--        let n = fromIntegral $ length xss
--            mus = recip (fromIntegral n) * sum xss
--            sgma = recip (fromIntegral $ n - 1)
--                * sum (map (\xs -> let xs' = xs - mus in M.outer xs' xs') xss)
--        in  joinMultivariateNormal mus sgma


--- Linear Models ---

instance ( KnownNat n, KnownNat k)
  => Transition Natural Source (Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)) where
    transition nfa =
        let (mvn,nmtx) = split nfa
            (nmu,nsg) = splitNaturalMultivariateNormal mvn
            invsg = -2 * nsg
            ssg = S.inverse invsg
            smu = S.matrixVectorMultiply ssg nmu
            smvn = joinMultivariateNormal smu ssg
            smtx = S.matrixMatrixMultiply ssg $ toMatrix nmtx
         in join smvn $ fromMatrix smtx

instance ( KnownNat n, KnownNat k)
  => Transition Source Natural (Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)) where
    transition lmdl =
        let (smvn,smtx) = split lmdl
            (smu,ssg) = splitMultivariateNormal smvn
            invsg = S.inverse ssg
            nmu = S.matrixVectorMultiply invsg smu
            nsg = -0.5 * invsg
            nmtx = S.matrixMatrixMultiply invsg $ toMatrix smtx
            nmvn = joinNaturalMultivariateNormal nmu nsg
         in join nmvn $ fromMatrix nmtx

instance ( KnownNat n, KnownNat k)
  => Transition Natural Source (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)) where
      transition nfa =
          let (nnrms,nmtx) = split nfa
              (nmu,nsg) = splitReplicatedProduct nnrms
              nmvn = joinNaturalMultivariateNormal (coordinates nmu) $ S.diagonalMatrix (coordinates nsg)
              nlm :: Natural # LinearModel n k
              nlm = join nmvn nmtx
              (smvn,smtx) = split $ transition nlm
              (smu,ssg) = splitMultivariateNormal smvn
              snrms = joinReplicatedProduct (Point smu) (Point $ S.takeDiagonal ssg)
           in join snrms smtx

instance ( KnownNat n, KnownNat k)
  => Transition Source Natural (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)) where
      transition sfa =
          let (snrms,smtx) = split sfa
              (smu,ssg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates $ splitReplicated snrms
              smvn = joinMultivariateNormal smu $ S.diagonalMatrix ssg
              slm :: Source # LinearModel n k
              slm = join smvn smtx
              (nmvn,nmtx) = split $ transition slm
              (nmu,nsg) = splitNaturalMultivariateNormal nmvn
              nnrms = joinReplicated $ S.zipWith (curry fromTuple) nmu $ S.takeDiagonal nsg
           in join nnrms nmtx

instance ( KnownNat n, KnownNat k)
  => Transition Source Natural (Affine Tensor (MVNMean n) (IsotropicNormal n) (MVNMean k)) where
      transition spca =
          let (iso,cwmtx) = split spca
              (cmu,cvr) = split iso
              invsg = recip . S.head $ coordinates cvr
              thtmu = Point $ realToFrac invsg * coordinates cmu
              thtsg = singleton $ (-0.5) * invsg
              imtx = fromMatrix $ realToFrac invsg * toMatrix cwmtx
           in join (join thtmu thtsg) imtx


instance Transition Natural Source (Affine Tensor NormalMean Normal NormalMean) where
      transition nfa =
          let nfa' :: Natural # LinearModel 1 1
              nfa' = breakPoint nfa
              sfa' :: Source # LinearModel 1 1
              sfa' = transition nfa'
           in breakPoint sfa'

instance Transition Source Natural (Affine Tensor NormalMean Normal NormalMean) where
      transition sfa =
          let sfa' :: Source # LinearModel 1 1
              sfa' = breakPoint sfa
              nfa' :: Natural # LinearModel 1 1
              nfa' = transition sfa'
           in breakPoint nfa'



--instance ( KnownNat n, KnownNat k)
--  => Transition Natural Source (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)) where
--    transition nfa =
--        let (nnrms,nmtx) = split nfa
--            (nmu,nsg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates
--                $ splitReplicated nnrms
--            invsg = -2 * nsg
--            ssg = recip invsg
--            smu = nmu / invsg
--            snrms = joinReplicated $ S.zipWith (curry fromTuple) smu ssg
--            smtx = S.matrixMatrixMultiply (S.diagonalMatrix ssg) $ toMatrix nmtx
--         in join snrms $ fromMatrix smtx

--instance ( KnownNat n, KnownNat k)
--  => Transition Source Natural (Affine Tensor (MVNMean n) (Replicated n Normal) (MVNMean k)) where
--    transition sfa =
--        let (snrms,smtx) = split sfa
--            (smu,ssg) = S.toPair . S.toColumns . S.fromRows . S.map coordinates
--                $ splitReplicated snrms
--            invsg = recip ssg
--            nmu = invsg * smu
--            nsg = -0.5 * invsg
--            nmtx = S.matrixMatrixMultiply (S.diagonalMatrix invsg) $ toMatrix smtx
--            nnrms = joinReplicated $ S.zipWith (curry fromTuple) nmu nsg
--         in join nnrms $ fromMatrix nmtx



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
    , MultivariateNormal
    , SymmetricNormal
    , DiagonalNormal
    , IsotropicNormal
    , multivariateNormalCorrelations
    , bivariateNormalConfidenceEllipse
    -- * Linear Models
    , SimpleLinearModel
    , LinearModel
    , FactorAnalysis
    , PrincipleComponentAnalysis
    ) where

-- Package --

import Goal.Core
import Goal.Probability.Statistical
import Goal.Probability.Conditional
import Goal.Probability.ExponentialFamily
import Goal.Probability.Distributions

import Goal.Geometry

import qualified Goal.Core.Vector.Storable as S

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
type MultivariateNormal f (n :: Nat) = LocationShape (MVNMean n) (f (MVNMean n) (MVNMean n))

type SymmetricNormal n = MultivariateNormal Symmetric n
type DiagonalNormal n = MultivariateNormal Diagonal n
type IsotropicNormal n = MultivariateNormal Scale n

-- | Linear models are linear functions with additive Guassian noise.
type SimpleLinearModel = Affine Tensor NormalMean Normal NormalMean

type FactorAnalysis n k = Affine Tensor (MVNMean n) (DiagonalNormal n) (MVNMean k)
type PrincipleComponentAnalysis n k = Affine Tensor (MVNMean n) (IsotropicNormal n) (MVNMean k)

-- | Linear models are linear functions with additive Guassian noise.
type LinearModel n k = Affine Tensor (MVNMean n) (MultivariateNormal Symmetric n) (MVNMean k)

---- | Split a MultivariateNormal into its Means and Covariance matrix.
--splitMultivariateNormal
--    :: KnownNat n
--    => Source # MultivariateNormal n
--    -> (S.Vector n Double, S.Matrix n n Double)
--splitMultivariateNormal mvn =
--    let (mu,cvr) = split mvn
--     in (coordinates mu, S.fromLowerTriangular $ coordinates cvr)
--
---- | Join a covariance matrix into a MultivariateNormal.
--joinMultivariateNormal
--    :: KnownNat n
--    => S.Vector n Double
--    -> S.Matrix n n Double
--    -> Source # MultivariateNormal n
--joinMultivariateNormal mus sgma =
--    join (Point mus) (Point $ S.lowerTriangular sgma)
--
---- | Split a MultivariateNormal into its Means and Covariance matrix.
--splitMeanMultivariateNormal
--    :: KnownNat n
--    => Mean # MultivariateNormal n
--    -> (S.Vector n Double, S.Matrix n n Double)
--splitMeanMultivariateNormal mvn =
--    let (mu,cvr) = split mvn
--     in (coordinates mu, S.fromLowerTriangular $ coordinates cvr)
--
---- | Join a covariance matrix into a MultivariateNormal.
--joinMeanMultivariateNormal
--    :: KnownNat n
--    => S.Vector n Double
--    -> S.Matrix n n Double
--    -> Mean # MultivariateNormal n
--joinMeanMultivariateNormal mus sgma =
--    join (Point mus) (Point $ S.lowerTriangular sgma)
--

naturalSymmetricToPrecision
    :: forall n . KnownNat n
    => Natural # Symmetric (MVNMean n) (MVNMean n)
    -> Natural # Tensor (MVNMean n) (MVNMean n)
naturalSymmetricToPrecision trng =
    let tns = toTensor trng
        tns' = 2 /> tns
        diag :: Natural # Diagonal (MVNMean n) (MVNMean n)
        diag = fromTensor tns'
     in tns' + toTensor diag

naturalPrecisionToSymmetric
    :: forall n . KnownNat n
    => Natural # Tensor (MVNMean n) (MVNMean n)
    -> Natural # Symmetric (MVNMean n) (MVNMean n)
naturalPrecisionToSymmetric tns =
    let diag :: Natural # Diagonal (MVNMean n) (MVNMean n)
        diag = fromTensor tns
     in fromTensor $ 2 .> tns - toTensor diag


---- | Split a MultivariateNormal into the precision weighted means and (-0.5*)
---- Precision matrix. Note that this performs an easy to miss computation for
---- converting the natural parameters in our reduced representation of MVNs into
---- the full precision matrix.
--splitNaturalMultivariateNormal
--    :: KnownNat n
--    => Natural # MultivariateNormal n
--    -> (S.Vector n Double, S.Matrix n n Double)
--splitNaturalMultivariateNormal np =
--    let (nmu,cvrs) = split np
--        nmu0 = coordinates nmu
--        nsgma0' = (/2) . S.fromLowerTriangular $ coordinates cvrs
--        nsgma0 = nsgma0' + S.diagonalMatrix (S.takeDiagonal nsgma0')
--     in (nmu0, nsgma0)
--
---- | Joins a MultivariateNormal out of the precision weighted means and (-0.5)
---- Precision matrix. Note that this performs an easy to miss computation for
---- converting the full precision Matrix into the reduced, EF representation we use here.
--joinNaturalMultivariateNormal
--    :: KnownNat n
--    => S.Vector n Double
--    -> S.Matrix n n Double
--    -> Natural # MultivariateNormal n
--joinNaturalMultivariateNormal nmu0 nsgma0 =
--    let nmu = Point nmu0
--        diag = S.diagonalMatrix $ S.takeDiagonal nsgma0
--     in join nmu . Point . S.lowerTriangular $ 2*nsgma0 - diag
--
-- | Confidence elipses for bivariate normal distributions.
bivariateNormalConfidenceEllipse
    :: ( Square f (MVNMean 2) (MVNMean 2), Map Source f (MVNMean 2) (MVNMean 2) )
    => Int
    -> Double
    -> Source # MultivariateNormal f 2
    -> [(Double,Double)]
bivariateNormalConfidenceEllipse nstps prcnt mvn =
    let (mu,sgma) = split mvn
        mrt = prcnt .> matrixRoot sgma
        xs = range 0 (2*pi) nstps
        sxs = [ fromTuple (cos x, sin x) | x <- xs ]
     in S.toPair . coordinates . (mu +) <$> mrt >$> sxs

-- | Computes the correlation matrix of a 'MultivariateNormal' distribution.
multivariateNormalCorrelations
    :: forall f n . (KnownNat n, Square f (MVNMean n) (MVNMean n))
    => Source # MultivariateNormal f n
    -> Source # Tensor (MVNMean n) (MVNMean n)
multivariateNormalCorrelations mvn =
    let cvrs = toTensor . snd $ split mvn
        diag :: Source # Diagonal (MVNMean n) (MVNMean n)
        diag = fromTensor cvrs
        sds = breakPoint $ sqrt diag
        sdmtx = sds >.< sds
     in cvrs / sdmtx

multivariateNormalLogBaseMeasure
    :: forall f n . (KnownNat n)
    => Proxy (MultivariateNormal f n)
    -> S.Vector n Double
    -> Double
multivariateNormalLogBaseMeasure _ _ =
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
    :: (KnownNat n, Square f (MVNMean n) (MVNMean n), Map Source f (MVNMean n) (MVNMean n))
    => Int
    -> Source # MultivariateNormal f n
    -> Random [S.Vector n Double]
sampleMultivariateNormal n p = do
    let (mu,sgma) = split p
        rtsgma = matrixRoot sgma
    x0s <- replicateM n . S.replicateM $ Random (R.normal 0 1)
    return $ coordinates . (mu +) <$> rtsgma >$> (Point <$> x0s)

--isotropicNormalToFull
--    :: KnownNat n
--    => Natural # IsotropicNormal n
--    -> Natural # MultivariateNormal n
--isotropicNormalToFull iso =
--    let (mus,sgma0) = split iso
--        sgma = realToFrac . S.head $ coordinates sgma0
--     in joinNaturalMultivariateNormal (coordinates mus) $ sgma * S.matrixIdentity
--
--fullNormalToIsotropic
--    :: KnownNat n
--    => Mean # MultivariateNormal n
--    -> Mean # IsotropicNormal n
--fullNormalToIsotropic iso =
--    let (mus,sgma0) = splitMeanMultivariateNormal iso
--        sgma = S.sum $ S.takeDiagonal sgma0
--     in join (Point mus) $ singleton sgma
--
--diagonalNormalToFull
--    :: KnownNat n
--    => Natural # DiagonalNormal n
--    -> Natural # MultivariateNormal n
--diagonalNormalToFull diag =
--    let (mus,prcs) = split diag
--     in joinNaturalMultivariateNormal (coordinates mus) . S.diagonalMatrix $ coordinates prcs
--
--fullNormalToDiagonal
--    :: KnownNat n
--    => Mean # MultivariateNormal n
--    -> Mean # DiagonalNormal n
--fullNormalToDiagonal diag =
--    let (mus,sgma) = splitMeanMultivariateNormal diag
--     in join (Point mus) . Point $ S.takeDiagonal sgma

-- Restricted MVNs --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The 'Source'
-- coordinates are the (vector) mean and the covariance matrix. For the
-- coordinates of a multivariate normal distribution, the elements of the mean
-- come first, and then the elements of the covariance matrix in row major
-- order.
--type IsotropicNormal (n :: Nat) = LocationShape (MVNMean n) NormalVariance
--
--type DiagonalNormal (n :: Nat) = LocationShape (MVNMean n) (Replicated n NormalVariance)


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


--- MVNMean ---


instance KnownNat n => Manifold (MVNMean n) where
    type Dimension (MVNMean n) = n

instance (KnownNat n) => Statistical (MVNMean n) where
    type SamplePoint (MVNMean n) = S.Vector n Double

instance KnownNat n => ExponentialFamily (MVNMean n) where
    sufficientStatistic x = Point x
    logBaseMeasure = mvnMeanLogBaseMeasure

type instance PotentialCoordinates (MVNMean n) = Natural


--- MultivariateNormal ---


--- Transition Instances

type instance PotentialCoordinates (MultivariateNormal f n) = Natural

instance KnownNat n => Transition Source Natural (SymmetricNormal n) where
    transition p =
        let (mu,sgma) = split p
            invsgma = inverse $ toTensor sgma
         in join (breakPoint $ invsgma >.> mu) . naturalPrecisionToSymmetric
             . breakPoint $ (-0.5) * invsgma

instance KnownNat n => Transition Natural Source (SymmetricNormal n) where
    transition p =
        let (nmu,ntrng) = split p
            nsgma = naturalSymmetricToPrecision ntrng
            insgma = (-0.5) .> inverse nsgma
            trng :: Mean # Symmetric (MVNMean n) (MVNMean n)
            trng = fromTensor insgma
         in join (breakPoint $ insgma >.> nmu) $ breakPoint trng

instance KnownNat n => Transition Source Natural (DiagonalNormal n) where
    transition p =
        let (mu,sgma) = split p
            invsgma = inverse sgma
         in join (breakPoint $ invsgma >.> mu) . breakPoint $ (-0.5) * invsgma

instance KnownNat n => Transition Natural Source (DiagonalNormal n) where
    transition p =
        let (nmu,nsgma) = split p
            insgma = (-0.5) .> inverse nsgma
         in join (breakPoint $ insgma >.> nmu) $ breakPoint insgma

instance KnownNat n => Transition Source Natural (IsotropicNormal n) where
    transition p =
        let (mu,sgma) = split p
            invsgma = inverse sgma
         in join (breakPoint $ invsgma >.> mu) . breakPoint $ (-0.5) * invsgma

instance KnownNat n => Transition Natural Source (IsotropicNormal n) where
    transition p =
        let (nmu,nsgma) = split p
            insgma = (-0.5) .> inverse nsgma
         in join (breakPoint $ insgma >.> nmu) $ breakPoint insgma

instance ( KnownNat n, Square f (MVNMean n) (MVNMean n) )
  => Transition Source Mean (MultivariateNormal f n) where
    transition mvn =
        let (mu,sgma) = split mvn
            mmvn :: Source # MultivariateNormal f n
            mmvn = join mu $ sgma + (mu >.< mu)
         in breakPoint mmvn

instance ( KnownNat n, Square f (MVNMean n) (MVNMean n) )
  => Transition Mean Source (MultivariateNormal f n) where
    transition mmvn =
        let (mu,msgma) = split mmvn
            mvn :: Mean # MultivariateNormal f n
            mvn = join mu $ msgma - (mu >.< mu)
         in breakPoint mvn

instance KnownNat n => Transition Natural Mean (SymmetricNormal n) where
    transition = toMean . toSource

instance KnownNat n => Transition Natural Mean (DiagonalNormal n) where
    transition = toMean . toSource

instance KnownNat n => Transition Natural Mean (IsotropicNormal n) where
    transition = toMean . toSource

instance KnownNat n => Transition Mean Natural (SymmetricNormal n) where
    transition = toNatural . toSource

instance KnownNat n => Transition Mean Natural (DiagonalNormal n) where
    transition = toNatural . toSource

instance KnownNat n => Transition Mean Natural (IsotropicNormal n) where
    transition = toNatural . toSource

--- Basic Instances

instance (KnownNat n, Square f (MVNMean n) (MVNMean n), Map Source f (MVNMean n) (MVNMean n))
  => AbsolutelyContinuous Source (MultivariateNormal f n) where
      densities mvn xs = do
          let (mu,sgma) = split mvn
              n = fromIntegral $ natValInt (Proxy @n)
              scl = (2*pi)**(-n/2) * determinant sgma**(-1/2)
              isgma = inverse sgma
          x <- xs
          let dff = Point $ x - coordinates mu
              expval = dff <.> (isgma >.> dff)
          return $ scl * exp (-expval / 2)

instance ( KnownNat n, Square f (MVNMean n) (MVNMean n), Map Source f (MVNMean n) (MVNMean n)
         , Transition c Source (MultivariateNormal f n) )
  => Generative c (MultivariateNormal f n) where
    sample n = sampleMultivariateNormal n . toSource

--- Exponential Family Instances

instance (KnownNat n, Square f (MVNMean n) (MVNMean n), Bilinear f (MVNMean n) (MVNMean n))
  => ExponentialFamily (MultivariateNormal f n) where
    sufficientStatistic x =
        let mx = sufficientStatistic x
         in join mx $ mx >.< mx
    averageSufficientStatistic xs =
        let mxs = sufficientStatistic <$> xs
         in join (average mxs) $ mxs >$< mxs
    logBaseMeasure = multivariateNormalLogBaseMeasure

instance KnownNat n => Legendre (SymmetricNormal n) where
    potential p =
        let (nmu,nsgma) = split p
            (insgma,dtmnt,_) = inverseLogDeterminant . negate $ 2 * (naturalSymmetricToPrecision nsgma)
         in 0.5 * (nmu <.> (insgma >.> nmu)) -0.5 * dtmnt

instance KnownNat n => Legendre (DiagonalNormal n) where
    potential p =
        let (nmu,nsgma) = split p
            (insgma,dtmnt,_) = inverseLogDeterminant . negate $ 2 * nsgma
         in 0.5 * (nmu <.> (insgma >.> nmu)) -0.5 * dtmnt

instance KnownNat n => Legendre (IsotropicNormal n) where
    potential p =
        let (nmu,nsgma) = split p
            (insgma,dtmnt,_) = inverseLogDeterminant . negate $ 2 * nsgma
         in 0.5 * (nmu <.> (insgma >.> nmu)) -0.5 * dtmnt

instance ( KnownNat n, Square f (MVNMean n) (MVNMean n), Legendre (MultivariateNormal f n) )
  => DuallyFlat (MultivariateNormal f n) where
    dualPotential p =
        let sgma = snd . split $ toSource p
            n = natValInt (Proxy @n)
            (_,lndet0,_) = inverseLogDeterminant sgma
            lndet = fromIntegral n*log (2*pi*exp 1) + lndet0
         in -0.5 * lndet

instance ( KnownNat n, Legendre (MultivariateNormal f n), Square f (MVNMean n) (MVNMean n)
         , Transition Natural Mean (MultivariateNormal f n) )
  => LogLikelihood Natural (MultivariateNormal f n) (S.Vector n Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance ( KnownNat n, Legendre (MultivariateNormal f n), Square f (MVNMean n) (MVNMean n)
         , Transition Natural Mean (MultivariateNormal f n) )
  => AbsolutelyContinuous Natural (MultivariateNormal f n) where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n, Transition Mean c (MultivariateNormal f n), Square f (MVNMean n) (MVNMean n) )
  => MaximumLikelihood c (MultivariateNormal f n) where
    mle = transition . averageSufficientStatistic

----- Linear Models ---
--
--instance ( KnownNat n, KnownNat k)
--  => Transition Natural Source (Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)) where
--    transition nfa =
--        let (mvn,nmtx) = split nfa
--            (nmu,nsg) = splitNaturalMultivariateNormal mvn
--            invsg = -2 * nsg
--            ssg = S.inverse invsg
--            smu = S.matrixVectorMultiply ssg nmu
--            smvn = joinMultivariateNormal smu ssg
--            smtx = S.matrixMatrixMultiply ssg $ toMatrix nmtx
--         in join smvn $ fromMatrix smtx
--
--instance ( KnownNat n, KnownNat k)
--  => Transition Source Natural (Affine Tensor (MVNMean n) (MultivariateNormal n) (MVNMean k)) where
--    transition lmdl =
--        let (smvn,smtx) = split lmdl
--            (smu,ssg) = splitMultivariateNormal smvn
--            invsg = S.inverse ssg
--            nmu = S.matrixVectorMultiply invsg smu
--            nsg = -0.5 * invsg
--            nmtx = S.matrixMatrixMultiply invsg $ toMatrix smtx
--            nmvn = joinNaturalMultivariateNormal nmu nsg
--         in join nmvn $ fromMatrix nmtx
--
--instance ( KnownNat n, KnownNat k)
--  => Transition Natural Source (FactorAnalysis n k) where
--      transition nfa =
--          let (nnrms,nmtx) = split nfa
--              (nmu,nsg) = split nnrms
--              nmvn = joinNaturalMultivariateNormal (coordinates nmu) $ S.diagonalMatrix (coordinates nsg)
--              nlm :: Natural # LinearModel n k
--              nlm = join nmvn nmtx
--              (smvn,smtx) = split $ transition nlm
--              (smu,ssg) = splitMultivariateNormal smvn
--              snrms = join (Point smu) (Point $ S.takeDiagonal ssg)
--           in join snrms smtx
--
--instance ( KnownNat n, KnownNat k)
--  => Transition Source Natural (FactorAnalysis n k) where
--      transition sfa =
--          let (snrms,smtx) = split sfa
--              (smu,ssg) = split snrms
--              smvn = joinMultivariateNormal (coordinates smu) . S.diagonalMatrix $ coordinates ssg
--              slm :: Source # LinearModel n k
--              slm = join smvn smtx
--              (nmvn,nmtx) = split $ transition slm
--              (nmu,nsg) = splitNaturalMultivariateNormal nmvn
--              nnrms = join (Point nmu) . Point $ S.takeDiagonal nsg
--           in join nnrms nmtx
--
--instance ( KnownNat n, KnownNat k)
--  => Transition Source Natural (PrincipleComponentAnalysis n k) where
--      transition spca =
--          let (iso,cwmtx) = split spca
--              (cmu,cvr) = split iso
--              invsg = recip . S.head $ coordinates cvr
--              thtmu = Point $ realToFrac invsg * coordinates cmu
--              thtsg = singleton $ (-0.5) * invsg
--              imtx = fromMatrix $ realToFrac invsg * toMatrix cwmtx
--           in join (join thtmu thtsg) imtx
--
--
--instance Transition Natural Source (Affine Tensor NormalMean Normal NormalMean) where
--      transition nfa =
--          let nfa' :: Natural # LinearModel 1 1
--              nfa' = breakPoint nfa
--              sfa' :: Source # LinearModel 1 1
--              sfa' = transition nfa'
--           in breakPoint sfa'
--
--instance Transition Source Natural (Affine Tensor NormalMean Normal NormalMean) where
--      transition sfa =
--          let sfa' :: Source # LinearModel 1 1
--              sfa' = breakPoint sfa
--              nfa' :: Natural # LinearModel 1 1
--              nfa' = transition sfa'
--           in breakPoint nfa'



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



{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE UndecidableInstances #-}

-- | Various instances of statistical manifolds, with a focus on exponential
-- families. In the documentation we use \(X\) to indicate a random variable
-- with the distribution being documented.
module Goal.Probability.Distributions
    ( -- * Univariate
      Bernoulli
    , Binomial
    , Categorical
    , categoricalWeights
    , Poisson
    , Normal
    , NormalMean
    , NormalShape
    , VonMises
    -- * Multivariate
    , Dirichlet
    , MultivariateNormal
    , joinMultivariateNormal
    , splitMultivariateNormal
    , multivariateNormalCorrelations
    , bivariateNormalConfidenceEllipse
    -- * LocationShape
    , LocationShape (LocationShape)
    ) where

-- Package --

import Goal.Core
import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily

import Goal.Geometry
import System.Random.MWC.Probability hiding (sample)

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Numeric.GSL.Special.Bessel as GSL
import qualified Numeric.GSL.Special.Gamma as GSL
import qualified Numeric.GSL.Special.Psi as GSL

import Foreign.Storable

-- Location Shape --

-- | A 'LocationShape' 'Manifold' is a 'Product' of some location 'Manifold' and some shape 'Manifold'.
newtype LocationShape l s = LocationShape (l,s)

deriving instance (Manifold l, Manifold s) => Manifold (LocationShape l s)
deriving instance (Manifold l, Manifold s) => Product (LocationShape l s)

-- Uniform --

-- Bernoulli Distribution --

-- | The Bernoulli family with 'Bool'ean 'SamplePoint's. (because why not). The source coordinate is \(P(X = True)\).
data Bernoulli

-- Binomial Distribution --

-- | A distribution over the sum of 'True' realizations of @n@ 'Bernoulli'
-- random variables. The 'Source' coordinate is the probability of \(P(X = True)\)
-- for each 'Bernoulli' random variable.
data Binomial (n :: Nat)

-- | Returns the number of trials used to define this binomial distribution.
binomialTrials :: forall c n. KnownNat n => Point c (Binomial n) -> Int
binomialTrials _ = natValInt (Proxy :: Proxy n)

-- | Returns the number of trials used to define this binomial distribution.
binomialSampleSpace :: forall n . KnownNat n => Proxy (Binomial n) -> Int
binomialSampleSpace _ = natValInt (Proxy :: Proxy n)

-- Categorical Distribution --

-- | A 'Categorical' distribution where the probability of the first category
-- \(P(X = 0)\) is given by the normalization constraint.
data Categorical (n :: Nat)

-- | Takes a weighted list of elements representing a probability mass function, and
-- returns a sample from the Categorical distribution.
sampleCategorical :: KnownNat n => S.Vector n Double -> Random r Int
sampleCategorical ps = do
    let ps' = S.postscanl' (+) 0 ps
    p <- uniform
    let midx = (+1) . finiteInt <$> S.findIndex (> p) ps'
    return $ fromMaybe 0 midx

-- | Returns the probabilities over the whole sample space \((0 \ldots n)\) of the
-- given categorical distribution.
categoricalWeights
    :: Transition c Source (Categorical n)
    => c # Categorical n
    -> S.Vector (n+1) Double
categoricalWeights wghts0 =
    let wghts = coordinates $ toSource wghts0
     in S.cons (1-S.sum wghts) wghts

-- | A 'Dirichlet' manifold contains distributions over weights of a
-- 'Categorical' distribution.
data Dirichlet (k :: Nat)

-- Poisson Distribution --

-- | Returns a sample from a Poisson distribution with the given rate.
samplePoisson :: Double -> Random s Int
samplePoisson lmda = uniform >>= renew 0
    where l = exp (-lmda)
          renew k p
            | p <= l = return k
            | otherwise = do
                u <- uniform
                renew (k+1) (p*u)

-- | The 'Manifold' of 'Poisson' distributions. The 'Source' coordinate is the
-- rate of the Poisson distribution.
data Poisson

-- Poisson Distribution --

-- Normal Distribution --

-- | The Mean of a normal distribution.
data NormalMean

-- | The variance of a normal distribution.
data NormalShape
-- | The 'Manifold' of 'Normal' distributions. The 'Source' coordinates are the
-- mean and the variance.
type Normal = LocationShape NormalMean NormalShape


-- Multivariate Normal --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The 'Source'
-- coordinates are the (vector) mean and the covariance matrix. For the
-- coordinates of a multivariate normal distribution, the elements of the mean
-- come first, and then the elements of the covariance matrix in row major
-- order.
data MultivariateNormal (n :: Nat)

splitMultivariateNormal0
    :: KnownNat n
    => c # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitMultivariateNormal0 (Point xs) =
    let (mus,cvrs) = S.splitAt xs
     in (mus,S.fromLowerTriangular cvrs)

joinMultivariateNormal0
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> c # MultivariateNormal n
joinMultivariateNormal0 mus sgma =
    Point $ mus S.++ S.lowerTriangular sgma

splitNaturalMultivariateNormal
    :: KnownNat n
    => Natural # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitNaturalMultivariateNormal np =
    let (nmu,nsgma) = splitMultivariateNormal0 np
     in (nmu, (+ scaleMatrix 0.5 nsgma) . scaleMatrix 0.5 . S.diagonalMatrix $ S.takeDiagonal nsgma)

joinNaturalMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Natural # MultivariateNormal n
joinNaturalMultivariateNormal nmu nsgma =
    let nsgma' = (+ scaleMatrix 2 nsgma) . scaleMatrix (-1) . S.diagonalMatrix $ S.takeDiagonal nsgma
     in joinMultivariateNormal0 nmu nsgma'

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


--unsafeCholesky
--    :: (KnownNat n, Field x, Storable x)
--    => Matrix n n x
--    -> Matrix n n x
--unsafeCholesky =
--    transpose . fromHMatrix . H.chol . H.trustSym . toHMatrix


-- | Splits a 'MultivariateNormal' distribution in 'Source' coordinates into its
-- mean vector and covariance matrix.
splitMultivariateNormal
    :: KnownNat n
    => Source # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
splitMultivariateNormal = splitMultivariateNormal0

-- | Constructs a 'MultivariateNormal' distribution in 'Source' coordinates from
-- a mean vector and covariance matrix.
joinMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Source # MultivariateNormal n
joinMultivariateNormal mus sgma =
    Point $ mus S.++ S.lowerTriangular sgma

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
     in log $ pi**(-fromIntegral n/2)

-- | samples a multivariateNormal by way of a covariance matrix i.e. by taking
-- the square root.
sampleMultivariateNormal
    :: KnownNat n
    => Source # MultivariateNormal n
    -> Random s (S.Vector n Double)
sampleMultivariateNormal p = do
    let (mus,sgma) = splitMultivariateNormal p
    nrms <- S.replicateM $ normal 0 1
    let rtsgma = S.matrixRoot sgma
    return $ mus + S.matrixVectorMultiply rtsgma nrms

-- von Mises --

-- | The 'Manifold' of 'VonMises' distributions. The 'Source' coordinates are
-- the mean and concentration.
data VonMises


--- Internal ---


binomialLogBaseMeasure0 :: (KnownNat n) => Proxy n -> Proxy (Binomial n) -> Int -> Double
binomialLogBaseMeasure0 prxyn _ = logChoose (natValInt prxyn)


--- Instances ---


-- Bernoulli Distribution --

instance Manifold Bernoulli where
    type Dimension Bernoulli = 1

instance Statistical Bernoulli where
    type (SamplePoint Bernoulli) = Bool

instance Discrete Bernoulli where
    type Cardinality Bernoulli = 2
    sampleSpace _ = [True,False]

instance ExponentialFamily Bernoulli where
    logBaseMeasure _ _ = 0
    sufficientStatistic True = Point $ S.singleton 1
    sufficientStatistic False = Point $ S.singleton 0

type instance PotentialCoordinates Bernoulli = Natural

instance Legendre Bernoulli where
    potential p = log $ 1 + exp (S.head $ coordinates p)

--instance {-# OVERLAPS #-} KnownNat k => Legendre (Replicated k Bernoulli) where
--    potential p = S.sum . S.map (log . (1 +) .  exp) $ coordinates p

instance Transition Natural Mean Bernoulli where
    transition = Point . S.map logistic . coordinates

instance DuallyFlat Bernoulli where
    dualPotential p =
        let eta = S.head $ coordinates p
         in logit eta * eta - log (1 / (1 - eta))

instance Transition Mean Natural Bernoulli where
    transition = Point . S.map logit . coordinates

instance Riemannian Natural Bernoulli where
    metric p =
        let stht = logistic . S.head $ coordinates p
         in Point . S.singleton $ stht * (1-stht)
    flat p p' =
        let stht = logistic . S.head $ coordinates p
         in breakPoint $ (stht * (1-stht)) .> p'

instance {-# OVERLAPS #-} KnownNat k => Riemannian Natural (Replicated k Bernoulli) where
    metric = error "Do not call metric on a replicated manifold"
    flat p p' =
        let sthts = S.map ((\stht -> stht * (1-stht)) . logistic) $ coordinates p
            dp = S.zipWith (*) sthts $ coordinates p'
         in Point dp

instance {-# OVERLAPS #-} KnownNat k => Riemannian Mean (Replicated k Bernoulli) where
    metric = error "Do not call metric on a replicated manifold"
    sharp p dp =
        let sthts' = S.map (\stht -> stht * (1-stht)) $ coordinates p
            p' = S.zipWith (*) sthts' $ coordinates dp
         in Point p'

instance Transition Source Mean Bernoulli where
    transition = breakPoint

instance Transition Mean Source Bernoulli where
    transition = breakPoint

instance Transition Source Natural Bernoulli where
    transition = transition . toMean

instance Transition Natural Source Bernoulli where
    transition = transition . toMean

instance (Transition c Source Bernoulli) => Generative c Bernoulli where
    samplePoint = bernoulli . S.head . coordinates . toSource

instance Transition Mean c Bernoulli => MaximumLikelihood c Bernoulli where
    mle = transition . averageSufficientStatistic

instance LogLikelihood Natural Bernoulli Bool where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance AbsolutelyContinuous Source Bernoulli where
    densities sb bs =
        let p = S.head $ coordinates sb
         in [ if b then p else 1 - p | b <- bs ]

instance AbsolutelyContinuous Mean Bernoulli where
    densities = densities . toSource

instance AbsolutelyContinuous Natural Bernoulli where
    logDensities = exponentialFamilyLogDensities

-- Binomial Distribution --

instance KnownNat n => Manifold (Binomial n) where
    type Dimension (Binomial n) = 1

instance KnownNat n => Statistical (Binomial n) where
    type SamplePoint (Binomial n) = Int

instance KnownNat n => Discrete (Binomial n) where
    type Cardinality (Binomial n) = n + 1
    sampleSpace prx = [0..binomialSampleSpace prx]

instance KnownNat n => ExponentialFamily (Binomial n) where
    logBaseMeasure = binomialLogBaseMeasure0 Proxy
    sufficientStatistic = Point . S.singleton . fromIntegral

type instance PotentialCoordinates (Binomial n) = Natural

instance KnownNat n => Legendre (Binomial n) where
    potential p =
        let n = fromIntegral $ binomialTrials p
            tht = S.head $ coordinates p
         in n * log (1 + exp tht)

instance KnownNat n => Transition Natural Mean (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
         in Point . S.singleton $ n * logistic (S.head $ coordinates p)

instance KnownNat n => DuallyFlat (Binomial n) where
    dualPotential p =
        let n = fromIntegral $ binomialTrials p
            eta = S.head $ coordinates p
        in eta * log (eta / (n - eta)) - n * log (n / (n - eta))

instance KnownNat n => Transition Mean Natural (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
            eta = S.head $ coordinates p
         in Point . S.singleton . log $ eta / (n - eta)

instance KnownNat n => Transition Source Natural (Binomial n) where
    transition = transition . toMean

instance KnownNat n => Transition Natural Source (Binomial n) where
    transition = transition . toMean

instance KnownNat n => Transition Source Mean (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
         in breakPoint $ n .> p

instance KnownNat n => Transition Mean Source (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
         in breakPoint $ n /> p

instance (KnownNat n, Transition c Source (Binomial n)) => Generative c (Binomial n) where
    samplePoint p0 = do
        let p = toSource p0
            n = binomialTrials p
        bls <- replicateM n . bernoulli . S.head $ coordinates p
        return $ sum [ if bl then 1 else 0 | bl <- bls ]

instance KnownNat n => AbsolutelyContinuous Source (Binomial n) where
    densities p ks =
        let n = binomialTrials p
            c = S.head $ coordinates p
         in [ choose n k * c^k * (1 - c)^(n-k) | k <- ks ]

instance KnownNat n => AbsolutelyContinuous Mean (Binomial n) where
    densities = densities . toSource

instance KnownNat n => AbsolutelyContinuous Natural (Binomial n) where
    logDensities = exponentialFamilyLogDensities

instance (KnownNat n, Transition Mean c (Binomial n)) => MaximumLikelihood c (Binomial n) where
    mle = transition . averageSufficientStatistic

instance KnownNat n => LogLikelihood Natural (Binomial n) Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- Categorical Distribution --

instance KnownNat n => Manifold (Categorical n) where
    type Dimension (Categorical n) = n

instance KnownNat n => Statistical (Categorical n) where
    type SamplePoint (Categorical n) = Int

instance KnownNat n => Discrete (Categorical n) where
    type Cardinality (Categorical n) = n
    sampleSpace prx = [0..dimension prx]

instance KnownNat n => ExponentialFamily (Categorical n) where
    logBaseMeasure _ _ = 0
    sufficientStatistic e = Point $ S.generate (\i -> if finiteInt i == (fromEnum e-1) then 1 else 0)

type instance (PotentialCoordinates (Categorical n)) = Natural

instance KnownNat n => Legendre (Categorical n) where
    --potential (Point cs) = log $ 1 + S.sum (S.map exp cs)
    potential = logSumExp . B.cons 0 . boxCoordinates

instance KnownNat n => Transition Natural Mean (Categorical n) where
    transition p =
        let exps = S.map exp $ coordinates p
            nrm = 1 + S.sum exps
         in nrm /> Point exps

instance KnownNat n => DuallyFlat (Categorical n) where
    dualPotential (Point cs) =
        let sc = 1 - S.sum cs
         in S.sum (S.map entropyFun cs) + entropyFun sc
        where entropyFun 0 = 0
              entropyFun x = x * log x

instance KnownNat n => Transition Mean Natural (Categorical n) where
    transition (Point xs) =
        let nrm = 1 - S.sum xs
         in  Point . log $ S.map (/nrm) xs

instance Transition Source Mean (Categorical n) where
    transition = breakPoint

instance Transition Mean Source (Categorical n) where
    transition = breakPoint

instance KnownNat n => Transition Source Natural (Categorical n) where
    transition = transition . toMean

instance KnownNat n => Transition Natural Source (Categorical n) where
    transition = transition . toMean

instance (KnownNat n, Transition c Source (Categorical n))
  => Generative c (Categorical n) where
    samplePoint p0 =
        let p = toSource p0
         in sampleCategorical $ coordinates p

instance (KnownNat n, Transition Mean c (Categorical n))
  => MaximumLikelihood c (Categorical n) where
    mle = transition . averageSufficientStatistic

instance KnownNat n => LogLikelihood Natural (Categorical n) Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


instance KnownNat n => AbsolutelyContinuous Source (Categorical n) where
    densities (Point ps) es = do
        e <- es
        let ek = fromEnum e
            p0 = 1 - S.sum ps
        return $ if ek == 0
                    then p0
                    else S.unsafeIndex ps $ ek - 1

instance KnownNat n => AbsolutelyContinuous Mean (Categorical n) where
    densities = densities . toSource

instance KnownNat n => AbsolutelyContinuous Natural (Categorical n) where
    logDensities = exponentialFamilyLogDensities

-- Dirichlet Distribution --

instance KnownNat k => Manifold (Dirichlet k) where
    type Dimension (Dirichlet k) = k

instance KnownNat k => Statistical (Dirichlet k) where
    type SamplePoint (Dirichlet k) = S.Vector k Double

instance (KnownNat k, Transition c Source (Dirichlet k))
  => Generative c (Dirichlet k) where
    samplePoint p0 = do
        let alphs = boxCoordinates $ toSource p0
        G.convert <$> dirichlet alphs

instance KnownNat k => ExponentialFamily (Dirichlet k) where
    logBaseMeasure _ = negate . S.sum
    sufficientStatistic xs = Point $ S.map log xs

logMultiBeta :: KnownNat k => S.Vector k Double -> Double
logMultiBeta alphs =
    S.sum (S.map GSL.lngamma alphs) - GSL.lngamma (S.sum alphs)

logMultiBetaDifferential :: KnownNat k => S.Vector k Double -> S.Vector k Double
logMultiBetaDifferential alphs =
    S.map (subtract (GSL.psi $ S.sum alphs) . GSL.psi) alphs

type instance PotentialCoordinates (Dirichlet k) = Natural

instance KnownNat k => Legendre (Dirichlet k) where
    potential = logMultiBeta . coordinates

instance KnownNat k => Transition Natural Mean (Dirichlet k) where
    transition = Point . logMultiBetaDifferential . coordinates

instance KnownNat k => AbsolutelyContinuous Source (Dirichlet k) where
    densities p xss = do
        xs <- xss
        let alphs = coordinates p
            prds = S.product $ S.zipWith (**) xs $ S.map (subtract 1) alphs
        return $ prds / exp (logMultiBeta alphs)

instance KnownNat k => AbsolutelyContinuous Natural (Dirichlet k) where
    logDensities = exponentialFamilyLogDensities

instance KnownNat k => LogLikelihood Natural (Dirichlet k) (S.Vector k Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance KnownNat k => Transition Source Natural (Dirichlet k) where
    transition = breakPoint

instance KnownNat k => Transition Natural Source (Dirichlet k) where
    transition = breakPoint

-- Poisson Distribution --

instance Manifold Poisson where
    type Dimension Poisson = 1

instance Statistical Poisson where
    type SamplePoint Poisson = Int

instance ExponentialFamily Poisson where
    sufficientStatistic = Point . S.singleton . fromIntegral
    logBaseMeasure _ k = negate $ logFactorial k

type instance PotentialCoordinates Poisson = Natural

instance Legendre Poisson where
    potential = exp . S.head . coordinates

instance Transition Natural Mean Poisson where
    transition = Point . exp . coordinates

instance DuallyFlat Poisson where
    dualPotential (Point xs) =
        let eta = S.head xs
         in eta * log eta - eta

instance Transition Mean Natural Poisson where
    transition = Point . log . coordinates

instance Transition Source Natural Poisson where
    transition = transition . toMean

instance Transition Natural Source Poisson where
    transition = transition . toMean

instance Transition Source Mean Poisson where
    transition = breakPoint

instance Transition Mean Source Poisson where
    transition = breakPoint

instance (Transition c Source Poisson) => Generative c Poisson where
    samplePoint = samplePoisson . S.head . coordinates . toSource

instance AbsolutelyContinuous Source Poisson where
    densities (Point xs) ks = do
        k <- ks
        let lmda = S.head xs
        return $ lmda^k / factorial k * exp (-lmda)

instance AbsolutelyContinuous Mean Poisson where
    densities = densities . toSource

instance AbsolutelyContinuous Natural Poisson where
    logDensities = exponentialFamilyLogDensities

instance Transition Mean c Poisson => MaximumLikelihood c Poisson where
    mle = transition . averageSufficientStatistic

instance LogLikelihood Natural Poisson Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- Normal Distribution --

instance Manifold NormalMean where
    type Dimension NormalMean = 1

instance Manifold NormalShape where
    type Dimension NormalShape = 1

instance Statistical NormalMean where
    type SamplePoint NormalMean = Double

instance Statistical NormalShape where
    type SamplePoint NormalShape = Double

instance ExponentialFamily NormalMean where
    sufficientStatistic x = singleton x
    logBaseMeasure _ _ = -1/2 * log (2 * pi)

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
         in normal mu (sqrt vr)

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




-- Multivariate Normal --

scaleMatrix :: Double -> S.Matrix m n Double -> S.Matrix m n Double
scaleMatrix x = S.withMatrix (S.scale x)

instance (KnownNat n, KnownNat (Triangular n)) => Manifold (MultivariateNormal n) where
    type Dimension (MultivariateNormal n) = n + Triangular n

instance (KnownNat n, KnownNat (Triangular n)) => Statistical (MultivariateNormal n) where
    type SamplePoint (MultivariateNormal n) = S.Vector n Double

instance (KnownNat n, KnownNat (Triangular n))
  => AbsolutelyContinuous Source (MultivariateNormal n) where
      densities p xss = do
          let (mus,sgma) = splitMultivariateNormal p
              nrm = recip . sqrt . S.determinant $ scaleMatrix (2*pi) sgma
          xs <- xss
          let dff = xs - mus
              expval = S.dotProduct dff $ S.matrixVectorMultiply (S.pseudoInverse sgma) dff
          return $ nrm * exp (-expval / 2)

instance (KnownNat n, KnownNat (Triangular n), Transition c Source (MultivariateNormal n))
  => Generative c (MultivariateNormal n) where
    samplePoint = sampleMultivariateNormal . toSource

instance KnownNat n => Transition Source Natural (MultivariateNormal n) where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
            invsgma = S.pseudoInverse sgma
         in joinNaturalMultivariateNormal (S.matrixVectorMultiply invsgma mu) (scaleMatrix (-0.5) invsgma)

instance KnownNat n => Transition Natural Source (MultivariateNormal n) where
    transition p =
        let (nmu,nsgma) = splitNaturalMultivariateNormal p
            insgma = scaleMatrix (-0.5) $ S.pseudoInverse nsgma
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
            lndet = log $ S.determinant nsgma
         in -0.5 * ( 0.5 * S.dotProduct nmu (S.matrixVectorMultiply insgma nmu) + lndet )

instance (KnownNat n, KnownNat (Triangular n)) => Transition Natural Mean (MultivariateNormal n) where
    transition p =
        let (tmu,tsgma) = splitNaturalMultivariateNormal p
            itsgma = S.pseudoInverse tsgma
            mmu0 = S.matrixVectorMultiply itsgma tmu
            mmu = S.scale (-0.25) mmu0
            msgma1 = scaleMatrix (-0.25) $ S.outerProduct mmu0 mmu0
            msgma2 = scaleMatrix 0.5 itsgma
            msgma = msgma1 + msgma2
         in breakPoint $ joinMultivariateNormal0 mmu msgma

instance (KnownNat n, KnownNat (Triangular n)) => DuallyFlat (MultivariateNormal n) where
    dualPotential p =
        let sgma = snd . splitMultivariateNormal $ toSource p
            lndet = log . S.determinant $ scaleMatrix (2*pi*exp 1) sgma
         in -0.5 * lndet

instance (KnownNat n, KnownNat (Triangular n)) => Transition Mean Natural (MultivariateNormal n) where
    transition = breakPoint . toNatural . toSource

instance KnownNat n => Transition Source Mean (MultivariateNormal n) where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
            G.Matrix mumu = S.outerProduct mu mu
         in joinMultivariateNormal0 mu . G.Matrix $ S.add mumu (G.toVector sgma)

instance KnownNat n => Transition Mean Source (MultivariateNormal n) where
    transition p =
        let (mmu,msgma) = splitMultivariateNormal0 p
            G.Matrix mmumu = scaleMatrix (-1) $ S.outerProduct mmu mmu
         in joinMultivariateNormal mmu . G.Matrix $ S.add (G.toVector msgma) mmumu

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

-- VonMises --

instance Manifold VonMises where
    type Dimension VonMises = 2

instance Statistical VonMises where
    type SamplePoint VonMises = Double

instance Generative Source VonMises where
    samplePoint p@(Point cs) = do
        let (mu,kap0) = S.toPair cs
            kap = max kap0 1e-5
            tau = 1 + sqrt (1 + 4 * square kap)
            rho = (tau - sqrt (2*tau))/(2*kap)
            r = (1 + square rho) / (2 * rho)
        u1 <- uniform
        u2 <- uniform
        u3 <- uniform
        let z = cos (pi * u1)
            f = (1 + r * z)/(r + z)
            c = kap * (r - f)
        if log (c / u2) + 1 - c < 0
           then samplePoint p
           else return . toPi $ signum (u3 - 0.5) * acos f + mu

instance AbsolutelyContinuous Source VonMises where
    densities p xs = do
        let (mu,kp) = S.toPair $ coordinates p
        x <- xs
        return $ exp (kp * cos (x - mu)) / (2*pi * GSL.bessel_I0 kp)

instance LogLikelihood Natural VonMises Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

type instance PotentialCoordinates VonMises = Natural

instance Legendre VonMises where
    potential p =
        let kp = snd . S.toPair . coordinates $ toSource p
         in log $ GSL.bessel_I0 kp

instance Transition Natural Mean VonMises where
    transition p =
        let kp = snd . S.toPair . coordinates $ toSource p
         in breakPoint $ (GSL.bessel_I1 kp / (GSL.bessel_I0 kp * kp)) .> p

instance AbsolutelyContinuous Natural VonMises where
    logDensities = exponentialFamilyLogDensities

instance Generative Natural VonMises where
    samplePoint = samplePoint . toSource

instance ExponentialFamily VonMises where
    sufficientStatistic tht = Point $ S.doubleton (cos tht) (sin tht)
    logBaseMeasure _ _ = -log(2 * pi)

instance Transition Source Natural VonMises where
    transition (Point cs) =
        let (mu,kap) = S.toPair cs
         in Point $ S.doubleton (kap * cos mu) (kap * sin mu)

instance Transition Natural Source VonMises where
    transition (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in Point $ S.doubleton (toPi $ atan2 tht1 tht0) (sqrt $ square tht0 + square tht1)

instance Transition Source Mean VonMises where
    transition = toMean . toNatural


--- Location Scale ---

instance (Statistical l, Statistical s, SamplePoint l ~ SamplePoint s)
  => Statistical (LocationShape l s) where
    type SamplePoint (LocationShape l s) = SamplePoint l

instance (Manifold l, Manifold s) => Translation (LocationShape l s) l where
    (>+>) yz y' =
        let (y,z) = split yz
         in join (y + y') z
    anchor = fst . split

type instance PotentialCoordinates (LocationShape l s) = Natural

instance ( Statistical l, Statistical s , Product (LocationShape l s)
         , Storable (SamplePoint s), SamplePoint l ~ SamplePoint s
         , AbsolutelyContinuous c (LocationShape l s), KnownNat n)
  => AbsolutelyContinuous c (LocationShape (Replicated n l) (Replicated n s)) where
      logDensities lss xs =
          let (l,s) = split lss
              ls = splitReplicated l
              ss = splitReplicated s
              lss' :: c # Replicated n (LocationShape l s)
              lss' = joinReplicated $ S.zipWith join ls ss
           in logDensities lss' xs


instance (KnownNat n, Manifold l, Manifold s)
  => Translation (Replicated n (LocationShape l s)) (Replicated n l) where
      {-# INLINE (>+>) #-}
      (>+>) w z =
          let ws = splitReplicated w
              zs = splitReplicated z
           in joinReplicated $ S.zipWith (>+>) ws zs
      {-# INLINE anchor #-}
      anchor = mapReplicatedPoint anchor



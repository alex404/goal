{-# LANGUAGE UndecidableInstances #-}

-- | Various instances of statistical manifolds, with a focus on exponential families.
module Goal.Probability.Distributions
    ( -- * Exponential Families
      Bernoulli
    , Binomial
    , binomialTrials
    , Categorical
    , categoricalWeights
    , Dirichlet
    , Poisson
    , Normal
    , LogNormal
    , MeanNormal
    , StandardNormal
    , meanNormalVariance
    , meanNormalToNormal
    , VonMises
    , MultivariateNormal
    , joinMultivariateNormal
    , splitMultivariateNormal
    , multivariateNormalCorrelations
    , exponentialFamilyDensity
    ) where

-- Package --

import Goal.Core
import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily

import Goal.Geometry
import System.Random.MWC.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Numeric.GSL.Special.Bessel as GSL
import qualified Numeric.GSL.Special.Gamma as GSL
import qualified Numeric.GSL.Special.Psi as GSL

-- Uniform --

type StandardNormal = MeanNormal (1/1)

-- | A 'Uniform' distribution on a specified interval of the real line. This
-- distribution does not have interesting geometric properties, and does not
-- have coordinates.
--data Uniform mn mx

-- Bernoulli Distribution --

-- | The Bernoulli 'Family' with 'SampleSpace' 'Bernoulli' = 'Bool' (because why not).
data Bernoulli

-- Binomial Distribution --

-- | Models a number of coin flips, with a probability of tails given
-- by the parameter of the family.
data Binomial (n :: Nat)

-- | Returns the number of trials used to define this binomial distribution.
binomialTrials :: forall c n. KnownNat n => Point c (Binomial n) -> Int
{-# INLINE binomialTrials #-}
binomialTrials _ = natValInt (Proxy :: Proxy n)

-- | Returns the number of trials used to define this binomial distribution.
binomialSampleSpace :: forall n . KnownNat n => Proxy (Binomial n) -> Int
{-# INLINE binomialSampleSpace #-}
binomialSampleSpace _ = natValInt (Proxy :: Proxy n)

-- Categorical Distribution --

-- | A 'Categorical' distribution where the probability of the last category is
-- given by the normalization constraint.
data Categorical (n :: Nat)

-- | Takes a weighted list of elements representing a probability mass function, and
-- returns a sample from the Categorical distribution.
sampleCategorical :: KnownNat n => S.Vector n Double -> Random r Int
{-# INLINE sampleCategorical #-}
sampleCategorical ps = do
    let ps' = S.postscanl' (+) 0 ps
    p <- uniform
    let midx = (+1) . finiteInt <$> S.findIndex (> p) ps'
    return $ fromMaybe 0 midx

categoricalWeights
    :: Transition c Source (Categorical n)
    => c # Categorical n
    -> S.Vector (n+1) Double
{-# INLINE categoricalWeights #-}
categoricalWeights wghts0 =
    let wghts = coordinates $ toSource wghts0
     in S.cons (1-S.sum wghts) wghts

-- | A 'Dirichlet' manifold contains distributions over histogram weights.
data Dirichlet (k :: Nat)

-- Poisson Distribution --

-- | Returns a sample from a Poisson distribution with the given rate.
samplePoisson :: Double -> Random s Int
{-# INLINE samplePoisson #-}
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

-- Normal Distribution --

-- | The 'Manifold' of 'Normal' distributions. The standard coordinates are the
-- mean and the variance.
data Normal

-- Normal Distribution --

-- | The 'Manifold' of 'LogNormal' distributions.
data LogNormal

-- MeanNormal Distribution --

-- | The 'Manifold' of 'Normal' distributions with known variance. The standard
-- coordinate is simply the mean.
data MeanNormal v

-- | Returns the known variance of the given 'MeanNormal' distribution.
meanNormalVariance :: forall n d c . (KnownNat n, KnownNat d)
                   => Point c (MeanNormal (n/d)) -> Double
{-# INLINE meanNormalVariance #-}
meanNormalVariance _ = realToFrac $ ratVal (Proxy :: Proxy (n/d))

-- | Returns the known variance of the given 'MeanNormal' distribution.
meanNormalToNormal :: forall n d . (KnownNat n, KnownNat d)
                   => Source # MeanNormal (n/d) -> Source # Normal
{-# INLINE meanNormalToNormal #-}
meanNormalToNormal p = Point $ coordinates p S.++ S.singleton (meanNormalVariance p)


-- Multivariate Normal --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The standard coordinates are the
-- (vector) mean and the covariance matrix. When building a multivariate normal
-- distribution using e.g. 'fromList', the elements of the mean come first, and
-- then the elements of the covariance matrix in row major order.
data MultivariateNormal (n :: Nat)

splitMultivariateNormal0
    :: KnownNat n
    => c # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
{-# INLINE splitMultivariateNormal0 #-}
splitMultivariateNormal0 (Point xs) =
    let (mus,cvrs) = S.splitAt xs
     in (mus,S.fromLowerTriangular cvrs)

joinMultivariateNormal0
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> c # MultivariateNormal n
{-# INLINE joinMultivariateNormal0 #-}
joinMultivariateNormal0 mus sgma =
    Point $ mus S.++ S.lowerTriangular sgma

splitNaturalMultivariateNormal
    :: KnownNat n
    => Natural # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
{-# INLINE splitNaturalMultivariateNormal #-}
splitNaturalMultivariateNormal np =
    let (nmu,nsgma) = splitMultivariateNormal0 np
     in (nmu, addMatrix (scaleMatrix 0.5 nsgma) . scaleMatrix 0.5 . S.diagonalMatrix $ S.takeDiagonal nsgma)

joinNaturalMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Natural # MultivariateNormal n
{-# INLINE joinNaturalMultivariateNormal #-}
joinNaturalMultivariateNormal nmu nsgma =
    let nsgma' = addMatrix (scaleMatrix 2 nsgma) . scaleMatrix (-1) . S.diagonalMatrix $ S.takeDiagonal nsgma
     in joinMultivariateNormal0 nmu nsgma'

splitMultivariateNormal
    :: KnownNat n
    => Source # MultivariateNormal n
    -> (S.Vector n Double, S.Matrix n n Double)
{-# INLINE splitMultivariateNormal #-}
splitMultivariateNormal = splitMultivariateNormal0

joinMultivariateNormal
    :: KnownNat n
    => S.Vector n Double
    -> S.Matrix n n Double
    -> Source # MultivariateNormal n
{-# INLINE joinMultivariateNormal #-}
joinMultivariateNormal mus sgma =
    Point $ mus S.++ S.lowerTriangular sgma

multivariateNormalCorrelations
    :: KnownNat k
    => Source # MultivariateNormal k
    -> S.Matrix k k Double
{-# INLINE multivariateNormalCorrelations #-}
multivariateNormalCorrelations mnrm =
    let cvrs = snd $ splitMultivariateNormal mnrm
        sds = S.map sqrt $ S.takeDiagonal cvrs
        sdmtx = S.outerProduct sds sds
     in G.Matrix $ S.zipWith (/) (G.toVector cvrs) (G.toVector sdmtx)

multivariateNormalBaseMeasure :: forall n . (KnownNat n)
                               => Proxy (MultivariateNormal n) -> S.Vector n Double -> Double
multivariateNormalBaseMeasure _ _ =
    let n = natValInt (Proxy :: Proxy n)
     in pi**(-fromIntegral n/2)

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


---- | The unnormalized log-density of an arbitrary exponential family distribution.
--unnormalizedLogDensity :: forall x s . ExponentialFamily x s => Natural # x -> s -> Double
--unnormalizedLogDensity p x =
--    p <.> sufficientStatistic x  + log (baseMeasure (Proxy @ x) x)


binomialBaseMeasure0 :: (KnownNat n) => Proxy n -> Proxy (Binomial n) -> Int -> Double
{-# INLINE binomialBaseMeasure0 #-}
binomialBaseMeasure0 prxyn _ = choose (natValInt prxyn)

meanNormalBaseMeasure0 :: (KnownNat n, KnownNat d) => Proxy (n/d) -> Proxy (MeanNormal (n/d)) -> Double -> Double
{-# INLINE meanNormalBaseMeasure0 #-}
meanNormalBaseMeasure0 prxyr _ x =
    let vr = realToFrac $ ratVal prxyr
     in (exp . negate $ 0.5 * square x / vr) / sqrt (2*pi*vr)

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
    baseMeasure _ _ = 1
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic True = Point $ S.singleton 1
    sufficientStatistic False = Point $ S.singleton 0

instance Legendre Bernoulli where
    type PotentialCoordinates Bernoulli = Natural
    {-# INLINE potential #-}
    potential p = log $ 1 + exp (S.head $ coordinates p)

--instance {-# OVERLAPS #-} KnownNat k => Legendre (Replicated k Bernoulli) where
--    {-# INLINE potential #-}
--    potential p = S.sum . S.map (log . (1 +) .  exp) $ coordinates p

instance Transition Natural Mean Bernoulli where
    {-# INLINE transition #-}
    transition = Point . S.map logistic . coordinates

instance DuallyFlat Bernoulli where
    {-# INLINE dualPotential #-}
    dualPotential p =
        let eta = S.head $ coordinates p
         in logit eta * eta - log (1 / (1 - eta))

instance Transition Mean Natural Bernoulli where
    {-# INLINE transition #-}
    transition = Point . S.map logit . coordinates

instance Riemannian Natural Bernoulli where
    {-# INLINE metric #-}
    metric p =
        let stht = logistic . S.head $ coordinates p
         in Point . S.singleton $ stht * (1-stht)
    {-# INLINE flat #-}
    flat p p' =
        let stht = logistic . S.head $ coordinates p
         in breakPoint $ (stht * (1-stht)) .> p'

instance {-# OVERLAPS #-} KnownNat k => Riemannian Natural (Replicated k Bernoulli) where
    {-# INLINE metric #-}
    metric = error "Do not call metric on a replicated manifold"
    {-# INLINE flat #-}
    flat p p' =
        let sthts = S.map ((\stht -> stht * (1-stht)) . logistic) $ coordinates p
            dp = S.zipWith (*) sthts $ coordinates p'
         in Point dp

instance {-# OVERLAPS #-} KnownNat k => Riemannian Mean (Replicated k Bernoulli) where
    {-# INLINE metric #-}
    metric = error "Do not call metric on a replicated manifold"
    {-# INLINE sharp #-}
    sharp p dp =
        let sthts' = S.map (\stht -> stht * (1-stht)) $ coordinates p
            p' = S.zipWith (*) sthts' $ coordinates dp
         in Point p'

instance Transition Source Mean Bernoulli where
    {-# INLINE transition #-}
    transition = breakPoint

instance Transition Mean Source Bernoulli where
    {-# INLINE transition #-}
    transition = breakPoint

instance Transition Source Natural Bernoulli where
    {-# INLINE transition #-}
    transition = transition . toMean

instance Transition Natural Source Bernoulli where
    {-# INLINE transition #-}
    transition = transition . toMean

instance (Transition c Source Bernoulli) => Generative c Bernoulli where
    {-# INLINE samplePoint #-}
    samplePoint = bernoulli . S.head . coordinates . toSource

instance Transition Mean c Bernoulli => MaximumLikelihood c Bernoulli where
    mle = transition . sufficientStatisticT

instance LogLikelihood Natural Bernoulli Bool where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance AbsolutelyContinuous Source Bernoulli where
    density (Point p) True = S.head p
    density (Point p) False = 1 - S.head p

instance AbsolutelyContinuous Mean Bernoulli where
    density = density . toSource

instance AbsolutelyContinuous Natural Bernoulli where
    density = exponentialFamilyDensity

-- Binomial Distribution --

instance KnownNat n => Manifold (Binomial n) where
    type Dimension (Binomial n) = 1

instance KnownNat n => Statistical (Binomial n) where
    type SamplePoint (Binomial n) = Int

instance KnownNat n => Discrete (Binomial n) where
    type Cardinality (Binomial n) = n + 1
    sampleSpace prx = [0..binomialSampleSpace prx]

instance KnownNat n => ExponentialFamily (Binomial n) where
    baseMeasure = binomialBaseMeasure0 Proxy
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic = Point . S.singleton . fromIntegral

instance KnownNat n => Legendre (Binomial n) where
    type PotentialCoordinates (Binomial n) = Natural
    {-# INLINE potential #-}
    potential p =
        let n = fromIntegral $ binomialTrials p
            tht = S.head $ coordinates p
         in n * log (1 + exp tht)

instance KnownNat n => Transition Natural Mean (Binomial n) where
    {-# INLINE transition #-}
    transition p =
        let n = fromIntegral $ binomialTrials p
         in Point . S.singleton $ n * logistic (S.head $ coordinates p)

instance KnownNat n => DuallyFlat (Binomial n) where
    {-# INLINE dualPotential #-}
    dualPotential p =
        let n = fromIntegral $ binomialTrials p
            eta = S.head $ coordinates p
        in eta * log (eta / (n - eta)) - n * log (n / (n - eta))

instance KnownNat n => Transition Mean Natural (Binomial n) where
    {-# INLINE transition #-}
    transition p =
        let n = fromIntegral $ binomialTrials p
            eta = S.head $ coordinates p
         in Point . S.singleton . log $ eta / (n - eta)

instance KnownNat n => Transition Source Natural (Binomial n) where
    {-# INLINE transition #-}
    transition = transition . toMean

instance KnownNat n => Transition Natural Source (Binomial n) where
    {-# INLINE transition #-}
    transition = transition . toMean

instance KnownNat n => Transition Source Mean (Binomial n) where
    {-# INLINE transition #-}
    transition p =
        let n = fromIntegral $ binomialTrials p
         in breakPoint $ n .> p

instance KnownNat n => Transition Mean Source (Binomial n) where
    {-# INLINE transition #-}
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
    density p k =
        let n = binomialTrials p
            c = S.head $ coordinates p
         in choose n k * c^k * (1 - c)^(n-k)

instance KnownNat n => AbsolutelyContinuous Mean (Binomial n) where
    density = density . toSource

instance KnownNat n => AbsolutelyContinuous Natural (Binomial n) where
    density = exponentialFamilyDensity

instance (KnownNat n, Transition Mean c (Binomial n)) => MaximumLikelihood c (Binomial n) where
    mle = transition . sufficientStatisticT

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
    baseMeasure _ _ = 1
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic e = Point $ S.generate (\i -> if finiteInt i == (fromEnum e-1) then 1 else 0)

instance KnownNat n => Legendre (Categorical n) where
    type (PotentialCoordinates (Categorical n)) = Natural
    {-# INLINE potential #-}
    --potential (Point cs) = log $ 1 + S.sum (S.map exp cs)
    potential = logSumExp . B.cons 0 . boxCoordinates

instance KnownNat n => Transition Natural Mean (Categorical n) where
    {-# INLINE transition #-}
    transition p =
        let exps = S.map exp $ coordinates p
            nrm = 1 + S.sum exps
         in nrm /> Point exps

instance KnownNat n => DuallyFlat (Categorical n) where
    {-# INLINE dualPotential #-}
    dualPotential (Point cs) =
        let sc = 1 - S.sum cs
         in S.sum (S.map entropyFun cs) + entropyFun sc
        where entropyFun 0 = 0
              entropyFun x = x * log x

instance KnownNat n => Transition Mean Natural (Categorical n) where
    {-# INLINE transition #-}
    transition (Point xs) =
        let nrm = 1 - S.sum xs
         in  Point . log $ S.map (/nrm) xs

instance Transition Source Mean (Categorical n) where
    {-# INLINE transition #-}
    transition = breakPoint

instance Transition Mean Source (Categorical n) where
    {-# INLINE transition #-}
    transition = breakPoint

instance KnownNat n => Transition Source Natural (Categorical n) where
    {-# INLINE transition #-}
    transition = transition . toMean

instance KnownNat n => Transition Natural Source (Categorical n) where
    {-# INLINE transition #-}
    transition = transition . toMean

instance (KnownNat n, Transition c Source (Categorical n))
  => Generative c (Categorical n) where
    {-# INLINE samplePoint #-}
    samplePoint p0 =
        let p = toSource p0
         in sampleCategorical $ coordinates p

instance (KnownNat n, Transition Mean c (Categorical n))
  => MaximumLikelihood c (Categorical n) where
    mle = transition . sufficientStatisticT

instance KnownNat n => LogLikelihood Natural (Categorical n) Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


instance KnownNat n => AbsolutelyContinuous Source (Categorical n) where
    density (Point ps) e =
        let ek = fromEnum e
         in if ek == 0
               then 1 - S.sum ps
               else S.unsafeIndex ps $ ek - 1

instance KnownNat n => AbsolutelyContinuous Mean (Categorical n) where
    density = density . toSource

instance KnownNat n => AbsolutelyContinuous Natural (Categorical n) where
    density = exponentialFamilyDensity

-- Dirichlet Distribution --

instance KnownNat k => Manifold (Dirichlet k) where
    type Dimension (Dirichlet k) = k

instance KnownNat k => Statistical (Dirichlet k) where
    type SamplePoint (Dirichlet k) = S.Vector k Double

instance (KnownNat k, Transition c Source (Dirichlet k))
  => Generative c (Dirichlet k) where
    {-# INLINE samplePoint #-}
    samplePoint p0 = do
        let alphs = boxCoordinates $ toSource p0
        G.convert <$> dirichlet alphs

instance KnownNat k => ExponentialFamily (Dirichlet k) where
    {-# INLINE baseMeasure #-}
    baseMeasure _ = recip . S.product
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = Point $ S.map log xs

logMultiBeta :: KnownNat k => S.Vector k Double -> Double
{-# INLINE logMultiBeta #-}
logMultiBeta alphs =
    S.sum (S.map GSL.lngamma alphs) - GSL.lngamma (S.sum alphs)

logMultiBetaDifferential :: KnownNat k => S.Vector k Double -> S.Vector k Double
{-# INLINE logMultiBetaDifferential #-}
logMultiBetaDifferential alphs =
    S.map (subtract (GSL.psi $ S.sum alphs) . GSL.psi) alphs

instance KnownNat k => Legendre (Dirichlet k) where
    type PotentialCoordinates (Dirichlet k) = Natural
    {-# INLINE potential #-}
    potential = logMultiBeta . coordinates

instance KnownNat k => Transition Natural Mean (Dirichlet k) where
    {-# INLINE transition #-}
    transition = Point . logMultiBetaDifferential . coordinates

instance KnownNat k => AbsolutelyContinuous Source (Dirichlet k) where
    density p xs =
        let alphs = coordinates p
            prds = S.product $ S.zipWith (**) xs $ S.map (subtract 1) alphs
         in prds / exp (logMultiBeta alphs)

instance KnownNat k => AbsolutelyContinuous Natural (Dirichlet k) where
    density = exponentialFamilyDensity

instance KnownNat k => LogLikelihood Natural (Dirichlet k) (S.Vector k Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance KnownNat k => Transition Source Natural (Dirichlet k) where
    {-# INLINE transition #-}
    transition = breakPoint

instance KnownNat k => Transition Natural Source (Dirichlet k) where
    {-# INLINE transition #-}
    transition = breakPoint

-- Poisson Distribution --

instance Manifold Poisson where
    type Dimension Poisson = 1

instance Statistical Poisson where
    type SamplePoint Poisson = Int

instance ExponentialFamily Poisson where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic = Point . S.singleton . fromIntegral
    baseMeasure _ k = recip $ factorial k

instance Legendre Poisson where
    type PotentialCoordinates Poisson = Natural
    {-# INLINE potential #-}
    potential = exp . S.head . coordinates

instance Transition Natural Mean Poisson where
    {-# INLINE transition #-}
    transition = Point . exp . coordinates

instance DuallyFlat Poisson where
    {-# INLINE dualPotential #-}
    dualPotential (Point xs) =
        let eta = S.head xs
         in eta * log eta - eta

instance Transition Mean Natural Poisson where
    {-# INLINE transition #-}
    transition = Point . log . coordinates

instance Transition Source Natural Poisson where
    {-# INLINE transition #-}
    transition = transition . toMean

instance Transition Natural Source Poisson where
    {-# INLINE transition #-}
    transition = transition . toMean

instance Transition Source Mean Poisson where
    {-# INLINE transition #-}
    transition = breakPoint

instance Transition Mean Source Poisson where
    {-# INLINE transition #-}
    transition = breakPoint

instance (Transition c Source Poisson) => Generative c Poisson where
    {-# INLINE samplePoint #-}
    samplePoint = samplePoisson . S.head . coordinates . toSource

instance AbsolutelyContinuous Source Poisson where
    density (Point xs) k =
        let lmda = S.head xs
         in  lmda^k / factorial k * exp (-lmda)

instance AbsolutelyContinuous Mean Poisson where
    density = density . toSource

instance AbsolutelyContinuous Natural Poisson where
    density = exponentialFamilyDensity

instance Transition Mean c Poisson => MaximumLikelihood c Poisson where
    mle = transition . sufficientStatisticT

instance LogLikelihood Natural Poisson Int where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- Normal Distribution --

instance Manifold Normal where
    type Dimension Normal = 2

instance Statistical Normal where
    type SamplePoint Normal = Double

instance ExponentialFamily Normal where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic x =
         Point . S.doubleton x $ x**2
    {-# INLINE baseMeasure #-}
    baseMeasure _ _ = recip . sqrt $ 2 * pi

instance Legendre Normal where
    type PotentialCoordinates Normal = Natural
    {-# INLINE potential #-}
    potential (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in -(square tht0 / (4*tht1)) - 0.5 * log(-2*tht1)

instance Transition Natural Mean Normal where
    {-# INLINE transition #-}
    transition p =
        let (tht0,tht1) = S.toPair $ coordinates p
            dv = tht0/tht1
         in Point $ S.doubleton (-0.5*dv) (0.25 * square dv - 0.5/tht1)

instance DuallyFlat Normal where
    {-# INLINE dualPotential #-}
    dualPotential (Point cs) =
        let (eta0,eta1) = S.toPair cs
         in -0.5 * log(eta1 - square eta0) - 1/2

instance Transition Mean Natural Normal where
    {-# INLINE transition #-}
    transition p =
        let (eta0,eta1) = S.toPair $ coordinates p
            dff = eta1 - square eta0
         in Point $ S.doubleton (eta0 / dff) (-0.5 / dff)

instance Riemannian Natural Normal where
    {-# INLINE metric #-}
    metric p =
        let (tht0,tht1) = S.toPair $ coordinates p
            d00 = -1/(2*tht1)
            d01 = tht0/(2*square tht1)
            d11 = 0.5*(1/square tht1 - square tht0 / (tht1^(3 :: Int)))
         in Point $ S.doubleton d00 d01 S.++ S.doubleton d01 d11

instance Riemannian Mean Normal where
    {-# INLINE metric #-}
    metric p =
        let (eta0,eta1) = S.toPair $ coordinates p
            eta02 = square eta0
            dff2 = square $ eta1 - eta02
            d00 = (dff2 + 2 * eta02) / dff2
            d01 = -eta0 / dff2
            d11 = 0.5 / dff2
         in Point $ S.doubleton d00 d01 S.++ S.doubleton d01 d11

-- instance Riemannian Source Normal where
--     {-# INLINE metric #-}
--     metric p =
--         let (_,vr) = S.toPair $ coordinates p
--          in Point $ S.doubleton (recip vr) 0 S.++ S.doubleton 0 (recip $ 2*square vr)

instance Transition Source Mean Normal where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (mu,vr) = S.toPair cs
         in Point . S.doubleton mu $ vr + square mu

instance Transition Mean Source Normal where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (eta0,eta1) = S.toPair cs
         in Point . S.doubleton eta0 $ eta1 - square eta0

instance Transition Source Natural Normal where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (mu,vr) = S.toPair cs
         in Point $ S.doubleton (mu / vr) (negate . recip $ 2 * vr)

instance Transition Natural Source Normal where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in Point $ S.doubleton (-0.5 * tht0 / tht1) (negate . recip $ 2 * tht1)

instance (Transition c Source Normal) => Generative c Normal where
    {-# INLINE samplePoint #-}
    samplePoint p =
        let (Point cs) = toSource p
            (mu,vr) = S.toPair cs
         in normal mu (sqrt vr)

instance AbsolutelyContinuous Source Normal where
    density (Point cs) x =
        let (mu,vr) = S.toPair cs
         in recip (sqrt $ vr*2*pi) * exp (negate $ (x - mu) ** 2 / (2*vr))

instance AbsolutelyContinuous Mean Normal where
    density = density . toSource

instance AbsolutelyContinuous Natural Normal where
    density = exponentialFamilyDensity

instance Transition Mean c Normal => MaximumLikelihood c Normal where
    mle = transition . sufficientStatisticT

instance LogLikelihood Natural Normal Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- LogNormal Distribution --

instance Manifold LogNormal where
    type Dimension LogNormal = 2

instance Statistical LogNormal where
    type SamplePoint LogNormal = Double

instance ExponentialFamily LogNormal where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic x =
         Point . S.doubleton (log x) $ log x**2
    {-# INLINE baseMeasure #-}
    baseMeasure _ x = recip $ x * sqrt (2 * pi)

toNaturalNormal :: Natural # LogNormal -> Natural # Normal
toNaturalNormal = breakPoint

toMeanNormal :: Mean # LogNormal -> Mean # Normal
toMeanNormal = breakPoint

toSourceNormal :: Source # LogNormal -> Source # Normal
toSourceNormal = breakPoint


instance Legendre LogNormal where
    type PotentialCoordinates LogNormal = Natural
    {-# INLINE potential #-}
    potential = potential . toNaturalNormal

instance Transition Natural Mean LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toMean . toNaturalNormal

instance DuallyFlat LogNormal where
    {-# INLINE dualPotential #-}
    dualPotential = dualPotential . toMeanNormal

instance Transition Mean Natural LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toNatural . toMeanNormal

instance Riemannian Natural LogNormal where
    {-# INLINE metric #-}
    metric = breakPoint . metric . toNaturalNormal

instance Riemannian Mean LogNormal where
    {-# INLINE metric #-}
    metric = breakPoint . metric . toMeanNormal

--instance Riemannian Source LogNormal where
--    {-# INLINE metric #-}
--    metric = breakPoint . metric . toSourceNormal

instance Transition Source Mean LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toMean . toSourceNormal

instance Transition Mean Source LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toSource . toMeanNormal

instance Transition Source Natural LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toNatural . toSourceNormal

instance Transition Natural Source LogNormal where
    {-# INLINE transition #-}
    transition = breakPoint . toSource . toNaturalNormal

instance (Transition c Source LogNormal) => Generative c LogNormal where
    {-# INLINE samplePoint #-}
    samplePoint p = do
        let nrm = toSourceNormal $ toSource p
        exp <$> samplePoint nrm

instance AbsolutelyContinuous Source LogNormal where
    density (Point cs) x =
        let (mu,vr) = S.toPair cs
         in recip (x * sqrt (vr*2*pi)) * exp (negate $ (log x - mu) ** 2 / (2*vr))

instance AbsolutelyContinuous Mean LogNormal where
    density = density . toSource

instance AbsolutelyContinuous Natural LogNormal where
    density = exponentialFamilyDensity

instance Transition Mean c LogNormal => MaximumLikelihood c LogNormal where
    mle = transition . sufficientStatisticT

instance LogLikelihood Natural LogNormal Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential



-- MeanNormal Distribution --

instance Manifold (MeanNormal v) where
    type Dimension (MeanNormal v) = 1

instance Statistical (MeanNormal v) where
    type SamplePoint (MeanNormal v) = Double

instance (KnownNat n, KnownNat d) => ExponentialFamily (MeanNormal (n / d)) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic = Point . S.singleton
    baseMeasure = meanNormalBaseMeasure0 Proxy

instance (KnownNat n, KnownNat d) => Legendre (MeanNormal (n/d)) where
    type PotentialCoordinates (MeanNormal (n/d)) = Natural
    {-# INLINE potential #-}
    potential p =
        let vr = meanNormalVariance p
            mu = S.head $ coordinates p
         in 0.5 * vr * square mu

instance (KnownNat n, KnownNat d) => Transition Natural Mean (MeanNormal (n/d)) where
    {-# INLINE transition #-}
    transition p =
        let vr = meanNormalVariance p
         in Point . S.singleton $ vr * S.head (coordinates p)

instance (KnownNat n, KnownNat d) => DuallyFlat (MeanNormal (n/d)) where
    {-# INLINE dualPotential #-}
    dualPotential p =
        let vr = meanNormalVariance p
            mu = S.head $ coordinates p
         in 0.5 / vr * square mu

instance (KnownNat n, KnownNat d) => Transition Mean Natural (MeanNormal (n/d)) where
    {-# INLINE transition #-}
    transition p =
        let vr = meanNormalVariance p
         in Point . S.singleton $ S.head (coordinates p) / vr

instance Transition Source Mean (MeanNormal v) where
    {-# INLINE transition #-}
    transition = breakPoint

instance Transition Mean Source (MeanNormal v) where
    {-# INLINE transition #-}
    transition = breakPoint

instance (KnownNat n, KnownNat d) => Transition Source Natural (MeanNormal (n/d)) where
    {-# INLINE transition #-}
    transition = transition . toMean

instance (KnownNat n, KnownNat d) => Transition Natural Source (MeanNormal (n/d)) where
    {-# INLINE transition #-}
    transition = toSource . toMean

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Source (MeanNormal (n/d)) where
    density p =
        let vr = meanNormalVariance p
            mu = S.head $ coordinates p
            nrm :: Double -> Double -> Point Source Normal
            nrm x y = Point $ S.doubleton x y
         in density $ nrm mu vr

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Mean (MeanNormal (n/d)) where
    density = density . toSource

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Natural (MeanNormal (n/d)) where
    density = exponentialFamilyDensity

instance (KnownNat n, KnownNat d, Transition Mean c (MeanNormal (n/d)))
  => MaximumLikelihood c (MeanNormal (n/d)) where
    mle = transition . sufficientStatisticT

instance (KnownNat n, KnownNat d, Transition c Source (MeanNormal (n/d)))
  => Generative c (MeanNormal (n/d)) where
    samplePoint p =
        let (Point cs) = toSource p
            mu = S.head cs
            vr = meanNormalVariance p
         in normal mu (sqrt vr)

instance (KnownNat n, KnownNat d) => LogLikelihood Natural (MeanNormal (n/d)) Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


-- Multivariate Normal --

scaleMatrix :: Double -> S.Matrix m n Double -> S.Matrix m n Double
scaleMatrix x = S.withMatrix (S.scale x)

addMatrix :: S.Matrix m n Double -> S.Matrix m n Double -> S.Matrix m n Double
addMatrix (G.Matrix xs) (G.Matrix ys) = G.Matrix $ S.add xs ys

instance (KnownNat n, KnownNat (S.Triangular n)) => Manifold (MultivariateNormal n) where
    type Dimension (MultivariateNormal n) = n + S.Triangular n

instance (KnownNat n, KnownNat (S.Triangular n)) => Statistical (MultivariateNormal n) where
    type SamplePoint (MultivariateNormal n) = S.Vector n Double

instance (KnownNat n, KnownNat (S.Triangular n))
  => AbsolutelyContinuous Source (MultivariateNormal n) where
    density p xs =
        let (mus,sgma) = splitMultivariateNormal p
            nrm = recip . sqrt . S.determinant $ scaleMatrix (2*pi) sgma
            dff = S.add xs (S.scale (-1) mus)
            expval = S.dotProduct dff $ S.matrixVectorMultiply (S.inverse sgma) dff
         in nrm * exp (-expval / 2)

instance (KnownNat n, KnownNat (S.Triangular n), Transition c Source (MultivariateNormal n))
  => Generative c (MultivariateNormal n) where
    samplePoint = sampleMultivariateNormal . toSource

instance KnownNat n => Transition Source Natural (MultivariateNormal n) where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
            invsgma = S.inverse sgma
         in joinNaturalMultivariateNormal (S.matrixVectorMultiply invsgma mu) (scaleMatrix (-0.5) invsgma)

instance KnownNat n => Transition Natural Source (MultivariateNormal n) where
    transition p =
        let (nmu,nsgma) = splitNaturalMultivariateNormal p
            insgma = scaleMatrix (-0.5) $ S.inverse nsgma
         in joinMultivariateNormal (S.matrixVectorMultiply insgma nmu) insgma

instance KnownNat n => LogLikelihood Natural (MultivariateNormal n) (S.Vector n Double) where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential


instance (KnownNat n, KnownNat (S.Triangular n)) => ExponentialFamily (MultivariateNormal n) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs = Point $ xs S.++ S.lowerTriangular (S.outerProduct xs xs)
    baseMeasure = multivariateNormalBaseMeasure

instance (KnownNat n, KnownNat (S.Triangular n)) => Legendre (MultivariateNormal n) where
    type PotentialCoordinates (MultivariateNormal n) = Natural
    {-# INLINE potential #-}
    potential p =
        let (nmu,nsgma) = splitNaturalMultivariateNormal p
            (insgma,lndet,_) = S.inverseLogDeterminant nsgma
         in -0.5 * ( 0.5 * S.dotProduct nmu (S.matrixVectorMultiply insgma nmu) + lndet )

instance (KnownNat n, KnownNat (S.Triangular n)) => Transition Natural Mean (MultivariateNormal n) where
    {-# INLINE transition #-}
    transition p =
        let (tmu,tsgma) = splitNaturalMultivariateNormal p
            itsgma = S.inverse tsgma
            mmu0 = S.matrixVectorMultiply itsgma tmu
            mmu = S.scale (-0.25) mmu0
            msgma1 = scaleMatrix (-0.25) $ S.outerProduct mmu0 mmu0
            msgma2 = scaleMatrix 0.5 itsgma
            msgma = addMatrix msgma1 msgma2
         in breakPoint $ joinMultivariateNormal0 mmu msgma

instance (KnownNat n, KnownNat (S.Triangular n)) => DuallyFlat (MultivariateNormal n) where
    {-# INLINE dualPotential #-}
    dualPotential p =
        let sgma = snd . splitMultivariateNormal $ toSource p
            (_,lndet,_) = S.inverseLogDeterminant $ scaleMatrix (2*pi*exp 1) sgma
         in -0.5 * lndet

instance (KnownNat n, KnownNat (S.Triangular n)) => Transition Mean Natural (MultivariateNormal n) where
    {-# INLINE transition #-}
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

instance (KnownNat n, KnownNat (S.Triangular n)) => AbsolutelyContinuous Natural (MultivariateNormal n) where
    density = exponentialFamilyDensity

instance (KnownNat n, Transition Mean c (MultivariateNormal n))
  => MaximumLikelihood c (MultivariateNormal n) where
    mle = transition . sufficientStatisticT


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
    {-# INLINE samplePoint #-}
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
    density p x =
        let (mu,kp) = S.toPair $ coordinates p
         in exp (kp * cos (x - mu)) / (2*pi * GSL.bessel_I0 kp)

instance LogLikelihood Natural VonMises Double where
    logLikelihood = exponentialFamilyLogLikelihood
    logLikelihoodDifferential = exponentialFamilyLogLikelihoodDifferential

instance Legendre VonMises where
    type PotentialCoordinates VonMises = Natural
    {-# INLINE potential #-}
    potential p =
        let kp = snd . S.toPair . coordinates $ toSource p
         in log $ GSL.bessel_I0 kp

instance Transition Natural Mean VonMises where
    {-# INLINE transition #-}
    transition p =
        let kp = snd . S.toPair . coordinates $ toSource p
         in breakPoint $ (GSL.bessel_I1 kp / (GSL.bessel_I0 kp * kp)) .> p

instance AbsolutelyContinuous Natural VonMises where
    density = exponentialFamilyDensity

instance Generative Natural VonMises where
    samplePoint = samplePoint . toSource

instance ExponentialFamily VonMises where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic tht = Point $ S.doubleton (cos tht) (sin tht)
    {-# INLINE baseMeasure #-}
    baseMeasure _ _ = recip $ 2 * pi

instance Transition Source Natural VonMises where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (mu,kap) = S.toPair cs
         in Point $ S.doubleton (kap * cos mu) (kap * sin mu)

instance Transition Natural Source VonMises where
    {-# INLINE transition #-}
    transition (Point cs) =
        let (tht0,tht1) = S.toPair cs
         in Point $ S.doubleton (toPi $ atan2 tht1 tht0) (sqrt $ square tht0 + square tht1)

instance Transition Source Mean VonMises where
    {-# INLINE transition #-}
    transition = toMean . toNatural

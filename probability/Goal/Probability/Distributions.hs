{-# LANGUAGE UndecidableInstances #-}

-- | Various instances of statistical manifolds, with a focus on exponential families.
module Goal.Probability.Distributions
    ( -- * Exponential Families
      Bernoulli
    , Binomial
    , binomialTrials
    , Categorical
    , categories
    , Poisson
    , Normal
    , MeanNormal
    , meanNormalVariance
    , VonMises
    ) where

-- Package --

import Goal.Core
import Goal.Probability.Statistical
import Goal.Probability.ExponentialFamily

import Goal.Geometry
import System.Random.MWC.Probability


-- Uniform --

-- | A 'Uniform' distribution on a specified interval of the real line. This
-- distribution does not have interesting geometric properties, and does not
-- have coordinates.
data Uniform

-- Bernoulli Distribution --

-- | The Bernoulli 'Family' with 'SampleSpace' 'Bernoulli' = 'Bool' (because why not).
data Bernoulli

-- Binomial Distribution --

-- | Models a number of coin flips, with a probability of tails given
-- by the parameter of the family.
data Binomial (n :: Nat)

binomialTrials :: (KnownNat n, RealFloat x) => Point c (Binomial n) x -> Int
binomialTrials = binomialTrials0 Proxy

categories :: (1 <= n, KnownNat n, Enum e) => Point c (Categorical e n) x -> Vector n e
categories = categories0 Proxy

-- Categorical Distribution --

-- | A 'Categorical' distribution where the probability of the last category is
-- given by the normalization constraint.
data Categorical e (n :: Nat)

-- | Takes a weighted list of elements representing a probability mass function, and
-- returns a sample from the Categorical distribution.
generateCategorical :: (Dense x, KnownNat n, 1 <= n) => Vector n a -> Vector (n-1) x -> Random s a
generateCategorical as ps = do
    let (as',an) = splitV as
        aps' = zipV as' . scanl1V' (+) $ realToFrac <$> ps
    p <- uniform
    let ma = fst <$> find (\(_,p') -> (p :: Double) < p') aps'
    return $ fromMaybe (headV an) ma

-- Curved Categorical Distribution --

-- | A 'CurvedCategorical' distribution is a 'Categorical' distribution where
-- each probability is explicitly represented.
data CurvedCategorical s

-- Poisson Distribution --

-- | Returns a sample from a Poisson distribution with the given rate.
generatePoisson :: Dense x => x -> Random s Int
generatePoisson lmda = uniform >>= renew 0
    where l = realToFrac $ exp (-lmda)
          renew k p
            | p <= l = return k
            | otherwise = do
                u <- uniform
                renew (k+1) (p*(u :: Double))

-- | The 'Manifold' of 'Poisson' distributions. The 'Source' coordinate is the
-- rate of the Poisson distribution.
data Poisson

-- Normal Distribution --

-- | The 'Manifold' of 'Normal' distributions. The standard coordinates are the
-- mean and the variance.
data Normal

-- MeanNormal Distribution --

-- | The 'Manifold' of 'Normal' distributions with known variance. The standard
-- coordinate is simply the mean.
data MeanNormal v

meanNormalVariance :: (KnownNat n, KnownNat d, RealFloat x) => Point c (MeanNormal (n/d)) x -> Rational
meanNormalVariance = meanNormalVariance0 Proxy
-- Multivariate Normal --

-- | The 'Manifold' of 'MultivariateNormal' distributions. The standard coordinates are the
-- (vector) mean and the covariance matrix. When building a multivariate normal
-- distribution using e.g. 'fromList', the elements of the mean come first, and
-- then the elements of the covariance matrix in row major order.
data MultivariateNormal (n :: Nat)

splitMultivariateNormal :: KnownNat n => Point c (MultivariateNormal n) x -> (Vector n x, Matrix n n x)
splitMultivariateNormal (Point xs) =
    let (mus,cvrs) = splitV xs
     in (mus,Matrix cvrs)

{-
-- | Samples from a multivariate Normal.
generateMultivariateNormal :: C.Vector Double -> M.Matrix Double -> RandST s (C.Vector Double)
generateMultivariateNormal mus rtsgma = do
    nrms <- C.replicateM n $ normal 0 1
    return $ mus + (M.#>) rtsgma nrms
    where n = C.length mus

-- | Generates a multivariateNormal by way of a covariance matrix i.e. by taking
-- the square root.
joinMultivariateNormal :: C.Vector Double -> M.Matrix Double -> c :#: MultivariateNormal
joinMultivariateNormal mus sgma =
    fromCoordinates (MultivariateNormal $ C.length mus) $ mus C.++ M.flatten sgma

     -}

-- von Mises --

-- | The 'Manifold' of 'VonMises' distributions. The 'Source' coordinates are
-- the mean and concentration.
data VonMises




--- Internal ---

binomialTrials0 :: (KnownNat n, RealFloat x) => Proxy n -> Point c (Binomial n) x -> Int
binomialTrials0 prxy _ = natValInt prxy

meanNormalVariance0 :: (KnownNat n, KnownNat d, RealFloat x) => Proxy (n/d) -> Point c (MeanNormal (n/d)) x -> Rational
meanNormalVariance0 prxy _ = ratVal prxy

categories0 :: (1 <= n, KnownNat n, Enum e) => Proxy (Categorical e n) -> Point c (Categorical e n) x -> Vector n e
categories0 prxy _ = sampleSpace prxy

binomialBaseMeasure0 :: (KnownNat n, Dense x) => Proxy n -> Proxy (Binomial n) -> Sample (Binomial n) -> x
binomialBaseMeasure0 prxyn _ = realToFrac . choose (natValInt prxyn)

meanNormalBaseMeasure0 :: (KnownNat n, KnownNat d, Dense x) => Proxy (n/d) -> Proxy (MeanNormal (n/d)) -> Sample (MeanNormal (n/d)) -> x
meanNormalBaseMeasure0 prxyr _ x0 =
    let x = realToFrac x0
        vr = realToFrac $ ratVal prxyr
     in (exp . negate $ 0.5 * x^2 / vr) / sqrt (2*pi*vr)

multivariateNormalBaseMeasure0 :: (KnownNat n, Dense x) => Proxy n -> Proxy (MultivariateNormal n) -> Vector n Double -> x
multivariateNormalBaseMeasure0 prxyn _ _ =
    let n = natValInt prxyn
     in (2*pi)**(-fromIntegral n/2)
--- Instances ---


-- Uniform --

{-
instance Manifold Uniform where
    type Dimension Uniform = 0

instance Statistical Uniform where
    type SampleSpace Uniform = Continuum
    sampleSpace _ = Continuum

instance Generative Source Uniform where
    generate p =
        let (Uniform a b) = manifold p
         in uniformR (a,b)

instance AbsolutelyContinuous Source Uniform where
    density p x =
        let (Uniform a b) = manifold p
         in if x >= a && x <= b
               then realToFrac $ recip $ b - a
               else 0
-}

-- Bernoulli Distribution --

instance Manifold Bernoulli where
    type Dimension Bernoulli = 1

instance Statistical Bernoulli where
    type Sample Bernoulli = Bool

instance Discrete Bernoulli where
    type Cardinality Bernoulli = 2
    sampleSpace _ = doubleton True False

instance ExponentialFamily Bernoulli where
    baseMeasure _ _ = 1
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic True = Point $ singleton 1
    sufficientStatistic False = Point $ singleton 0

instance Legendre Natural Bernoulli where
    {-# INLINE potential #-}
    potential p = log $ 1 + exp (headV $ coordinates p)

instance Legendre Mean Bernoulli where
    {-# INLINE potential #-}
    potential p =
        let eta = headV $ coordinates p
         in logit eta * eta - log (1 / (1 - eta))

instance Transition Source Mean Bernoulli where
    transition = breakChart

instance Transition Mean Source Bernoulli where
    transition = breakChart

instance Transition Source Natural Bernoulli where
    transition = dualTransition . toMean

instance Transition Natural Source Bernoulli where
    transition = transition . dualTransition

instance (Transition c Source Bernoulli) => Generative c Bernoulli where
    {-# INLINE generate #-}
    generate = bernoulli . realToFrac . headV . coordinates . toSource

instance Transition Mean c Bernoulli => MaximumLikelihood c Bernoulli where
    mle = transition . sufficientStatisticT

instance AbsolutelyContinuous Source Bernoulli where
    density (Point p) True = headV p
    density (Point p) False = 1 - headV p

instance AbsolutelyContinuous Mean Bernoulli where
    density = density . toSource

instance AbsolutelyContinuous Natural Bernoulli where
    density = exponentialFamilyDensity


-- Binomial Distribution --

instance KnownNat n => Manifold (Binomial n) where
    type Dimension (Binomial n) = 1

instance KnownNat n => Statistical (Binomial n) where
    type Sample (Binomial n) = Int

instance KnownNat n => Discrete (Binomial n) where
    type Cardinality (Binomial n) = n + 1
    sampleSpace _ = generateV id

instance KnownNat n => ExponentialFamily (Binomial n) where
    baseMeasure = binomialBaseMeasure0 Proxy
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic = Point . singleton . fromIntegral

instance KnownNat n => Legendre Natural (Binomial n) where
    {-# INLINE potential #-}
    potential p =
        let n = fromIntegral $ binomialTrials p
            tht = headV $ coordinates p
         in n * log (1 + exp tht)

instance KnownNat n => Legendre Mean (Binomial n) where
    {-# INLINE potential #-}
    potential p =
        let n = fromIntegral $ binomialTrials p
            eta = headV $ coordinates p
        in eta * log (eta / (n - eta)) - n * log (n / (n - eta))

instance KnownNat n => Transition Source Natural (Binomial n) where
    transition = dualTransition . toMean

instance KnownNat n => Transition Natural Source (Binomial n) where
    transition = transition . dualTransition

instance KnownNat n => Transition Source Mean (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
         in breakChart $ n .> p

instance KnownNat n => Transition Mean Source (Binomial n) where
    transition p =
        let n = fromIntegral $ binomialTrials p
         in breakChart $ n /> p


instance (KnownNat n, Transition c Source (Binomial n)) => Generative c (Binomial n) where
    generate p0 = do
        let p = toSource p0
            n = binomialTrials p
        bls <- replicateM n . bernoulli . realToFrac . headV $ coordinates p
        return $ sum [ if bl then 1 else 0 | bl <- bls ]

instance KnownNat n => AbsolutelyContinuous Source (Binomial n) where
    density p k =
        let n = binomialTrials p
            c = headV $ coordinates p
         in realToFrac (choose n k) * c^k * (1 - c)^(n-k)

instance KnownNat n => AbsolutelyContinuous Mean (Binomial n) where
    density = density . toSource

instance KnownNat n => AbsolutelyContinuous Natural (Binomial n) where
    density = exponentialFamilyDensity

instance (KnownNat n, Transition Mean c (Binomial n)) => MaximumLikelihood c (Binomial n) where
    mle = transition . sufficientStatisticT

-- Categorical Distribution --

instance (KnownNat n, 1 <= n) => Manifold (Categorical e n) where
    type Dimension (Categorical e n) = n - 1

instance (KnownNat n, 1 <= n) => Statistical (Categorical e n) where
    type Sample (Categorical e n) = e

instance (Enum e, KnownNat n, 1 <= n) => Discrete (Categorical e n) where
    type Cardinality (Categorical e n) = n
    sampleSpace _ = generateV toEnum

instance (Enum e, KnownNat n, 1 <= n) => ExponentialFamily (Categorical e n) where
    baseMeasure _ _ = 1
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic k = Point $ generateV (\i -> if i == fromEnum k then 1 else 0)

instance (Enum e, KnownNat n, 1 <= n) => Legendre Natural (Categorical e n) where
    {-# INLINE potential #-}
    potential (Point cs) = log $ 1 + sum (exp <$> cs)

instance (Enum e, KnownNat n, 1 <= n) => Legendre Mean (Categorical e n) where
    {-# INLINE potential #-}
    potential (Point cs) =
        let scs = 1 - sum cs
         in sum (zipWithV (*) cs $ log <$> cs) + scs * log scs

instance Transition Source Mean (Categorical e n) where
    transition = breakChart

instance Transition Mean Source (Categorical e n) where
    transition = breakChart

instance (Enum e, KnownNat n, 1 <= n) => Transition Source Natural (Categorical e n) where
    transition = dualTransition . toMean

instance (Enum e, KnownNat n, 1 <= n) => Transition Natural Source (Categorical e n) where
    transition = transition . dualTransition

instance (Enum e, KnownNat n, 1 <= n, Transition c Source (Categorical e n)) => Generative c (Categorical e n) where
    generate p0 =
        let p = toSource p0
         in generateCategorical (categories p) (coordinates p)

instance (KnownNat n, 1 <= n, Enum e, Transition Mean c (Categorical e n)) => MaximumLikelihood c (Categorical e n) where
    mle = transition . sufficientStatisticT

instance (Enum e, KnownNat n, 1 <= n) => AbsolutelyContinuous Source (Categorical e n) where
    density (Point ps) e =
        let k = fromEnum e
            vi = generateV (\i -> if i == k then 1 else 0)
         in sum $ zipWithV (*) vi ps

instance (KnownNat n, 1 <= n, Enum e) => AbsolutelyContinuous Mean (Categorical e n) where
    density = density . toSource

instance (KnownNat n, 1 <= n, Enum e) => AbsolutelyContinuous Natural (Categorical e n) where
    density = exponentialFamilyDensity


{-

-- Curved Categorical Distribution --

instance Finite s => Manifold (CurvedCategorical s) where
    dimension = length . samples

instance Finite s => Statistical (CurvedCategorical s) where
    type SampleSpace (CurvedCategorical s) = s
    sampleSpace (CurvedCategorical s) = s

instance Finite s => Generative Source (CurvedCategorical s) where
    generate p = generateCategorical (samples $ manifold p) (coordinates p)

instance Finite s => AbsolutelyContinuous Source (CurvedCategorical s) where
    density p k = cs C.! idx
          where ks = samples $ manifold p
                cs = coordinates p
                idx = fromMaybe (error "attempted to calculate density of non-categorical element")
                    $ elemIndex k ks
                    -}

-- Poisson Distribution --

instance Manifold Poisson where
    type Dimension Poisson = 1

instance Statistical Poisson where
    type Sample Poisson = Int

instance ExponentialFamily Poisson where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic = Point . singleton . fromIntegral
    baseMeasure _ k = recip . realToFrac $ factorial k

instance Legendre Natural Poisson where
    {-# INLINE potential #-}
    potential = exp . headV . coordinates

instance Legendre Mean Poisson where
    {-# INLINE potential #-}
    potential (Point xs) =
        let eta = headV xs
         in eta * log eta - eta

instance Transition Source Natural Poisson where
    transition = transition . toMean

instance Transition Natural Source Poisson where
    transition = transition . dualTransition

instance Transition Source Mean Poisson where
    transition = breakChart

instance Transition Mean Source Poisson where
    transition = breakChart

instance (Transition c Source Poisson) => Generative c Poisson where
    generate = generatePoisson . headV . coordinates . toSource

instance AbsolutelyContinuous Source Poisson where
    density (Point xs) k =
        let lmda = headV xs
         in  lmda^k / realToFrac (factorial k) * exp (-lmda)

instance AbsolutelyContinuous Mean Poisson where
    density = density . toSource

instance AbsolutelyContinuous Natural Poisson where
    density = exponentialFamilyDensity

instance Transition Mean c Poisson => MaximumLikelihood c Poisson where
    mle = transition . sufficientStatisticT

-- Normal Distribution --

instance Manifold Normal where
    type Dimension Normal = 2

instance Statistical Normal where
    type Sample Normal = Double

instance ExponentialFamily Normal where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic x = fmap realToFrac . Point . doubleton x $ x**2
    baseMeasure _ _ = recip . sqrt $ 2 * pi

instance Legendre Natural Normal where
    {-# INLINE potential #-}
    potential (Point cs) =
        let (tht0,tht1) = toPair cs
         in -(tht0^2 / (4*tht1)) - 0.5 * log(-2*tht1)

instance Legendre Mean Normal where
    {-# INLINE potential #-}
    potential (Point cs) =
        let (eta0,eta1) = toPair cs
         in -0.5 * log(eta1 - eta0^2) - 1/2

instance Riemannian Natural Normal where
    metric = legendreMetric

instance Riemannian Mean Normal where
    metric = legendreMetric

instance Transition Source Mean Normal where
    transition (Point cs) =
        let (mu,vr) = toPair cs
         in Point . doubleton mu $ vr + mu^2

instance Transition Mean Source Normal where
    transition (Point cs) =
        let (eta0,eta1) = toPair cs
         in Point . doubleton eta0 $ eta1 - eta0^2

instance Transition Source Natural Normal where
    transition (Point cs) =
        let (mu,vr) = toPair cs
         in Point $ doubleton (mu / vr) (negate . recip $ 2 * vr)

instance Transition Natural Source Normal where
    transition (Point cs) =
        let (tht0,tht1) = toPair cs
         in Point $ doubleton (-0.5 * tht0 / tht1) (negate . recip $ 2 * tht1)

instance (Transition c Source Normal) => Generative c Normal where
    generate p =
        let (Point cs) = toSource p
            (mu,vr) = toPair cs
         in normal (realToFrac mu) (realToFrac $ sqrt vr)

instance AbsolutelyContinuous Source Normal where
    density (Point cs) x =
        let (mu,vr) = toPair cs
         in recip (sqrt $ vr*2*pi) * exp (negate $ (realToFrac x - mu) ** 2 / (2*vr))

instance AbsolutelyContinuous Mean Normal where
    density = density . toSource

instance AbsolutelyContinuous Natural Normal where
    density = exponentialFamilyDensity

instance Transition Mean c Normal => MaximumLikelihood c Normal where
    mle = transition . sufficientStatisticT


{-
instance Riemannian Source Normal where
    metric (Point xs) =
         in  [recip vr,0,0,recip $ 2*vr^2]
         -}

-- MeanNormal Distribution --

instance Manifold (MeanNormal v) where
    type Dimension (MeanNormal v) = 1

instance Statistical (MeanNormal v) where
    type Sample (MeanNormal v) = Double

instance (KnownNat n, KnownNat d) => ExponentialFamily (MeanNormal (n / d)) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic x = fmap realToFrac . Point $ singleton x
    baseMeasure = meanNormalBaseMeasure0 Proxy

instance (KnownNat n, KnownNat d) => Legendre Natural (MeanNormal (n/d)) where
    {-# INLINE potential #-}
    potential p =
        let vr = realToFrac $ meanNormalVariance p
            mu = headV $ coordinates p
         in 0.5 * vr * mu^2

instance (KnownNat n, KnownNat d) => Legendre Mean (MeanNormal (n/d)) where
    {-# INLINE potential #-}
    potential p =
        let vr = realToFrac $ meanNormalVariance p
            mu = headV $ coordinates p
         in 0.5 / vr * mu^2

instance Transition Source Mean (MeanNormal v) where
    transition = breakChart

instance Transition Mean Source (MeanNormal v) where
    transition = breakChart

instance (KnownNat n, KnownNat d) => Transition Source Natural (MeanNormal (n/d)) where
    transition = dualTransition . toMean

instance (KnownNat n, KnownNat d) => Transition Natural Source (MeanNormal (n/d)) where
    transition = toSource . dualTransition

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Source (MeanNormal (n/d)) where
    density p =
        let vr = realToFrac $ meanNormalVariance p
            mu = headV $ coordinates p
            nrm :: x -> x -> Point Source Normal x
            nrm x y = Point $ doubleton x y
         in density $ nrm mu vr

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Mean (MeanNormal (n/d)) where
    density = density . toSource

instance (KnownNat n, KnownNat d) => AbsolutelyContinuous Natural (MeanNormal (n/d)) where
    density = exponentialFamilyDensity

instance (KnownNat n, KnownNat d, Transition Mean c (MeanNormal (n/d))) => MaximumLikelihood c (MeanNormal (n/d)) where
    mle = transition . sufficientStatisticT

instance (KnownNat n, KnownNat d, Transition c Source (MeanNormal (n/d))) => Generative c (MeanNormal (n/d)) where
    generate p =
        let (Point cs) = toSource p
            mu = headV cs
            vr = realToFrac $ meanNormalVariance p
         in normal (realToFrac mu) (realToFrac $ sqrt vr)

-- Multivariate Normal --

instance KnownNat n => Manifold (MultivariateNormal n) where
    type Dimension (MultivariateNormal n) = n + n * n

instance KnownNat n => Statistical (MultivariateNormal n) where
    type Sample (MultivariateNormal n) = Vector n Double

instance KnownNat n => ExponentialFamily (MultivariateNormal n) where
    {-# INLINE sufficientStatistic #-}
    sufficientStatistic xs =
        let Matrix cvrs = matrixMatrixMultiply (columnVector xs) (rowVector xs)
         in fmap realToFrac . Point $ joinV xs cvrs
    baseMeasure = multivariateNormalBaseMeasure0 Proxy

{-
instance Legendre Natural MultivariateNormal where
    potential p =
        let (tmu,tsgma) = splitMultivariateNormal p
            invtsgma = matrixInverse tsgma
         in -0.25 * dotProduct tmu (matrixVectorMultiply invtsgma tmu) - 0.5 * log(M.det $ M.scale (-2) tsgma)

instance Legendre Mean MultivariateNormal where
    potential p =
        let (mmu,msgma) = splitMultivariateNormal p
         in -0.5 * (1 + M.dot mmu (M.pinv msgma M.#> mmu)) - 0.5 * log (M.det msgma)

instance Transition Source Natural MultivariateNormal where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
            invsgma = M.pinv sgma
         in fromCoordinates (manifold p) $ (invsgma M.#> mu) C.++ M.flatten (M.scale (-0.5) invsgma)

instance Transition Natural Source MultivariateNormal where
    transition p =
        let (emu,esgma) = splitMultivariateNormal p
            invesgma = M.scale (-0.5) $ M.pinv esgma
         in fromCoordinates (manifold p) $ (invesgma M.#> emu) C.++ M.flatten invesgma

instance Transition Source Mean MultivariateNormal where
    transition p =
        let (mu,sgma) = splitMultivariateNormal p
         in fromCoordinates (manifold p) $ mu C.++ M.flatten (sgma + M.outer mu mu)

instance Transition Mean Source MultivariateNormal where
    transition p =
        let (mmu,msgma) = splitMultivariateNormal p
         in fromCoordinates (manifold p) $ mmu C.++ M.flatten (msgma -M.outer mmu mmu)

instance Generative Source MultivariateNormal where
    generate p =
        let n = sampleSpaceDimension $ manifold p
            (mus,sds) = C.splitAt n $ coordinates p
         in generateMultivariateNormal mus $ M.reshape n sds

instance AbsolutelyContinuous Source MultivariateNormal where
    density p xs =
        let n = sampleSpaceDimension $ manifold p
            (mus,sgma) = splitMultivariateNormal p
         in recip ((2*pi)**(fromIntegral n / 2) * sqrt (M.det sgma))
            * exp (-0.5 * ((M.tr (M.pinv sgma) M.#> C.zipWith (-) xs mus) `M.dot` C.zipWith (-) xs mus))

instance MaximumLikelihood Source MultivariateNormal where
    mle _ xss =
        let n = fromIntegral $ length xss
            mus = recip (fromIntegral n) * sum xss
            sgma = recip (fromIntegral $ n - 1)
                * sum (map (\xs -> let xs' = xs - mus in M.outer xs' xs') xss)
        in  joinMultivariateNormal mus sgma
        -}

-- VonMises --

instance Manifold VonMises where
    type Dimension VonMises = 2

instance Statistical VonMises where
    type Sample VonMises = Double

instance Generative Source VonMises where
    generate p@(Point cs) = do
        let (mu,kap) = toPair $ realToFrac <$> cs
            tau = 1 + sqrt (1 + 4 * kap^2)
            rho = (tau - sqrt (2*tau))/(2*kap)
            r = (1 + rho^2) / (2 * rho)
        [u1,u2,u3] <- replicateM 3 uniform
        let z = cos (pi * u1)
            f = (1 + r * z)/(r + z)
            c = kap * (r - f)
        if log (c / u2) + 1 - c < 0
           then generate p
           else return . toPi $ signum (u3 - 0.5) * acos f + mu

instance ExponentialFamily VonMises where
    sufficientStatistic tht = fmap realToFrac . Point $ doubleton (cos tht) (sin tht)
    baseMeasure _ _ = recip $ 2 * pi

instance Transition Source Natural VonMises where
    transition (Point cs) =
        let (mu,kap) = toPair cs
         in Point $ doubleton (kap * cos mu) (kap * sin mu)

instance Transition Natural Source VonMises where
    transition (Point cs) =
        let (tht0,tht1) = toPair cs
         in Point $ doubleton (atan2 tht1 tht0) (sqrt $ tht0^2 + tht1^2)

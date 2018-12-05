{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeOperators
    #-}

module Goal.NeuralData
    ( -- * Types
      prjnm
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    -- * General
    , stimulusResponseMap
    , meanSDInliers
    -- * Subsampling
    , generateIndices
    , subSampleResponses
    -- * Empirical Analysis
    , empiricalTuningCurves
    , subsampleEmpiricalTuningCurves
    , empiricalPosterior0
    -- * Independent Poisson Analysis
    , getFittedIPLikelihood
    , strengthenIPLikelihood
    , fitIPLikelihood
    , fitLinearDecoder
    , subIPLikelihood
    , subsampleIPLikelihood
    -- ** Partition Functions
    , conditionalIPLogPartitionFunction
    , linearConditionalIPLogPartitionFunction
    , affineConditionalIPLogPartitionFunction
    -- ** Inference
    , numericalVonMisesPPCPosterior
    , approximateVonMisesPPCPosterior
    , correctedVonMisesPPCPosterior
    -- ** Statistics
    , fisherInformation
    , averageLogFisherInformation
    -- * Mixture GLM Analysis
    , fitMixtureLikelihood
    , getFittedMixtureLikelihood
    , strengthenMixtureLikelihood
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Numeric.Sum as NS

-- Qualified --

import qualified Data.Map as M


--- Types ---


--- CSV ---

prjnm :: String
prjnm = "neural-data"

--- Inference ---


-- Under the assumption of a flat prior
unnormalizedEmpiricalPPCLogPosterior0
    :: KnownNat k
    => Bool
    -> Mean # Neurons k
    -> Response k
    -> Double
unnormalizedEmpiricalPPCLogPosterior0 True mz z =
     dualTransition mz <.> sufficientStatistic z - NS.sum NS.kbn (listCoordinates mz)
unnormalizedEmpiricalPPCLogPosterior0 False mz z =
     dualTransition mz <.> sufficientStatistic z

-- Under the assumption of a flat prior
empiricalPosterior0
    :: (Ord s, KnownNat k)
    => Bool -> M.Map s (Mean # Neurons k) -> Response k -> M.Map s Double
empiricalPosterior0 nrmb xzmp z =
    let uldns = flip (unnormalizedEmpiricalPPCLogPosterior0 nrmb) z <$> xzmp
        avg = NS.sum NS.kbn uldns / fromIntegral (length uldns)
        udns = exp . subtract avg <$> uldns
        nrm = traceGiven $ NS.sum NS.kbn udns
     in (/nrm) <$> udns

nstms :: Int
nstms = 100

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) (nstms+1)

errbnd :: Double
errbnd = 1e-12

linearConditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
linearConditionalIPLogPartitionFunction lkl z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
        logupst x = sz <.> (nzx >.>* x) - log (2*pi)
        mx = maximum $ logupst <$> xsmps
        upst0 x = exp $ logupst x - mx
     in (nz <.> sz +) . (+ mx) . log1p . subtract 1 . fst $ integrate errbnd upst0 0 (2*pi)

affineConditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Response k
    -> Double
affineConditionalIPLogPartitionFunction lkl rprms z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
     in (nz <.> sz +) . potential $ sz <.< nzx <-> rprms

conditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
conditionalIPLogPartitionFunction lkl z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
        logupst x = sz <.> (nzx >.>* x) - potential (lkl >.>* x) - log (2*pi)
        mx = maximum $ logupst <$> xsmps
        upst0 x = exp $ logupst x - mx
     in (nz <.> sz +) . (+ mx) . log1p . subtract 1 . fst $ integrate errbnd upst0 0 (2*pi)

-- Under the assumption of a flat prior
numericalVonMisesPPCPosterior
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
numericalVonMisesPPCPosterior ppc z =
    let logupst x = sufficientStatistic z <.> (ppc >.>* x) - potential (ppc >.>* x)
        nrm = conditionalIPLogPartitionFunction ppc z
     in \x -> exp $ logupst x - nrm

-- Under the assumption of a flat prior
approximateVonMisesPPCPosterior
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
approximateVonMisesPPCPosterior ppc z =
    let logupst x = sufficientStatistic z <.> (ppc >.>* x)
        nrm = linearConditionalIPLogPartitionFunction ppc z
     in \x -> exp $ logupst x - nrm

-- Under the assumption of a flat prior
correctedVonMisesPPCPosterior
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Response k
    -> Double
    -> Double
correctedVonMisesPPCPosterior ppc rprms z =
    let logupst x = sufficientStatistic z <.> (ppc >.>* x) - rprms <.> sufficientStatistic x
        nrm = linearConditionalIPLogPartitionFunction ppc z
     in \x -> exp $ logupst x - nrm

getNeuralData :: Read s => String -> Dataset -> IO (NatNumber,[([Int], s)])
getNeuralData expnm dst = read <$> goalReadDataset prjnm expnm dst

getFittedIPLikelihood
    :: String
    -> Dataset
    -> IO (Int,[Double])
getFittedIPLikelihood expnm dst = do
    (k,xs) <- read <$> goalReadDataset prjnm expnm dst
    return (k,xs)

getFittedMixtureLikelihood
    :: String
    -> Dataset
    -> IO (Int,Int,[Double])
getFittedMixtureLikelihood expnm dst = do
    (k,n,xs) <- read <$> goalReadDataset prjnm expnm dst
    return (k,n,xs)

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . B.fromList <$> ks) ss

strengthenIPLikelihood
    :: KnownNat k
    => [Double]
    -> Mean #> Natural # Neurons k <* VonMises
strengthenIPLikelihood xs = Point . fromJust $ S.fromList xs

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n, 1 <= n)
    => [Double]
    -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---


meanSDInliers :: [Double] -> (Double,Double)
meanSDInliers xs =
    let (mu,vr) = estimateMeanVariance xs
        xs' = filter (\x -> square (x-mu) < 4*vr) xs
        (mu',vr') = estimateMeanVariance xs'
     in (mu',sqrt vr')


stimulusResponseMap :: Ord s => [(Response k, s)] -> M.Map s [Response k]
stimulusResponseMap zxs = M.fromListWith (++) [(x,[z]) | (z,x) <- zxs]

empiricalTuningCurves :: (Ord s, KnownNat k) => M.Map s [Response k] -> M.Map s (Mean # Neurons k)
empiricalTuningCurves zxmp = sufficientStatisticT <$> zxmp

ppcStimulusDerivatives
    :: KnownNat k => Mean #> Natural # Neurons k <* VonMises
    -> SamplePoint VonMises
    -> S.Vector k Double
ppcStimulusDerivatives ppc x =
    let fxs = coordinates . dualTransition $ ppc >.> mx
        tcs = toRows . snd $ splitAffine ppc
     in S.zipWith zipper fxs tcs
    where mx = sufficientStatistic x
          (cx,sx) = S.toPair $ coordinates mx
          zipper fx (Point cs) =
              let (tht1,tht2) = S.toPair cs
               in fx*(cx * tht2 - sx * tht1)

fisherInformation
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Double
    -> Double
fisherInformation ppc x =
    let fxs2' = S.map square $ ppcStimulusDerivatives ppc x
        fxs = coordinates . dualTransition $ ppc >.>* x
     in S.sum $ S.zipWith (/) fxs2' fxs

averageLogFisherInformation
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Double
averageLogFisherInformation ppc =
    average $ log . (/(2*pi*exp 1)) . fisherInformation ppc <$> xsmps


--- Subsampling ---

fitIPLikelihood
    :: forall r k . KnownNat k
    => [(Response k,Double)]
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
fitIPLikelihood xzs = do
    let eps = -0.1
        nepchs = 500
    kps <- S.replicateM $ uniformR (0.2,0.6)
    let sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps $ S.range 0 (2*pi)
    gns' <- Point <$> S.replicateM (uniformR (0,2))
    let gns0 = transition . sufficientStatisticT $ fst <$> xzs
        gns = gns0 <+> gns'
        ppc0 = vonMisesPopulationEncoder True (Right gns) sps
        (zs,xs) = unzip xzs
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
    return (vanillaGradientSequence backprop eps defaultAdamPursuit ppc0 !! nepchs)

-- NB: Actually affine, not linear
fitLinearDecoder
    :: forall k . KnownNat k
    => [(Response k,Double)]
    -> Mean #> Natural # VonMises <* Neurons k
fitLinearDecoder xzs =
    let eps = -0.1
        nepchs = 500
        sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        nxz = transpose . fromRows $ S.map toNatural sps
        nx = Point $ S.fromTuple (0,0.5)
        aff0 = joinAffine nx nxz
        (zs,xs) = unzip xzs
        backprop aff = joinTangentPair aff $ stochasticConditionalCrossEntropyDifferential zs xs aff
     in (vanillaGradientSequence backprop eps defaultAdamPursuit aff0 !! nepchs)

fitMixtureLikelihood
    :: forall k n r . (KnownNat k, KnownNat n, 1 <= n)
    => [(Response k,Double)]
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises) -- ^ Function
fitMixtureLikelihood xzs = do
    let eps = -0.05
        nepchs = 500
    kps <- S.replicateM $ uniformR (0.2,0.6)
    let sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps $ S.range 0 (2*pi)
        wghts0 :: Natural # Categorical Int n
        wghts0 = zero
        wghts = toSource wghts0
        gnss0 = S.replicate . transition . sufficientStatisticT $ fst <$> xzs
    gnss' <- S.replicateM . fmap Point . S.replicateM $ uniformR (0,2)
    let gnss = S.zipWith (<+>) gnss0 gnss'
        ppc0 = vonMisesMixturePopulationEncoder True wghts gnss sps
        (zs,xs) = unzip xzs
        backprop mglm = joinTangentPair mglm $ mixtureStochasticConditionalCrossEntropyDifferential xs zs mglm
    return (vanillaGradientSequence backprop eps defaultAdamPursuit ppc0 !! nepchs)


subIPLikelihood
    :: forall k k' . (KnownNat k, KnownNat k', k' <= k)
    => Mean #> Natural # Neurons k <* VonMises
    ->  Mean #> Natural # Neurons k' <* VonMises
subIPLikelihood ppc =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . S.take . S.toRows $ toMatrix tns
        (bs',_ :: S.Vector (k - k') Double) = S.splitAt $ coordinates bs
     in joinAffine (Point bs') tns'

subsampleIPLikelihood
    :: (KnownNat k, KnownNat k')
    => Mean #> Natural # Neurons k <* VonMises
    -> S.Vector k' Int
    ->  Mean #> Natural # Neurons k' <* VonMises
subsampleIPLikelihood ppc idxs =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . flip S.backpermute idxs . S.toRows $ toMatrix tns
        bs' = Point . flip S.backpermute idxs $ coordinates bs
     in joinAffine bs' tns'

subsampleEmpiricalTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subsampleEmpiricalTuningCurves nzxmp idxs =
     Point . flip S.backpermute idxs . coordinates <$> nzxmp

subSampleResponses
    :: (Ord s, KnownNat k)
    => M.Map s [Response k]
    -> B.Vector k' Int
    -> M.Map s [Response k']
subSampleResponses zxmp idxs =
     map (`B.backpermute` idxs) <$> zxmp

generateIndices
    :: forall k m r . (KnownNat k, KnownNat m)
    => Proxy (k + m)
    -> Random r (B.Vector k Int)
generateIndices _ = do
    let idxs :: B.Vector (k + m) Int
        idxs = B.generate finiteInt
    subsampleVector idxs

--indexChain
--    :: (KnownNat k, KnownNat k', k' <= k)
--    => Int
--    -> Proxy k
--    -> Random r (Chain ([B.Vector k' Int]))
--indexChain nsmps prxk =
--    accumulateRandomChain (replicateM nsmps $ generateIndices prxk)

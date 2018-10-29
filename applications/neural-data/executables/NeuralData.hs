{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE GADTs,ScopedTypeVariables,DataKinds,TypeOperators #-}

module NeuralData
    ( -- * Types
    NeuralModel
    , Neurons
    , Response
    , prjnm
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    , getPopulationSize
    -- * General
    , stimulusResponseMap
    , meanSDInliers
    -- * Subsampling
    , generateIndices
    , subSampleResponses
    -- * Empirical Analysis
    , empiricalTuningCurves
    , subSampleEmpiricalTuningCurves
    , empiricalPPCPosterior0
    -- * Von Mises Analysis
    , fitPPC
    , subPPC
    , subSamplePPC
    , approximateConditionalLogPartitionFunction
    , correctedConditionalLogPartitionFunction
    , conditionalLogPartitionFunction
    , numericalVonMisesPPCPosterior
    , approximateVonMisesPPCPosterior
    , correctedVonMisesPPCPosterior
    , fisherInformation
    , averageLogFisherInformation
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


type NeuralModel s k = Harmonium Tensor (Replicated k Poisson) s

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
empiricalPPCPosterior0
    :: (Ord s, KnownNat k)
    => Bool -> M.Map s (Mean # Neurons k) -> Response k -> M.Map s Double
empiricalPPCPosterior0 nrmb xzmp z =
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

approximateConditionalLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
approximateConditionalLogPartitionFunction lkl z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
        logupst x = sz <.> (nzx >.>* x)
        mx = maximum $ logupst <$> xsmps
        upst0 x = exp $ logupst x - mx
     in (nz <.> sz +) . (+ mx) . log1p . subtract 1 . fst $ integrate errbnd upst0 0 (2*pi)

correctedConditionalLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Response k
    -> Double
correctedConditionalLogPartitionFunction lkl rprms z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
        logupst x = sz <.> (nzx >.>* x) - sufficientStatistic x <.> rprms
        mx = maximum $ logupst <$> xsmps
        upst0 x = exp $ logupst x - mx
     in (nz <.> sz +) . (+ mx) . log1p . subtract 1 . fst $ integrate errbnd upst0 0 (2*pi)

conditionalLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
conditionalLogPartitionFunction lkl z =
    let (nz,nzx) = splitAffine lkl
        sz = sufficientStatistic z
        logupst x = sz <.> (nzx >.>* x) - potential (lkl >.>* x)
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
        nrm = conditionalLogPartitionFunction ppc z
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
        nrm = approximateConditionalLogPartitionFunction ppc z
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
        nrm = approximateConditionalLogPartitionFunction ppc z
     in \x -> exp $ logupst x - nrm

getNeuralData :: Read s => String -> Dataset -> IO [([Int], s)]
getNeuralData expnm dst =
    read <$> goalReadDataset prjnm expnm dst

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . B.fromList <$> ks) ss

getPopulationSize :: [([Int], s)] -> Int
getPopulationSize = length . fst . head


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

fitPPC
    :: forall k . KnownNat k
    => [(Response k,Double)]
    -> Mean #> Natural # Neurons k <* VonMises
fitPPC xzs =
    let eps = -0.1
        nepchs = 500
        sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        ppc0 = vonMisesPopulationEncoder True (Left 1) sps
        (zs,xs) = unzip xzs
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
     in (vanillaGradientSequence backprop eps defaultAdamPursuit ppc0 !! nepchs)


subPPC
    :: forall k k' . (KnownNat k, KnownNat k', k' <= k)
    => Mean #> Natural # Neurons k <* VonMises
    ->  Mean #> Natural # Neurons k' <* VonMises
subPPC ppc =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . S.take . S.toRows $ toMatrix tns
        (bs',_ :: S.Vector (k - k') Double) = S.splitAt $ coordinates bs
     in joinAffine (Point bs') tns'

subSamplePPC
    :: (KnownNat k, KnownNat k')
    => Mean #> Natural # Neurons k <* VonMises
    -> S.Vector k' Int
    ->  Mean #> Natural # Neurons k' <* VonMises
subSamplePPC ppc idxs =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . flip S.backpermute idxs . S.toRows $ toMatrix tns
        bs' = Point . flip S.backpermute idxs $ coordinates bs
     in joinAffine bs' tns'

subSampleEmpiricalTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subSampleEmpiricalTuningCurves nzxmp idxs =
     Point . flip S.backpermute idxs . coordinates <$> nzxmp

subSampleResponses
    :: (Ord s, KnownNat k)
    => M.Map s [Response k]
    -> B.Vector k' Int
    -> M.Map s [Response k']
subSampleResponses zxmp idxs =
     map (`B.backpermute` idxs) <$> zxmp

generateIndices
    :: forall k k' r . (KnownNat k, KnownNat k', k' <= k)
    => Proxy k
    -> Random r (B.Vector k' Int)
generateIndices _ = do
    let idxs :: B.Vector k Int
        idxs = B.generate finiteInt
    subsampleVector idxs

--indexChain
--    :: (KnownNat k, KnownNat k', k' <= k)
--    => Int
--    -> Proxy k
--    -> Random r (Chain ([B.Vector k' Int]))
--indexChain nsmps prxk =
--    accumulateRandomChain (replicateM nsmps $ generateIndices prxk)

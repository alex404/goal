{-# LANGUAGE GADTs,ScopedTypeVariables,DataKinds,TypeOperators,DeriveGeneric #-}

module NeuralData
    ( -- * Types
    NeuralModel
    , Neurons
    , Response
    -- ** CSV
    , CoefficientsOfVariation (CoefficientsOfVariation)
    , VonMisesInformations (VonMisesInformations)
    , DivergenceStatistics (DivergenceStatistics,posteriorDivergence)
    , Posteriors (Posteriors)
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    , getPopulationSize
    -- * General
    , stimulusResponseMap
    -- * Subsampling
    , generateIndices
    , subSampleResponses
    -- * Empirical Analysis
    , empiricalTuningCurves
    , subSampleEmpiricalTuningCurves
    , empiricalPPCPosterior0
    -- * Von Mises Analysis
    , fitPPC
    , subSamplePPC
    , numericalVonMisesPPCPosterior
    , approximateVonMisesPPCPosterior
    , examplePosteriors
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
type Neurons k = Replicated k Poisson
type Response k = SamplePoint (Neurons k)

--- CSV ---


data CoefficientsOfVariation = CoefficientsOfVariation
    { meanResponseCV :: Double
    , sdResponseCV :: Double
    , meanTuningCurveCV :: Double
    , sdTuningCurveCV :: Double }
    deriving (Generic, Show)

instance FromNamedRecord CoefficientsOfVariation
instance ToNamedRecord CoefficientsOfVariation
instance DefaultOrdered CoefficientsOfVariation
instance NFData CoefficientsOfVariation

data VonMisesInformations = VonMisesInformations
    { averagedFisherInformation :: Double
    , minimalCVFisherInformation :: Double
    , maximalCVFisherInformation :: Double }
    deriving (Generic, Show)

instance FromNamedRecord VonMisesInformations
instance ToNamedRecord VonMisesInformations
instance DefaultOrdered VonMisesInformations
instance NFData VonMisesInformations

newtype DivergenceStatistics = DivergenceStatistics
    { posteriorDivergence :: Double }
    deriving (Generic, Show)

instance FromNamedRecord DivergenceStatistics
instance ToNamedRecord DivergenceStatistics
instance DefaultOrdered DivergenceStatistics
instance NFData DivergenceStatistics

data Posteriors = Posteriors
    { stimulusSamples :: Double
    , numericalPosterior :: Double
    , approximatePosterior :: Double }
    deriving (Generic, Show)

instance FromNamedRecord Posteriors
instance ToNamedRecord Posteriors
instance DefaultOrdered Posteriors
instance NFData Posteriors


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

-- Under the assumption of a flat prior
numericalVonMisesPPCPosterior
    :: KnownNat k
    => Double
    -> Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
numericalVonMisesPPCPosterior err ppc z =
    let nx0 = z *<.< snd (splitAffine ppc)
        upst x0 = exp $ nx0 <.> sufficientStatistic x0 - head (sumOfTuningCurves ppc [x0])
        nrm = fst $ integrate err upst 0 (2*pi)
     in (/nrm) . upst

-- Under the assumption of a flat prior
approximateVonMisesPPCPosterior
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
approximateVonMisesPPCPosterior ppc z =
    density (z *<.< snd (splitAffine ppc))

-- Under the assumption of a flat prior
examplePosteriors
    :: KnownNat k
    => Double
    -> [Double]
    -> Int
    -> Mean #> Natural # Neurons k <* VonMises
    -> Random r [[Posteriors]]
examplePosteriors err xs nplts ppc = do
    let plts = range 0 (2*pi) nplts
    zs <- mapM samplePoint $ ppc >$>* xs
    return $ do
        z <- zs
        let npstrs = numericalVonMisesPPCPosterior err ppc z <$> plts
            apstrs = approximateVonMisesPPCPosterior ppc z <$> plts
        return $ zipWith3 Posteriors plts npstrs apstrs

getNeuralData :: Read s => Collection -> Dataset -> IO [([Int], s)]
getNeuralData (Collection clcstr) (Dataset dststr) =
    read <$> readFile ("projects/" ++ clcstr ++ "/data/" ++ dststr ++ ".dat")

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . B.fromList <$> ks) ss

getPopulationSize :: [([Int], s)] -> Int
getPopulationSize = length . fst . head


--- Analysis ---


stimulusResponseMap :: Ord s => [(Response k, s)] -> M.Map s [Response k]
stimulusResponseMap zxs = M.fromListWith (++) [(x,[z]) | (z,x) <- zxs]

empiricalTuningCurves :: (Ord s, KnownNat k) => M.Map s [Response k] -> M.Map s (Mean # Neurons k)
empiricalTuningCurves zxmp = sufficientStatisticT <$> zxmp


--- Subsampling ---

fitPPC
    :: forall k . KnownNat k
    => Int
    -> Double
    -> [(Response k,Double)]
    -> Mean #> Natural # Neurons k <* VonMises
fitPPC nepchs eps xzs =
    let sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        ppc0 = vonMisesPopulationEncoder True (Left 1) sps
        (zs,xs) = unzip xzs
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
     in (vanillaGradientSequence backprop eps defaultAdamPursuit ppc0 !! nepchs)


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

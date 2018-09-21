{-# LANGUAGE GADTs,ScopedTypeVariables,DataKinds,TypeOperators,DeriveGeneric #-}

module NeuralData
    ( -- * Types
    NeuralModel
    , Neurons
    , Response
    -- ** CSV
    , CoefficientsOfVariation (CoefficientsOfVariation)
    , VonMisesInformations (VonMisesInformations)
    , averageCoVs
    , averageVMIs
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    , getPopulationSize
    -- * Analysis
    , stimulusResponseMap
    , empiricalTuningCurves
    , responseStatistics
    -- * Inference
    , empiricalPPCPosterior0
    , numericalVonMisesPPCPosterior
    , approximateVonMisesPPCPosterior
    -- * Subsampling
    , generateIndices
    , subSampleTuningCurves
    , subSampleResponses
    , subSamplePPC
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Numeric.Sum as NS

-- Qualified --

import qualified Data.Map as M
import qualified Data.List as L


--- Types ---


type NeuralModel s k = Harmonium Tensor (Replicated k Poisson) s
type Neurons k = Replicated k Poisson
type Response k = SamplePoint (Neurons k)


--- CSV ---


data CoefficientsOfVariation = CoefficientsOfVariation
    { responseCV :: Double
    , minimumCV :: Double
    , averageCV :: Double
    , maximumCV :: Double }
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

averageCoVs :: [[CoefficientsOfVariation]] -> [CoefficientsOfVariation]
averageCoVs cvsss = mapfun . foldr1 foldfun <$> L.transpose (map (take k) cvsss)
    where
        n = L.genericLength cvsss
        k = length $ head cvsss
        foldfun (CoefficientsOfVariation a1 b1 c1 d1) (CoefficientsOfVariation a2 b2 c2 d2) =
            CoefficientsOfVariation (a1+a2) (b1+b2) (c1+c2) (d1+d2)
        mapfun (CoefficientsOfVariation a1 b1 c1 d1) =
            CoefficientsOfVariation (a1/n) (b1/n) (c1/n) (d1/n)

averageVMIs :: [[VonMisesInformations]] -> [VonMisesInformations]
averageVMIs vmisss = mapfun . foldr1 foldfun <$> L.transpose (map (take k) vmisss)
    where
        n = L.genericLength vmisss
        k = length $ head vmisss
        foldfun (VonMisesInformations a1 b1 c1) (VonMisesInformations a2 b2 c2) =
            VonMisesInformations (a1+a2) (b1+b2) (c1+c2)
        mapfun (VonMisesInformations a1 b1 c1) = VonMisesInformations (a1/n) (b1/n) (c1/n)


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
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
numericalVonMisesPPCPosterior ppc z =
    let nx0 = z *<.< snd (splitAffine ppc)
        upst x0 = nx0 <.> sufficientStatistic x0 - head (sumOfTuningCurves ppc [x0])
        nrm = integrate 1e-10 upst 0 (2*pi)
     in (/nrm) . upst

-- Under the assumption of a flat prior
approximateVonMisesPPCPosterior
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
    -> Double
approximateVonMisesPPCPosterior ppc z x =
    (z *<.< snd (splitAffine ppc)) <.> sufficientStatistic x


getNeuralData :: Read s => Collection -> Dataset -> IO [([Int], s)]
getNeuralData (Collection clcstr) (Dataset dststr) =
    read <$> readFile (clcstr ++ "/data/" ++ dststr ++ ".dat")

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

tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)

responseStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => [(Response k,s)]
    -> Int
    -> Proxy k'
    -> Random r ( CoefficientsOfVariation, ([S.Vector k' Int], S.Vector k' Int, S.Vector k' Int))
responseStatistics zxss n _ = do
    let zxmp = stimulusResponseMap zxss
        nzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let sidxss = G.convert <$> idxss
        subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleTuningCurves nzxmp <$> sidxss
        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcss
        sidxcvs = zip scvs sidxss
        (mxcv,mxsidxs) = maximumBy (comparing fst) sidxcvs
        (mncv,mnsidxs) = minimumBy (comparing fst) sidxcvs
        cvs = CoefficientsOfVariation (average rcvs) mncv (average scvs) mxcv
    return (cvs, (sidxss,mnsidxs,mxsidxs))


--- Subsampling ---


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

subSampleTuningCurves
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> S.Vector k' Int
    -> M.Map s (Mean # Neurons k')
subSampleTuningCurves nzxmp idxs =
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

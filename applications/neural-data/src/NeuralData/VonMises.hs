{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    DeriveGeneric,
    TypeOperators
    #-}

module NeuralData.VonMises
    ( -- * Independent Poisson Analysis
      getFittedIPLikelihood
    , strengthenIPLikelihood
    , fitIPLikelihood
    , fitLinearDecoder
    , subIPLikelihood
    , subsampleIPLikelihood
    , linearDecoderDivergence
    , analyzeTuningCurves
    -- ** Partition Functions
    , conditionalIPLogPartitionFunction
    , affineConditionalIPLogPartitionFunction
    -- ** Statistics
    , fisherInformation
    , averageLogFisherInformation
    , populationParameterHistogram
    , PopulationParameterCounts
    ) where


--- Imports ---


-- Goal --

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import qualified Data.List as L


--- Types ---


--- CSV ---


--- Inference ---


-- Test this to see if I can just replace the application with *<.< on the
-- affine transformation.
affineConditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Response k
    -> Double
affineConditionalIPLogPartitionFunction lkl rprms z =
     potential $ z *<.< snd (splitAffine lkl) <-> rprms

conditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
conditionalIPLogPartitionFunction lkl z =
    let sz = sufficientStatistic z
        logupst x = sz <.> (snd (splitAffine lkl) >.>* x) - potential (lkl >.>* x) - log (2*pi)
     in logIntegralExp 1e-6 logupst 0 (2*pi) (tail $ range 0 (2*pi) 100)

-- Under the assumption of a flat prior
linearDecoderDivergence
    :: KnownNat k
    => Mean #> Natural # VonMises <* Neurons k
    -> (Double -> Double) -- ^ True Density
    -> Response k
    -> Double
linearDecoderDivergence dcd trudns z =
    let dcddns x = density (dcd >.>* z) x
        dv0 x = trudns x * log (trudns x / dcddns x)
     in fst $ integrate 1e-3 dv0 mnx mxx

getFittedIPLikelihood
    :: String
    -> String
    -> IO (NatNumber,[Double])
getFittedIPLikelihood expnm dst =
    read <$> goalReadDataset (Experiment prjnm expnm) dst

strengthenIPLikelihood
    :: KnownNat k
    => [Double]
    -> Mean #> Natural # Neurons k <* VonMises
strengthenIPLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---

analyzeTuningCurves
    :: forall k . KnownNat k
    => Sample VonMises
    -> Mean #> Natural # Neurons k <* VonMises
    -> [[Double]]
analyzeTuningCurves xsmps lkl =
    let nzs = lkl >$>* xsmps
        tcss = listCoordinates . dualTransition <$> nzs
        stcs = potential <$> nzs
        (rho0,rprms) = regressRectificationParameters lkl xsmps
        rcrv = rectificationCurve rho0 rprms xsmps
        mxtcs = maximum <$> tcss
     in zipWith (++) (L.transpose (xsmps:stcs:rcrv:[mxtcs])) tcss

ppcStimulusDerivatives
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
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
    average $ log . (/(2*pi*exp 1)) . fisherInformation ppc <$> tail (range 0 (2*pi) (101))

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
    :: forall s k . KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Sample VonMises
    -> Random s (Mean #> Natural # VonMises <* Neurons k)
fitLinearDecoder lkl xs = do
    zs <- mapM samplePoint (lkl >$>* xs)
    let eps = -0.1
        nepchs = 500
        sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        nxz = transpose . fromRows $ S.map toNatural sps
        nx = Point $ S.fromTuple (0,0.5)
        aff0 = joinAffine nx nxz
        backprop aff = joinTangentPair aff $ stochasticConditionalCrossEntropyDifferential zs xs aff
    return (vanillaGradientSequence backprop eps defaultAdamPursuit aff0 !! nepchs)

subIPLikelihood
    :: forall k m . (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k + m) <* VonMises
    ->  Mean #> Natural # Neurons k <* VonMises
subIPLikelihood ppc =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . S.take . S.toRows $ toMatrix tns
        bs' = S.take $ coordinates bs
     in joinAffine (Point bs') tns'

subsampleIPLikelihood
    :: (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k+m) <* VonMises
    -> S.Vector k Int
    -> Mean #> Natural # Neurons k <* VonMises
subsampleIPLikelihood ppc idxs =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . flip S.backpermute idxs . S.toRows $ toMatrix tns
        bs' = Point . flip S.backpermute idxs $ coordinates bs
     in joinAffine bs' tns'


--- CSV ---


data PopulationParameterCounts = PopulationParameterCounts
    { binCentre :: Double
    , parameterCount :: Int }
    deriving (Generic, Show)

instance FromNamedRecord PopulationParameterCounts
instance ToNamedRecord PopulationParameterCounts
instance DefaultOrdered PopulationParameterCounts

populationParameterHistogram
    :: KnownNat k
    => Int
    -> Mean #> Natural # Neurons k <* VonMises
    -> [[PopulationParameterCounts]]
populationParameterHistogram nbns lkl =
    let (nz,nxs) = splitVonMisesPopulationEncoder True lkl
        gns = listCoordinates $ toSource nz
        (mus,kps) = unzip $ S.toPair . coordinates . toSource <$> S.toList nxs
     in do
         prms <- [gns,mus,kps]
--         let (mu,vr) = estimateMeanVariance prms
--             mbnds = if vr < 0.1
--                        then Just (mu-0.5,mu+0.5)
--                        else Nothing
         let (bns,[cnts]) = histogram nbns Nothing [prms]
         return $ zipWith PopulationParameterCounts bns cnts





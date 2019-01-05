{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
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
    -- ** Partition Functions
    , conditionalIPLogPartitionFunction
    , affineConditionalIPLogPartitionFunction
    -- ** Statistics
    , fisherInformation
    , averageLogFisherInformation
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Types ---


--- CSV ---

prjnm :: String
prjnm = "neural-data"

--- Inference ---


-- Test this to see if I can just replace the application with *<.< on the
-- affine transformation.
affineConditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Response k
    -> Double
{-# INLINE affineConditionalIPLogPartitionFunction #-}
affineConditionalIPLogPartitionFunction lkl rprms z =
     potential $ z *<.< snd (splitAffine lkl) <-> rprms

conditionalIPLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
{-# INLINE conditionalIPLogPartitionFunction #-}
conditionalIPLogPartitionFunction lkl z =
    let sz = sufficientStatistic z
        logupst x = sz <.> (snd (splitAffine lkl) >.>* x) - potential (lkl >.>* x) - log (2*pi)
     in logIntegralExp 1e-9 logupst 0 (2*pi) (tail $ range 0 (2*pi) 100)

-- Under the assumption of a flat prior
linearDecoderDivergence
    :: KnownNat k
    => Mean #> Natural # VonMises <* Neurons k
    -> Mean #> Natural # Neurons k <* VonMises
    -> Double -- ^ Log partition function of true posterior
    -> Response k
    -> Double
{-# INLINE linearDecoderDivergence #-}
linearDecoderDivergence dcd lkl nrm z =
    let logpst x = sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x) - potential (lkl >.>* x) - log (2*pi) - nrm
        pst = exp . logpst
        logdcd x = log $ density (dcd >.>* z) x
        dv0 x = pst x * (logpst x - logdcd x)
     in fst $ integrate 1e-9 dv0 0 (2*pi)

getFittedIPLikelihood
    :: String
    -> String
    -> IO (NatNumber,[Double])
{-# INLINE getFittedIPLikelihood #-}
getFittedIPLikelihood expnm dst = do
    (k,xs) <- read <$> goalReadDataset (Experiment prjnm expnm) dst
    return (k,xs)

strengthenIPLikelihood
    :: KnownNat k
    => [Double]
    -> Mean #> Natural # Neurons k <* VonMises
{-# INLINE strengthenIPLikelihood #-}
strengthenIPLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---


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
{-# INLINE fitIPLikelihood #-}
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
{-# INLINE fitLinearDecoder #-}
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

subIPLikelihood
    :: forall k m . (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k + m) <* VonMises
    ->  Mean #> Natural # Neurons k <* VonMises
{-# INLINE subIPLikelihood #-}
subIPLikelihood ppc =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . S.take . S.toRows $ toMatrix tns
        bs' = S.take $ coordinates bs
     in joinAffine (Point bs') tns'

subsampleIPLikelihood
    :: (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k+m) <* VonMises
    -> S.Vector k Int
    ->  Mean #> Natural # Neurons k <* VonMises
{-# INLINE subsampleIPLikelihood #-}
subsampleIPLikelihood ppc idxs =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . flip S.backpermute idxs . S.toRows $ toMatrix tns
        bs' = Point . flip S.backpermute idxs $ coordinates bs
     in joinAffine bs' tns'

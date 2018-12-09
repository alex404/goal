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
    -- ** Partition Functions
    , conditionalIPLogPartitionFunction
    , affineConditionalIPLogPartitionFunction
    -- ** Inference
    , numericalVonMisesPPCPosterior
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


nstms :: Int
nstms = 100

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) (nstms+1)

errbnd :: Double
errbnd = 1e-12

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

getFittedIPLikelihood
    :: String
    -> Dataset
    -> IO (Int,[Double])
getFittedIPLikelihood expnm dst = do
    (k,xs) <- read <$> goalReadDataset prjnm expnm dst
    return (k,xs)

strengthenIPLikelihood
    :: KnownNat k
    => [Double]
    -> Mean #> Natural # Neurons k <* VonMises
strengthenIPLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---


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

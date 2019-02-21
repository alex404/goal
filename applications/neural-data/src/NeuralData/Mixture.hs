{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    Arrows,
    TypeApplications,
    TypeOperators
    #-}

module NeuralData.Mixture
    ( -- * Mixture
      fitMixtureLikelihood
    , getFittedMixtureLikelihood
    , strengthenMixtureLikelihood
    , randomMixtureLikelihood
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import NeuralData
import NeuralData.VonMises

--- Types ---


--- Inference ---


getFittedMixtureLikelihood
    :: String
    -> String
    -> IO (NatNumber,NatNumber,[Double])
getFittedMixtureLikelihood expnm dst = do
    (k,n,xs) <- read . fromJust <$> goalReadDataset (Experiment prjnm expnm) dst
    return (k,n,xs)

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs

randomMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Source # LogNormal
    -> Source # VonMises
    -> Source # LogNormal
    -> Random r (Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises)
randomMixtureLikelihood rmxs sgns sprf sprcs = do
    mxs <- samplePoint rmxs
    gns <- S.replicateM $ randomGains sgns
    tcs <- randomTuningCurves sprf sprcs
    let nctgl = toNatural . Point @ Source $ S.init mxs
    return $ vonMisesMixturePopulationEncoder True nctgl (S.map toNatural gns) tcs


--- Analysis ---


fitMixtureLikelihood
    :: forall r k n . (KnownNat k, KnownNat n)
    => Double -- ^ Learning Rate
    -> Int -- ^ Batch size
    -> Int -- ^ Number of epochs
    -> Natural # Dirichlet (n+1) -- ^ Initial Mixture Weights Distribution
    -> Natural # LogNormal -- ^ Initial Precision Distribution
    -> [(Response k,Double)]
    -> Random r [Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises]
fitMixtureLikelihood eps nbtch nepchs rmxs rkp zxs = do
    kps <- S.replicateM $ samplePoint rkp
    mxs <- samplePoint rmxs
    let (zs,xs) = unzip zxs
        mus = S.generate $ \fnt ->
            let zis = fromIntegral . (`B.index` fnt) <$> zs
             in weightedCircularAverage $ zip zis xs
        sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps mus
    let gns0 = transition . sufficientStatisticT $ fst <$> zxs
        mtx = snd . splitAffine $ vonMisesPopulationEncoder True (Right gns0) sps
    gnss' <- S.replicateM $ Point <$> S.replicateM (uniformR (-10,0))
    let gnss = S.map (gns0 <+>) gnss'
        nctgl = toNatural . Point @ Source $ S.init mxs
        mxmdl = buildMixtureModel gnss nctgl
        mxlkl0 = joinBottomSubLinear mxmdl mtx
        grdcrc = loopCircuit' mxlkl0 $ proc (zxs',mxlkl) -> do
            let (zs',xs') = unzip zxs'
                dmxlkl = vanillaGradient
                    $ mixtureStochasticConditionalCrossEntropyDifferential xs' zs' mxlkl
            mxlkl' <- gradientCircuit eps defaultAdamPursuit -< joinTangentPair mxlkl dmxlkl
            returnA -< mxlkl'
    streamCircuit grdcrc . take nepchs . breakEvery nbtch $ cycle zxs

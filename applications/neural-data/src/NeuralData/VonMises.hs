{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    TypeOperators,
    BangPatterns,
    DeriveGeneric,
    TypeApplications
    #-}

module NeuralData.VonMises
    ( -- * Parsing
      getFittedIPLikelihood
    , strengthenIPLikelihood
    -- * Indexing
    , subIPLikelihood
    , subsampleIPLikelihood
    -- * Fitting
    , fitIPLikelihood
    , fitLinearDecoder
    -- * Generating
    , randomLikelihood
    -- * Algorithms
    , linearDecoderDivergence
    , fisherInformation
    , averageLogFisherInformation
    -- * Analyses
    , analyzeTuningCurves
    , populationParameters
    -- ** Divergence Estimation
    , estimateInformations
    , informationSubsamplingAnalysis
    , informationResamplingAnalysis
    , informationsToRatios
    , normalInformationStatistics
    , histogramInformationStatistics
    -- ** CSV
    , Informations (Informations)
    , InformationRatios (InformationRatios)
    , ParameterCounts (ParameterCounts)
    , ParameterDensities (ParameterDensities)
    , ParameterDensityParameters (ParameterDensityParameters)
    , NormalInformations (NormalInformations)
    , InformationCounts (InformationCounts)
    ) where


--- Imports ---


-- Goal --

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Data.List as L


--- Types ---


--- Generating ---

randomPreferredStimuli :: KnownNat k => Random r (S.Vector k Double)
randomPreferredStimuli = S.replicateM $ uniformR (mnx,mxx)

randomPrecisions :: KnownNat k => Source # LogNormal -> Random r (S.Vector k Double)
randomPrecisions = S.replicateM . samplePoint

randomGains :: KnownNat k => Source # LogNormal -> Random r (Source # Neurons k)
randomGains = initialize

randomTuningCurves :: KnownNat k => Source # LogNormal -> Random r (S.Vector k (Natural # VonMises))
randomTuningCurves sx = do
    mus <- randomPreferredStimuli
    kps <- randomPrecisions sx
    let mukps = S.zipWith S.doubleton mus kps
    return $ S.map (toNatural . Point @ Source) mukps

randomLikelihood
    :: KnownNat k
    => Source # LogNormal
    -> Source # LogNormal
    -> Random r (Mean #> Natural # Neurons k <* VonMises)
randomLikelihood sxgn sxprc = do
    gns <- randomGains sxgn
    tcs <- randomTuningCurves sxprc
    return $ vonMisesPopulationEncoder True (Right $ toNatural gns) tcs


--- Inference ---


-- Under the assumption of a flat prior
linearDecoderDivergence
    :: KnownNat k
    => Mean #> Natural # VonMises <* Neurons k
    -> (Double -> Double) -- ^ True Density
    -> Response k
    -> Double
linearDecoderDivergence dcd trudns z =
    let dcddns = density (dcd >.>* z)
        dv0 x = trudns x * log (trudns x / dcddns x)
     in fst $ integrate 1e-2 dv0 mnx mxx

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

-- | Returns x axis samples, and then y axis sum of tuning curves, rectification
-- curve fit, and individual tuning curves.
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

liePotential :: Natural # VonMises -> Double
liePotential nvm =
    logIntegralExp 1e-6  (unnormalizedLogDensity nvm) 0 (2*pi) (range 0 (2*pi) 100)

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
    average $ log . (/(2*pi*exp 1)) . fisherInformation ppc <$> tail (range 0 (2*pi) 101)

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


data ParameterCounts = ParameterCounts
    { binCentre :: Double
    , parameterCount :: Int
    , parameterAverage :: Double }
    deriving (Show,Generic)

instance ToNamedRecord ParameterCounts where
    toNamedRecord = goalCSVNamer

instance DefaultOrdered ParameterCounts where
    headerOrder = goalCSVOrder

data ParameterDensities = ParameterDensities
    { parameterValue :: Double
    , parameterDensity :: Double }
    deriving (Show,Generic)

instance ToNamedRecord ParameterDensities where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered ParameterDensities where
    headerOrder = goalCSVOrder

data ParameterDensityParameters = ParameterDensityParameters
    { parameterMean :: Double
    , parameterShape :: Double }
    deriving (Show,Generic)

instance ToNamedRecord ParameterDensityParameters where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered ParameterDensityParameters where
    headerOrder = goalCSVOrder

data Informations = Informations
    { mutualInformation :: Double
    , linearDivergence :: Double
    , affineDivergence :: Double
    , decoderDivergence :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord Informations where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered Informations where
    headerOrder = goalCSVOrder

data InformationRatios = InformationRatios
    { linearDivergenceRatio :: Double
    , affineDivergenceRatio :: Double
    , decoderDivergenceRatio :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord InformationRatios where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered InformationRatios where
    headerOrder = goalCSVOrder

data InformationCounts = InformationCounts
    { informationValue :: Double
    , linearDivergenceRatioDensity :: Double
    , affineDivergenceRatioDensity :: Double
    , decoderRatioDensity :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord InformationCounts where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered InformationCounts where
    headerOrder = goalCSVOrder

data NormalInformations = NormalInformations
    { mutualInformationMean :: Double
    , mutualInformationSD :: Double
    , linearDivergenceMean :: Double
    , linearDivergenceSD :: Double
    , linearDivergenceRatioMean :: Double
    , linearDivergenceRatioSD :: Double
    , affineDivergenceMean :: Double
    , affineDivergenceSD :: Double
    , affineDivergenceRatioMean :: Double
    , affineDivergenceRatioSD :: Double
    , decoderDivergenceMean :: Maybe Double
    , decoderDivergenceSD :: Maybe Double
    , meanDecoderDivergenceRatio :: Maybe Double
    , sdDecoderDivergenceRatio :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord NormalInformations where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered NormalInformations where
    headerOrder = goalCSVOrder


--- Statistics ---


informationsToRatios :: Informations -> InformationRatios
informationsToRatios (Informations mi lndvg affdvg mdcddvg) =
    InformationRatios (lndvg/mi) (affdvg/mi) ((/mi) <$> mdcddvg)

populationParameters
    :: KnownNat k
    => Int
    -> Mean #> Natural # Neurons k <* VonMises
    -> [ ( [ParameterCounts]
         , [ParameterDensities]
         , ParameterDensityParameters ) ]
populationParameters nbns lkl =
    let (nz,nxs) = splitVonMisesPopulationEncoder True lkl
        gns = listCoordinates $ toSource nz
        (mus,kps) = unzip $ S.toPair . coordinates . toSource <$> S.toList nxs
     in do
         (bl,prms) <- zip [False,True,False] [gns,mus,kps]
         let (bns,[cnts],[wghts]) = histograms nbns Nothing [prms]
             dx = head (tail bns) - head bns
             (ppds,dprms) = if bl
               then let backprop vm' = joinTangentPair vm' $ stochasticCrossEntropyDifferential prms vm'
                        vm0 = Point $ S.doubleton 0.01 0.01
                        vm :: Natural # VonMises
                        vm = vanillaGradientSequence backprop (-0.1) defaultAdamPursuit vm0 !! 500
                        xs = range mnx mxx 1000
                        dnss = density vm <$> xs
                        (mu,prcs) = S.toPair . coordinates $ toSource vm
                     in ( zipWith ParameterDensities xs dnss
                        , ParameterDensityParameters mu prcs )
               else let lgnrm :: Natural # LogNormal
                        lgnrm = mle $ filter (/= 0) prms
                        xs = range 0 (last bns + dx/2) 1000
                        dnss = density lgnrm <$> xs
                        (mu,sd) = S.toPair . coordinates $ toSource lgnrm
                     in ( zipWith ParameterDensities xs dnss
                        , ParameterDensityParameters mu sd )
         return (zipWith3 ParameterCounts bns cnts wghts,ppds,dprms)


--- Divergence Estimations ---


histogramInformationStatistics
    :: Int
    -> [Informations]
    -> [InformationCounts]
histogramInformationStatistics nbns ppcinfs =
    let (mis,lndvgs,affdvgs,mdcddvgs) = L.unzip4 [ (mi,lndvg,affdvg,mdcddvg)
          | Informations mi lndvg affdvg mdcddvg <- ppcinfs ]
        dcddvgs = if null (head mdcddvgs)
                      then []
                      else fromJust <$> mdcddvgs
        lnrtos = zipWith (/) lndvgs mis
        affrtos = zipWith (/) affdvgs mis
        dcdrtos = zipWith (/) dcddvgs mis
        tracer = concat
            [ "\nMimimal MI: ", show $ minimum mis, "; Maximal MI: ", show $ maximum mis
            , "; Minimal Affine KL: ", show $ minimum affdvgs
            , "; Maximal Affine KL: ", show $ maximum affdvgs, "\n" ]
        filterFun x = not (isInfinite x) && not (isNaN x)
        (bns,_,[lndns,affdns,dcddns0])
          = trace tracer . histograms nbns Nothing $ filter filterFun . map log <$> [lnrtos, affrtos,dcdrtos]
        dcddns = if null dcddns0 then repeat Nothing else Just <$> dcddns0
     in L.zipWith4 InformationCounts (exp <$> bns) lndns affdns dcddns

normalInformationStatistics
    :: [Informations]
    -> NormalInformations
normalInformationStatistics ppcinfs =
    let dvgss = do
            Informations mi lndvg affdvg mdcddvg <- ppcinfs
            return (mi,lndvg,lndvg/mi,affdvg,affdvg/mi,mdcddvg,(/mi) <$> mdcddvg)
        (mis,lndvgs,lnrtos,affdvgs,affrtos,mdcddvgs,mdcdrtos) = L.unzip7 dvgss
        [ (mimu,misd),(lnmu,lnsd),(lnrtomu,lnrtosd),(affmu,affsd),(affrtomu,affrtosd)]
            = meanSDInliers <$> [mis,lndvgs,lnrtos,affdvgs,affrtos]
        (mdcdmu,mdcdsd,mdcdrtomu,mdcdrtosd) =
            if isNothing (head mdcddvgs)
               then (Nothing,Nothing,Nothing,Nothing)
               else let (dvgmu,dvgsd) = meanSDInliers $ fromJust <$> mdcddvgs
                        (rtomu,rtosd) = meanSDInliers $ fromJust <$> mdcdrtos
                     in (Just dvgmu, Just dvgsd, Just rtomu, Just rtosd)
     in NormalInformations
        mimu misd lnmu lnsd lnrtomu lnrtosd affmu affsd affrtomu affrtosd mdcdmu mdcdsd mdcdrtomu mdcdrtosd


informationResamplingAnalysis
    :: forall k r . KnownNat k
    => Int -- ^ Number of regression/rectification samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Int -- ^ Number of population samples
    -> Source # LogNormal -- ^ Gain model
    -> Source # LogNormal -- ^ Precision model
    -> Proxy k -- ^ Population size
    -> Random r [Informations] -- ^ Divergence Statistics
{-# INLINE informationResamplingAnalysis #-}
informationResamplingAnalysis nrct ncntr nmcmc mndcd npop rgns rprcs _ = do
    (lkls :: [Mean #> Natural # Neurons k <* VonMises])
        <- replicateM npop $ randomLikelihood rgns rprcs
    mapM (estimateInformations nrct ncntr nmcmc mndcd) lkls

informationSubsamplingAnalysis
    :: forall k m r . (KnownNat k, KnownNat m)
    => Int -- ^ Number of regression/rectification samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Int -- ^ Number of subpopulation samples
    -> Mean #> Natural # Neurons (k+m+1) <* VonMises -- ^ Complete likelihood
    -> Proxy k -- ^ Subpopulation size
    -> Random r [Informations] -- ^ Divergence Statistics
{-# INLINE informationSubsamplingAnalysis #-}
informationSubsamplingAnalysis nrct ncntr nmcmc mndcd nsub lkl _ = do
    lkls <- replicateM nsub $ do
            (idxs :: B.Vector (k+1) Int) <- generateIndices (Proxy @ (k+m+1))
            return . subsampleIPLikelihood lkl $ G.convert idxs
    mapM (estimateInformations nrct ncntr nmcmc mndcd) lkls

estimateInformations
    :: forall k r . KnownNat k
    => Int -- ^ Number of regression/rectification samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Mean #> Natural # Neurons k <* VonMises -- ^ Complete likelihood
    -> Random r Informations -- ^ Divergence Statistics
{-# INLINE estimateInformations #-}
estimateInformations nrct ncntr nmcmc mndcd lkl = do
    let [rctsmps,cntrsmps,mcmcsmps] = tail . range mnx mxx . (+1) <$> [nrct,ncntr,nmcmc]
        mdcdsmps = tail . range mnx mxx . (+1) <$> mndcd
    mdcd <- case mdcdsmps of
        Just smps -> Just <$> fitLinearDecoder lkl smps
        Nothing -> return Nothing
    estimateConditionalInformations mdcd rctsmps cntrsmps nmcmc mcmcsmps lkl

-- Assumes a uniform prior over stimuli
estimateConditionalInformations
    :: KnownNat k
    => Maybe (Mean #> Natural # VonMises <* Neurons k)
    -> Sample VonMises
    -> Sample VonMises
    -> Int
    -> Sample VonMises
    -> Mean #> Natural # Neurons k <* VonMises
    -> Random r Informations
{-# INLINE estimateConditionalInformations #-}
estimateConditionalInformations mdcd rctsmps cntrsmps nmcmc mcmcsmps lkl = do
    let (rho0,rprms) = regressRectificationParameters lkl rctsmps
    (truprt0,ptnl0,lnprt0,affprt0,mdcddvg0)
        <- foldM (informationsFolder mdcd cntrsmps lkl rprms) (0,0,0,0,Just 0) mcmcsmps
    let k' = fromIntegral nmcmc
        (truprt,ptnl,lnprt,affprt,mdcddvg)
          = (truprt0/k',ptnl0/k',lnprt0/k',affprt0/k',(/k') <$> mdcddvg0)
        !lndvg = lnprt - truprt - rho0
        !affdvg = affprt - truprt - rho0
        !mi = ptnl - truprt - rho0
    return $ Informations mi lndvg affdvg mdcddvg

informationsFolder
    :: KnownNat k
    => Maybe (Mean #> Natural # VonMises <* Neurons k)
    -> Sample VonMises -- ^ centering samples
    -> Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> (Double,Double,Double,Double,Maybe Double)
    -> SamplePoint VonMises
    -> Random r (Double,Double,Double,Double,Maybe Double)
{-# INLINE informationsFolder #-}
informationsFolder mdcd cntrsmps lkl rprms (truprt,ptnl,lnprt,affprt,mdcddvg) x = do
    z <- samplePoint $ lkl >.>* x
    let (dns,truprt') = numericalRecursiveBayesianInference 1e-6 mnx mxx cntrsmps [lkl] [z] (const 1)
        lnprt' = liePotential . fromOneHarmonium $ rectifiedBayesRule zero lkl z zero
        affprt' = liePotential . fromOneHarmonium $ rectifiedBayesRule rprms lkl z zero
        ptnl' = sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x)
        mdcddvg' = do
            dcd <- mdcd
            dcddvg <- mdcddvg
            let dcddvg' = linearDecoderDivergence dcd dns z
            return $ dcddvg + dcddvg'
    return (truprt + truprt',ptnl + ptnl',lnprt + lnprt',affprt + affprt', mdcddvg')

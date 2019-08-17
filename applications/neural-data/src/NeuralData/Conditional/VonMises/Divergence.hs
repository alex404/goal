{-# LANGUAGE DeriveGeneric,TypeApplications #-}

module NeuralData.Conditional.VonMises.Divergence
    ( -- * Parsing
      strengthenIPLikelihood
    -- * Indexing
    , subIPLikelihood
    , subsampleIPLikelihood
    -- * Fitting
    , fitLinearDecoder
    -- * Generating
    , randomGains
    , randomTuningCurves
    , randomIPLikelihood
    -- * Algorithms
    , linearDecoderDivergence
    -- * Analyses
    , analyzeTuningCurves
    , populationParameters
    -- ** Divergence Estimation
    , informationResamplingAnalysis
    , informationSubsamplingAnalysis
    , estimateInformations
    , informationsToRatios
    , normalInformationStatistics
    , logNormalInformationStatistics
    , histogramInformationStatistics
    -- ** CSV
    , Informations (Informations)
    , InformationRatios (InformationRatios)
    , ParameterCounts (ParameterCounts)
    , ParameterDistributionFit (ParameterDistributionFit)
    , NormalInformations (NormalInformations)
    , LogNormalInformations (LogNormalInformations)
    , InformationCounts (InformationCounts)
    ) where


--- Imports ---


-- Goal --

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.List as L


--- Types ---


--- Generating ---

randomGains
    :: KnownNat k => Source # LogNormal -> Random r (Source # Neurons k)
randomGains = initialize

randomPreferredStimuli
    :: KnownNat k => Source # VonMises -> Random r (S.Vector k Double)
randomPreferredStimuli = S.replicateM . samplePoint

randomPrecisions
    :: KnownNat k => Source # LogNormal -> Random r (S.Vector k Double)
randomPrecisions = S.replicateM . samplePoint

randomTuningCurves
    :: KnownNat k
    => Source # VonMises
    -> Source # LogNormal
    -> Random r (S.Vector k (Natural # VonMises))
randomTuningCurves sprf sprcs = do
    mus <- randomPreferredStimuli sprf
    kps <- randomPrecisions sprcs
    let mukps = S.zipWith S.doubleton mus kps
    return $ S.map (toNatural . Point @ Source) mukps

randomIPLikelihood
    :: KnownNat k
    => Source # LogNormal
    -> Source # VonMises
    -> Source # LogNormal
    -> Random r (Natural #> Neurons k <* VonMises)
randomIPLikelihood sgns sprf sprcs = do
    gns <- randomGains sgns
    tcs <- randomTuningCurves sprf sprcs
    return $ joinVonMisesPopulationEncoder (Right $ toNatural gns) tcs


--- Inference ---


-- Under the assumption of a flat prior
linearDecoderDivergence
    :: KnownNat k
    => Natural #> VonMises <* Neurons k
    -> (Double -> Double) -- ^ True Density
    -> Response k
    -> Double
linearDecoderDivergence dcd trudns z =
    let dcddns = density (dcd >.>* z)
        dv0 x = trudns x * log (trudns x / dcddns x)
     in fst $ integrate 1e-2 dv0 mnx mxx

strengthenIPLikelihood
    :: KnownNat k
    => [Double]
    -> Natural #> Neurons k <* VonMises
strengthenIPLikelihood xs = Point . fromJust $ S.fromList xs


--- Analysis ---

-- | Returns x axis samples, and then y axis sum of tuning curves, conjugation
-- curve fit, and individual tuning curves.
analyzeTuningCurves
    :: forall k . KnownNat k
    => Sample VonMises
    -> Natural #> Neurons k <* VonMises
    -> [[Double]]
analyzeTuningCurves xsmps lkl =
    let nzs = lkl >$>* xsmps
        tcss = listCoordinates . toMean <$> nzs
        stcs = potential <$> nzs
        (rho0,rprms) = regressConjugationParameters lkl xsmps
        rcrv = conjugationCurve rho0 rprms xsmps
        mxtcs = maximum <$> tcss
     in zipWith (++) (L.transpose (xsmps:stcs:rcrv:[mxtcs])) tcss

liePotential :: Natural # VonMises -> Double
liePotential nvm = logIntegralExp 1e-6  (unnormalizedLogDensity nvm) 0 (2*pi) (range 0 (2*pi) 100)

-- NB: Actually affine, not linear
fitLinearDecoder
    :: forall s k . KnownNat k
    => Natural #> Neurons k <* VonMises
    -> Sample VonMises
    -> Random s (Natural #> VonMises <* Neurons k)
fitLinearDecoder lkl xs = do
    zs <- mapM samplePoint (lkl >$>* xs)
    let eps = 0.1
        nepchs = 500
        sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        nxz = transpose . fromRows $ S.map toNatural sps
        nx = Point $ S.fromTuple (0,0.5)
        aff0 = joinAffine nx nxz
        backprop = conditionalLogLikelihoodDifferential (zip xs zs)
    return (vanillaGradientSequence backprop eps defaultAdamPursuit aff0 !! nepchs)

subIPLikelihood
    :: forall k m . (KnownNat k, KnownNat m)
    => Natural #> Neurons (k + m) <* VonMises
    ->  Natural #> Neurons k <* VonMises
subIPLikelihood ppc =
    let (bs,tns) = splitAffine ppc
        tns' = fromMatrix . S.fromRows . S.take . S.toRows $ toMatrix tns
        bs' = S.take $ coordinates bs
     in joinAffine (Point bs') tns'

subsampleIPLikelihood
    :: (KnownNat k, KnownNat m)
    => Natural #> Neurons (k+m) <* VonMises
    -> S.Vector k Int
    -> Natural #> Neurons k <* VonMises
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

data ParameterDistributionFit = ParameterDistributionFit
    { parameterValue :: Double
    , parameterDensity :: Double }
    deriving (Show,Generic)

instance ToNamedRecord ParameterDistributionFit where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered ParameterDistributionFit where
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
    { meanMutualInformation :: Double
    , sdMutualInformation :: Double
    , meanLinearDivergence :: Double
    , sdLinearDivergence :: Double
    , meanLinearDivergenceRatio :: Double
    , sdLinearDivergenceRatio :: Double
    , meanAffineDivergence :: Double
    , sdAffineDivergence :: Double
    , meanAffineDivergenceRatio :: Double
    , sdAffineDivergenceRatio :: Double
    , meanDecoderDivergence :: Maybe Double
    , sdDecoderDivergence :: Maybe Double
    , meanDecoderDivergenceRatio :: Maybe Double
    , sdDecoderDivergenceRatio :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord NormalInformations where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered NormalInformations where
    headerOrder = goalCSVOrder

data LogNormalInformations = LogNormalInformations
    { logMeanHomogeneousDivergenceRatio :: Double
    , lowerLogBoundHomogeneousDivergenceRatio :: Double
    , upperLogBoundHomogeneousDivergenceRatio :: Double
    , logMeanConjugateDivergenceRatio :: Double
    , lowerLogBoundConjugateDivergenceRatio :: Double
    , upperLogBoundConjugateDivergenceRatio :: Double
    , logMeanDecoderDivergenceRatio :: Maybe Double
    , lowerLogBoundDecoderDivergenceRatio :: Maybe Double
    , upperLogBoundDecoderDivergenceRatio :: Maybe Double }
    deriving (Show,Generic)

instance ToNamedRecord LogNormalInformations where
    toNamedRecord = goalCSVNamer
instance DefaultOrdered LogNormalInformations where
    headerOrder = goalCSVOrder


--- Statistics ---


informationsToRatios :: Informations -> InformationRatios
informationsToRatios (Informations mi lndvg affdvg mdcddvg) =
    InformationRatios (lndvg/mi) (affdvg/mi) ((/mi) <$> mdcddvg)

populationParameters
    :: KnownNat k
    => Int
    -> Natural #> Neurons k <* VonMises
    -> ( [([ParameterCounts], [ParameterDistributionFit])]
       , (Source # LogNormal, Source # VonMises, Source # LogNormal) )
populationParameters nbns lkl =
    let (nz,nxs) = splitVonMisesPopulationEncoder lkl
        gns = listCoordinates $ toSource nz
        (mus,kps) = unzip $ S.toPair . coordinates . toSource <$> S.toList nxs
        (pcntss,pftdnss,[gncs,prfcs,prcscs]) = unzip3 $ do
          (bl,prms) <- zip [False,True,False] [gns,mus,kps]
          let (bns,[cnts],[wghts]) = histograms nbns Nothing [prms]
              dx = head (tail bns) - head bns
              (ppds,dprms) = if bl
                  then let backprop = logLikelihoodDifferential prms
                           vm0 = Point $ S.doubleton 0.01 0.01
                           vm :: Natural # VonMises
                           vm = vanillaGradientSequence backprop (-0.1) defaultAdamPursuit vm0 !! 500
                           xs = range mnx mxx 1000
                           dnss = density vm <$> xs
                        in ( zipWith ParameterDistributionFit xs dnss, coordinates $ toSource vm )
                  else let lgnrm :: Natural # LogNormal
                           lgnrm = mle $ filter (/= 0) prms
                           xs = range 0 (last bns + dx/2) 1000
                           dnss = density lgnrm <$> xs
                        in ( zipWith ParameterDistributionFit xs dnss, coordinates $ toSource lgnrm )
          return (zipWith3 ParameterCounts bns cnts wghts,ppds,dprms)
     in (zip pcntss pftdnss, (Point gncs, Point prfcs, Point prcscs))


--- Divergence Estimations ---


filterFun :: Double -> Bool
filterFun x = not (isInfinite x) && not (isNaN x)

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

logNormalInformationStatistics
    :: [Informations]
    -> LogNormalInformations
logNormalInformationStatistics ppcinfs =
    let dvgss = do
            Informations mi lndvg affdvg mdcddvg <- ppcinfs
            return (lndvg/mi,affdvg/mi,(/mi) <$> mdcddvg)
        (lnrtos,affrtos,mdcdrtos) = L.unzip3 dvgss
        (lnrtomu,lnrtosd) = meanSDInliers . filter filterFun $ log <$> lnrtos
        (affrtomu,affrtosd) = meanSDInliers . filter filterFun $ log <$> affrtos
        (lnrtolbnd,lnrtoubnd) = (exp $ lnrtomu - lnrtosd, exp $ lnrtomu + lnrtosd)
        (affrtolbnd,affrtoubnd) = (exp $ affrtomu - affrtosd, exp $ affrtomu + affrtosd)
        (mdcdrtomu,mdcdrtolbnd,mdcdrtoubnd) =
            if isNothing (head mdcdrtos)
               then (Nothing,Nothing,Nothing)
               else let (rtomu,rtosd) = meanSDInliers . filter filterFun $ log . fromJust <$> mdcdrtos
                     in (Just $ exp rtomu, Just . exp $ rtomu - rtosd, Just . exp $ rtomu + rtosd)
     in LogNormalInformations
        (exp lnrtomu) lnrtolbnd lnrtoubnd (exp affrtomu) affrtolbnd affrtoubnd mdcdrtomu mdcdrtolbnd mdcdrtoubnd

informationResamplingAnalysis
    :: forall n k r . (KnownNat k, KnownNat n)
    => Int -- ^ Number of regression/conjugation samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Int -- ^ Number of subpopulation samples
    -> Natural #> Neurons k <* VonMises -- ^ Complete likelihood
    -> Proxy n -- ^ Max population size
    -> Random r [Informations] -- ^ Divergence Statistics
informationResamplingAnalysis nrct ncntr nmcmc mndcd nsub lkl prxn = do
    let (sgns,sprfs,sprcs) = snd $ populationParameters 10 lkl
    informationResamplingAnalysis0 nrct ncntr nmcmc mndcd nsub sgns sprfs sprcs prxn

informationResamplingAnalysis0
    :: forall n r . KnownNat n
    => Int -- ^ Number of regression/conjugation samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Int -- ^ Number of population samples
    -> Source # LogNormal -- ^ Gain model
    -> Source # VonMises -- ^ Preferred Stimulus model
    -> Source # LogNormal -- ^ Precision model
    -> Proxy n -- ^ Population size
    -> Random r [Informations] -- ^ Divergence Statistics
{-# INLINE informationResamplingAnalysis0 #-}
informationResamplingAnalysis0 nrct ncntr nmcmc mndcd npop sgns sprf sprcs _ = do
    (lkls :: [Natural #> Neurons (n+1) <* VonMises])
        <- replicateM npop $ randomIPLikelihood sgns sprf sprcs
    mapM (estimateInformations nrct ncntr nmcmc mndcd) lkls

informationSubsamplingAnalysis
    :: forall k m r . (KnownNat k, KnownNat m)
    => Int -- ^ Number of regression/conjugation samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Int -- ^ Number of subpopulation samples
    -> Natural #> Neurons (k+m+1) <* VonMises -- ^ Complete likelihood
    -> Proxy k -- ^ Subpopulation size = k+1
    -> Random r [Informations] -- ^ Divergence Statistics
{-# INLINE informationSubsamplingAnalysis #-}
informationSubsamplingAnalysis nrct ncntr nmcmc mndcd nsub lkl _ = do
    lkls <- replicateM nsub $ do
            (idxs :: S.Vector (k+1) Int) <- generateIndices (Proxy @ (k+m+1))
            return . subsampleIPLikelihood lkl $ G.convert idxs
    mapM (estimateInformations nrct ncntr nmcmc mndcd) lkls

estimateInformations
    :: forall k r . KnownNat k
    => Int -- ^ Number of regression/conjugation samples
    -> Int -- ^ Number of numerical centering samples
    -> Int -- ^ Number of monte carlo integration samples
    -> Maybe Int -- ^ (Maybe) number of linear decoder samples
    -> Natural #> Neurons k <* VonMises -- ^ Complete likelihood
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
    => Maybe (Natural #> VonMises <* Neurons k)
    -> Sample VonMises
    -> Sample VonMises
    -> Int
    -> Sample VonMises
    -> Natural #> Neurons k <* VonMises
    -> Random r Informations
{-# INLINE estimateConditionalInformations #-}
estimateConditionalInformations mdcd rctsmps cntrsmps nmcmc mcmcsmps lkl = do
    let (rho0,rprms) = regressConjugationParameters lkl rctsmps
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
    => Maybe (Natural #> VonMises <* Neurons k)
    -> Sample VonMises -- ^ centering samples
    -> Natural #> Neurons k <* VonMises
    -> Natural # VonMises
    -> (Double,Double,Double,Double,Maybe Double)
    -> SamplePoint VonMises
    -> Random r (Double,Double,Double,Double,Maybe Double)
{-# INLINE informationsFolder #-}
informationsFolder mdcd cntrsmps lkl rprms (truprt,ptnl,lnprt,affprt,mdcddvg) x = do
    z <- samplePoint $ lkl >.>* x
    let (dns,truprt') = numericalRecursiveBayesianInference 1e-6 mnx mxx cntrsmps [lkl] [z] (const 1)
        lnprt' = liePotential . fromOneHarmonium $ conjugatedBayesRule zero lkl z zero
        affprt' = liePotential . fromOneHarmonium $ conjugatedBayesRule rprms lkl z zero
        ptnl' = sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x)
        mdcddvg' = do
            dcd <- mdcd
            dcddvg <- mdcddvg
            let dcddvg' = linearDecoderDivergence dcd dns z
            return $ dcddvg + dcddvg'
    return (truprt + truprt',ptnl + ptnl',lnprt + lnprt',affprt + affprt', mdcddvg')


--- Graveyard  ---


--ppcStimulusDerivatives
--    :: KnownNat k
--    => Natural #> Neurons k <* VonMises
--    -> SamplePoint VonMises
--    -> S.Vector k Double
--ppcStimulusDerivatives ppc x =
--    let fxs = coordinates . toMean $ ppc >.> mx
--        tcs = toRows . snd $ splitAffine ppc
--     in S.zipWith zipper fxs tcs
--    where mx = sufficientStatistic x
--          (cx,sx) = S.toPair $ coordinates mx
--          zipper fx (Point cs) =
--              let (tht1,tht2) = S.toPair cs
--               in fx*(cx * tht2 - sx * tht1)
--
--fisherInformation
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Double
--    -> Double
--fisherInformation ppc x =
--    let fxs2' = S.map square $ ppcStimulusDerivatives ppc x
--        fxs = coordinates . toMean $ ppc >.>* x
--     in S.sum $ S.zipWith (/) fxs2' fxs
--
--averageLogFisherInformation
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Double
--averageLogFisherInformation ppc =
--    average $ log . (/(2*pi*exp 1)) . fisherInformation ppc <$> tail (range 0 (2*pi) 101)

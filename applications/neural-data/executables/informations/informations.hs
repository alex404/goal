{-# LANGUAGE
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    BangPatterns,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData
import NeuralData.VonMises

import Paths_neural_data

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import Data.Semigroup ((<>))
import Data.List


--- CSV ---


data VonMisesInformations = VonMisesInformations
    { meanLinearDivergence :: Double
    , sdLinearDivergence :: Double
    , meanAffineDivergence :: Double
    , sdAffineDivergence :: Double
    , meanDecoderDivergence :: Double
    , sdDecoderDivergence :: Double
    , meanVonMisesMutualInformation :: Double
    , sdVonMisesMutualInformation :: Double
    , meanLinearDivergenceRatio :: Double
    , sdLinearDivergenceRatio :: Double
    , meanAffineDivergenceRatio :: Double
    , sdAffineDivergenceRatio :: Double
    , meanDecoderDivergenceRatio :: Double
    , sdDecoderDivergenceRatio :: Double }
    deriving (Generic, Show)

instance FromNamedRecord VonMisesInformations
instance ToNamedRecord VonMisesInformations
instance DefaultOrdered VonMisesInformations
instance NFData VonMisesInformations



--- Analysis ---

ananm :: String
ananm = "informations"

nstms :: Int
nstms = 1000

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) (nstms+1)

informationsFolder
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Mean #> Natural # VonMises <* Neurons k
    -> (Double,Double,Double,Double,Double)
    -> Double
    -> Random r (Double,Double,Double,Double,Double)
informationsFolder lkl rprms dcd (zpn,zqn,zcn,ptn,dcddvg) x = do
    z <- samplePoint $ lkl >.>* x
    let (dns,zpn') = numericalRecursiveBayesianInference 1e-6 0 (2*pi) 100 [lkl] [z] (const 1)
        zqn' = potential . fromOneHarmonium $ rectifiedBayesRule zero lkl z zero
        zcn' = potential . fromOneHarmonium $ rectifiedBayesRule rprms lkl z zero
        ptn' = sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x)
        dcddvg' = linearDecoderDivergence dcd dns z
    return (zpn + zpn',zqn + zqn',zcn + zcn',ptn + ptn',dcddvg + dcddvg')

-- Assumes a uniform prior over stimuli
estimateInformations
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Mean #> Natural # VonMises <* Neurons k
    -> Random r (Double,Double,Double,Double,Double,Double,Double)
estimateInformations lkl dcd = do
    let (rho0,rprms) = regressRectificationParameters lkl xsmps
    (zpnavg0,zqnavg0,zcnavg0,ptnavg0,dcddvg0) <- foldM (informationsFolder lkl rprms dcd) (0,0,0,0,0) xsmps
    let k' = fromIntegral nstms
        (zpnavg,zqnavg,zcnavg,ptnavg,!dcddvg) = (zpnavg0/k',zqnavg0/k',zcnavg0/k',ptnavg0/k',dcddvg0/k')
        !pq0dvg = zqnavg - zpnavg - rho0
        !pqdvg = zcnavg - zpnavg - rho0
        !mi = ptnavg - zpnavg - rho0
    return (pq0dvg,pqdvg,dcddvg,mi,pq0dvg/mi,pqdvg/mi,dcddvg/mi)

vonMisesInformationsStatistics
    :: forall k m . (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k+m+1) <* VonMises
    -> Int
    -> Proxy k
    -> IO VonMisesInformations
vonMisesInformationsStatistics lkl nsmps _ = do
    (ldvgs,advgs,dcdavgs,mis,nrmldvgs,nrmadvgs,nrmdcdavgs) <- realize . fmap unzip7 . replicateM nsmps $ do
        (idxs :: B.Vector (k+1) Int) <- generateIndices (Proxy @ (k+m+1))
        let sublkl = subsampleIPLikelihood lkl $ G.convert idxs
        dcd <- fitLinearDecoder sublkl xsmps
        estimateInformations sublkl dcd
    let [ (ldvgmu,ldvgsd)
        , (advgmu,advgsd)
        , (dcdmu,dcdsd)
        , (mimu,misd)
        , (nrmldvgmu,nrmldvgsd)
        , (nrmadvgmu,nrmadvgsd)
        , (nrmdcdmu,nrmdcdsd) ] = meanSDInliers <$> [ldvgs,advgs,dcdavgs,mis,nrmldvgs,nrmadvgs,nrmdcdavgs]
        stp = concat ["Step ", show . natValInt $ Proxy @ k
                     , "; ldvgmu: " ++ show ldvgmu
                     , "; advgmu: " ++ show advgmu
                     , "; dcdmu: " ++ show dcdmu
                     , "; mimu: " ++ show mimu ]
    putStrLn stp
    return $ VonMisesInformations
        ldvgmu ldvgsd advgmu advgsd dcdmu dcdsd mimu misd nrmldvgmu nrmldvgsd nrmadvgmu nrmadvgsd nrmdcdmu nrmdcdsd

fitAnalyzeInformations
    :: forall k . KnownNat k
    => Int
    -> [([Int], Double)]
    -> Proxy k
    -> IO [VonMisesInformations]
fitAnalyzeInformations nsmps zxss0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxss0

    lkl <- realize $ fitIPLikelihood zxs

    (alldvgs0 :: B.Vector k VonMisesInformations)
        <- B.generatePM $ vonMisesInformationsStatistics lkl nsmps

    return $ B.toList alldvgs0


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

vminfOpts :: Parser AnalysisOpts
vminfOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 's' <> value 10)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg nsmps) = do

    let expmnt = Experiment prjnm expnm

    dsts <- if null dstarg
               then fromJust <$> goalReadDatasetsCSV expmnt
               else return [dstarg]

    infgpi <- getDataFileName "informations/informations.gpi"

    forM_ dsts $ \dst -> do

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs = case someNatVal k of
                    SomeNat prxk -> fitAnalyzeInformations nsmps zxs prxk

        infs <- rinfs

        let msbexp = (Just $ SubExperiment ananm dst)

        goalWriteNamedAnalysis True expmnt msbexp infs

        runGnuplot expmnt msbexp defaultGnuplotOptions infgpi

--- Main ---


main :: IO ()
main = do

    let opts = info (vminfOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts

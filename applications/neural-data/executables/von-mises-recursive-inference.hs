{-# LANGUAGE
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData
import NeuralData.VonMises

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
ananm = "von-mises-informations"

nstms :: Int
nstms = 100

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
    let zpn' = snd $ numericalRecursiveBayesianInference 1e-12 0 (2*pi) 100 [lkl] [z] (const 1)
        zqn' = potential . fromOneHarmonium $ rectifiedBayesRule zero lkl z zero
        zcn' = potential . fromOneHarmonium $ rectifiedBayesRule rprms lkl z zero
        ptn' = sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x)
        dcddvg' = linearDecoderDivergence dcd lkl zpn' z
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
        (zpnavg,zqnavg,zcnavg,ptnavg,dcddvg) = (zpnavg0/k',zqnavg0/k',zcnavg0/k',ptnavg0/k',dcddvg0/k')
        pq0dvg = zqnavg - zpnavg - rho0
        pqdvg = zcnavg - zpnavg - rho0
        mi = ptnavg - zpnavg - rho0
    return (pq0dvg,pqdvg,dcddvg,mi,pq0dvg/mi,pqdvg/mi,dcddvg/mi)

vonMisesInformationsStatistics
    :: forall k m r . (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k+m+1) <* VonMises
    -> [(Response (k+m+1), Double)]
    -> Int
    -> Proxy k
    -> Random r VonMisesInformations
vonMisesInformationsStatistics lkl zxs nsmps _ = do
    (ldvgs,advgs,dcdavgs,mis,nrmldvgs,nrmadvgs,nrmdcdavgs) <- fmap unzip7 . replicateM nsmps $ do
        (idxs :: B.Vector (k+1) Int) <- generateIndices (Proxy @ (k+m+1))
        let sublkl = subsampleIPLikelihood lkl $ G.convert idxs
            (zs,xs) = unzip zxs
            sbzxs = flip zip xs $ (`B.backpermute` idxs) <$> zs
            dcd = fitLinearDecoder sbzxs
        estimateInformations sublkl dcd
    let [(ldvgmu,ldvgsd),(advgmu,advgsd),(dcdmu,dcdsd),(mimu,misd),(nrmldvgmu,nrmldvgsd),(nrmadvgmu,nrmadvgsd),(nrmdcdmu,nrmdcdsd) ] =
            meanSDInliers <$> [ldvgs,advgs,dcdavgs,mis,nrmldvgs,nrmadvgs,nrmdcdavgs]
        stp = concat ["Step ", show . natValInt $ Proxy @ k
                     , "; ldvgmu: " ++ show ldvgmu
                     , "; advgmu: " ++ show advgmu
                     , "; dcdmu: " ++ show dcdmu
                     , "; mimu: " ++ show mimu ]
    return . trace stp $ VonMisesInformations
        ldvgmu ldvgsd advgmu advgsd dcdmu dcdsd mimu misd nrmldvgmu nrmldvgsd nrmadvgmu nrmadvgsd nrmdcdmu nrmdcdsd

fitAnalyzeInformations
    :: forall k r . KnownNat k
    => Int
    -> [([Int], Double)]
    -> Proxy k
    -> Random r [VonMisesInformations]
fitAnalyzeInformations nsmps zxss0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxss0

    lkl <- fitIPLikelihood zxs

    (alldvgs0 :: B.Vector k VonMisesInformations)
        <- B.generatePM $ vonMisesInformationsStatistics lkl zxs nsmps

    return $ B.toList alldvgs0

--analyzeInformations
--    :: forall k r . KnownNat k
--    => Int
--    -> [Double]
--    -> Proxy k
--    -> Random r [VonMisesInformations]
--analyzeInformations nsmps css _ = undefined
--
--    let zxss :: [(Response k, Double)]
--        zxss = strengthenNeuralData zxss0
--
--    (alldvgs0 :: B.Vector k VonMisesInformations)
--        <- B.generatePM' $ vonMisesInformationsStatistics zxss nsmps
--
--    return $ B.toList alldvgs0


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

    forM_ dsts $ \dst -> do

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs = case someNatVal k of
                    SomeNat prxk -> fitAnalyzeInformations nsmps zxs prxk

        infs <- realize rinfs

        goalWriteNamedAnalysis True expmnt (Just $ SubExperiment ananm dst) infs

--    if take 4 expnm == "true"
--
--       then forM_ dsts $ \dst -> do
--
--                (k,cs) <- getFitPPC expnm dst
--
--                let wghts :: [Double]
--                    infs = withNat k (analyzeInformations cs)
--
--                goalWriteAnalysis prjnm expnm ananm (Just dst) wghts
--                goalAppendAnalysis prjnm expnm ananm (Just dst) stcs
--                mapM_ (goalAppendAnalysis prjnm expnm ananm (Just dst)) tcss
--
--runOpts :: AnalysisOpts -> IO ()
--runOpts (AnalysisOpts expnm dstarg nsmps) = do
--
--    let pth = "projects/" ++ expnm ++ "/analysis/vminf"
--
--    createDirectoryIfMissing True pth
--
--    dsts <- if dstarg == ""
--               then fromJust <$> maybeGetDatasets prjnm expnm
--               else return [Dataset dstarg]
--
--    (kzxss :: [(Int,[([Int], s)])]) <- mapM (getNeuralData expnm) dsts
--    csvss <- realize $ mapM (analyzeInformations nsmps) zxss
--
--    forM_ (zip csvss dsts) $ \(csvs, Dataset dstnm) ->
--        goalWriteAnalysis prjnm expnm "von-mises-informations" (Just dstnm)
--            $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (vminfOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



{-# LANGUAGE
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData

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
    , meanVonMisesMutualInformation :: Double
    , sdVonMisesMutualInformation :: Double
    , meanLinearDivergenceRatio :: Double
    , sdLinearDivergenceRatio :: Double
    , meanAffineDivergenceRatio :: Double
    , sdAffineDivergenceRatio :: Double }
    deriving (Generic, Show)

instance FromNamedRecord VonMisesInformations
instance ToNamedRecord VonMisesInformations
instance DefaultOrdered VonMisesInformations
instance NFData VonMisesInformations



--- Analysis ---

ananm :: String
ananm = "von-mises-informations"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) (nstms+1)

nstms :: Int
nstms = 100


informationsFolder
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> (Double,Double,Double,Double)
    -> Double
    -> Random r (Double,Double,Double,Double)
informationsFolder lkl rprms (zpn,zqn,zcn,ptn) x = do
    let nz = lkl >.>* x
    z <- samplePoint nz
    let zpn' = conditionalIPLogPartitionFunction lkl z
        zqn' = affineConditionalIPLogPartitionFunction lkl zero z
        zcn' = affineConditionalIPLogPartitionFunction lkl rprms z
        ptn' = sufficientStatistic z <.> nz
    return (zpn + zpn',zqn + zqn',zcn + zcn',ptn + ptn')

-- Assumes a uniform prior over stimuli
estimateInformations
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Random r (Double,Double,Double,Double,Double)
estimateInformations lkl = do
    let stcavg = average $ potential <$> lkl >$>* xsmps
        rprms = snd $ regressRectificationParameters lkl xsmps
        rctavg = average [rprms <.> sufficientStatistic x | x <- xsmps]
    (zpnavg0,zqnavg0,zcnavg0,ptnavg0) <- foldM (informationsFolder lkl rprms) (0,0,0,0) xsmps
    let k' = fromIntegral nstms
        (zpnavg,zqnavg,zcnavg,ptnavg) = (zpnavg0/k',zqnavg0/k',zcnavg0/k',ptnavg0/k')
        pq0dvg = zqnavg - zpnavg - stcavg
        pqdvg = zcnavg - zpnavg - stcavg + rctavg
        mi = ptnavg - stcavg - zpnavg
    return (pq0dvg,pqdvg,mi,pq0dvg/mi,pqdvg/mi)

vonMisesInformationsStatistics
    :: forall k m r . (KnownNat k, KnownNat m)
    => Mean #> Natural # Neurons (k+m) <* VonMises
    -> Int
    -> Proxy k
    -> Random r VonMisesInformations
vonMisesInformationsStatistics lkl nsmps _ = do
    (ldvgs,advgs,mis,nrmldvgs,nrmadvgs) <- fmap unzip5 . replicateM nsmps $ do
        (idxs :: B.Vector k Int) <- generateIndices (Proxy @ (k+m))
        let sublkl = subsampleIPLikelihood lkl $ G.convert idxs
        estimateInformations sublkl
    let [(ldvgmu,ldvgsd),(advgmu,advgsd),(mimu,misd),(nrmldvgmu,nrmldvgsd),(nrmadvgmu,nrmadvgsd)] =
            meanSDInliers <$> [ldvgs,advgs,mis,nrmldvgs,nrmadvgs]
        stp = concat ["Step ", show . natValInt $ Proxy @ k
                     , "; ldvgmu: " ++ show ldvgmu
                     , "; advgmu: " ++ show advgmu
                     , "; mimu: " ++ show mimu ]
    return . trace stp $ VonMisesInformations
        ldvgmu ldvgsd advgmu advgsd mimu misd nrmldvgmu nrmldvgsd nrmadvgmu nrmadvgsd

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
        <- B.generatePM $ vonMisesInformationsStatistics lkl nsmps

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

    dsts <- if dstarg == ""
               then fromJust <$> goalReadDatasetsCSV prjnm expnm
               else return [Dataset dstarg]

    forM_ dsts $ \dst -> do

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let foo = case someNatVal k of
                    SomeNat prxk -> fitAnalyzeInformations nsmps zxs prxk


        infs <- realize foo

        goalWriteNamedAnalysis prjnm expnm ananm (Just dst) infs

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



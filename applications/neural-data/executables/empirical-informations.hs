{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import Goal.Plot
import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV
import qualified Data.Map as M

import Data.Semigroup ((<>))


--- CSV ---

data EmpiricalInformations = EmpiricalInformations
    { meanEmpiricalDivergence :: Double
    , sdEmpiricalDivergence :: Double
    , meanEmpiricalDivergenceRatio :: Double
    , sdEmpiricalDivergenceRatio :: Double
    , meanEmpiricalMutualInformation :: Double
    , sdEmpiricalMutualInformation :: Double }
    deriving (Generic, Show)

instance CSV.FromNamedRecord EmpiricalInformations
instance CSV.ToNamedRecord EmpiricalInformations
instance CSV.DefaultOrdered EmpiricalInformations
instance NFData EmpiricalInformations



--- Analysis ---

partialDot :: Mean # Neurons k -> Natural # Neurons k -> Double
partialDot mzs nrts = sum $ do
    (mz,nrt) <- zip (listCoordinates mzs) (listCoordinates nrts)
    guard $ mz > 0
    return $ mz * nrt

empiricalConditionalLogPartitionFunction
    :: KnownNat k
    => [Mean # Neurons k]
    -> [Natural # Neurons k]
    -> Mean # Neurons k
    -> Double
empiricalConditionalLogPartitionFunction rtss nrtss mz =
    logSumExp [ partialDot mz nrts - sum (listCoordinates rts)
      | (rts,nrts) <- zip rtss nrtss ]

empiricalApproximateConditionalLogPartitionFunction
    :: KnownNat k
    => [Natural # Neurons k]
    -> Mean # Neurons k
    -> Double
empiricalApproximateConditionalLogPartitionFunction nrtss mz =
    logSumExp $ partialDot mz <$> nrtss

informationsFolder
    :: KnownNat k
    => [Mean # Neurons k]
    -> [Natural # Neurons k]
    -> (Double,Double,Double)
    -> (Mean # Neurons k, Natural # Neurons k)
    -> Random r (Double,Double,Double)
informationsFolder rtss nrtss (zpn,zqn,ptn) (rts,nrts) = do
    mz <- sufficientStatistic <$> samplePoint rts
    let zpn' = empiricalConditionalLogPartitionFunction rtss nrtss mz
        zqn' = empiricalApproximateConditionalLogPartitionFunction nrtss mz
        ptn' = partialDot mz nrts
    return (zpn + zpn',zqn + zqn',ptn + ptn')


-- Assumes a uniform prior over stimuli
estimateInformations
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> Random r (Double,Double,Double)
estimateInformations lkl = do
    let rtss = M.elems lkl
        k = 10 * length rtss
        nrtss = dualTransition <$> rtss
        cycrts = take k . cycle $ zip rtss nrtss
    (zpnavg0,zqnavg0,ptnavg0) <- foldM (informationsFolder rtss nrtss) (0,0,0) cycrts
    let k' = fromIntegral k
        (zpnavg,zqnavg,ptnavg) = (zpnavg0/k',zqnavg0/k',ptnavg0/k')
        stcavg = average $ sum . listCoordinates <$> rtss
        pqdvg = zqnavg - zpnavg - stcavg
        mi = ptnavg - stcavg - zpnavg + (log . fromIntegral . length $ M.keys lkl)
    return (pqdvg,pqdvg/mi,mi)

empiricalInformationStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => [(Response k,s)]
    -> Int
    -> Proxy k'
    -> Random r EmpiricalInformations
empiricalInformationStatistics zxss n _ = do
    let nzxmp = empiricalTuningCurves $ stimulusResponseMap zxss
    (dvgs,nrmdvgs,mis) <- fmap unzip3 . replicateM n $ do
        (idxs :: B.Vector k' Int) <- generateIndices (Proxy :: Proxy k)
        let sublkl = subSampleEmpiricalTuningCurves nzxmp $ G.convert idxs
        estimateInformations sublkl
    let [(dvgmu,dvgsd),(nrmdvgmu,nrmdvgsd),(mimu,misd)] = meanSDInliers <$> [dvgs,nrmdvgs,mis]
        stp = concat [ "Step ", show . natValInt $ (Proxy :: Proxy k')
                     , "; dvgmu: " ++ show dvgmu
                     , "; nrmdvgmu: " ++ show nrmdvgmu
                     , "; mimu: " ++ show mimu ]
    return . trace stp $ EmpiricalInformations dvgmu dvgsd nrmdvgmu nrmdvgsd mimu misd

analyzeInformation0
    :: forall k s r . (Ord s, Read s, KnownNat k)
    => Int
    -> [([Int], s)]
    -> Proxy k
    -> Random r [EmpiricalInformations]
analyzeInformation0 nsmps zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0

    (alldvgs0 :: B.Vector k EmpiricalInformations)
        <- B.generatePM' $ empiricalInformationStatistics zxss nsmps

    return $ B.toList alldvgs0

analyzeInformation
    :: (Ord x, Read x)
    => Int
    -> [([Int],x)]
    -> Random r [EmpiricalInformations]
analyzeInformation nsmps zxs =
    withNat (getPopulationSize zxs) $ analyzeInformation0 nsmps zxs


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

dvgOpts :: Parser AnalysisOpts
dvgOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 10)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts clcstr dststr nsmps) = do

    let pth = "projects/" ++ clcstr ++ "/analysis/inf"

    let clc = Collection clcstr

    createDirectoryIfMissing True pth

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    csvss <- case clcstr of
               "coen-cagli-2015" -> do
                   (zxss :: [[([Int],Int)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeInformation nsmps) zxss
               "patterson-2013" -> do
                   (zxss :: [[([Int],Double)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeInformation nsmps) zxss
               _ -> error "Invalid project"

    forM_ (zip csvss dsts) $ \(csvs, Dataset dststr') ->
        BS.writeFile (pth ++"/" ++ dststr' ++ ".csv") $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (dvgOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



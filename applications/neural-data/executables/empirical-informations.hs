{-# LANGUAGE TupleSections,DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

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


data Informations = Informations
    { meanDivergence :: Double
    , sdDivergence :: Double
    , meanNormalizedDivergence :: Double
    , sdNormalizedDivergence :: Double
    , meanMutualInformation :: Double
    , sdMutualInformation :: Double }
    deriving (Generic, Show)

instance FromNamedRecord Informations
instance ToNamedRecord Informations
instance DefaultOrdered Informations
instance NFData Informations


--- Analysis ---

partialDot :: Mean # Neurons k -> Natural # Neurons k -> Double
partialDot mzs nrts = sum $ do
    (mz,nrt) <- zip (listCoordinates mzs) (listCoordinates nrts)
    guard $ mz > 0
    return $ mz * nrt

conditionalLogPartitionFunction
    :: KnownNat k
    => [Mean # Neurons k]
    -> [Natural # Neurons k]
    -> Mean # Neurons k
    -> Double
conditionalLogPartitionFunction rtss nrtss mz =
    logSumExp [ partialDot mz nrts - sum (listCoordinates rts)
      | (rts,nrts) <- zip rtss nrtss ]

approximateConditionalLogPartitionFunction
    :: KnownNat k
    => [Natural # Neurons k]
    -> Mean # Neurons k
    -> Double
approximateConditionalLogPartitionFunction nrtss mz =
    logSumExp $ partialDot mz <$> nrtss

-- Assumes a uniform prior over stimuli
estimateInformations
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> Random r (Double,Double,Double)
estimateInformations lkl = do
    let rtss = M.elems lkl
        nrtss = dualTransition <$> rtss
        stcavg = average $ sum . listCoordinates <$> rtss
    zss <- mapM (sample 10) rtss
    let nrtmzs = concat [ (nrts,) . sufficientStatistic <$> zs0 | (nrts,zs0) <- zip nrtss zss ]
        mzs = snd <$> nrtmzs
        zpnavg = average $ conditionalLogPartitionFunction rtss nrtss <$> mzs
        zqnavg = average $ approximateConditionalLogPartitionFunction nrtss <$> mzs
        ptnavg = average [ partialDot mz nrts | (nrts,mz) <- nrtmzs ]
        pqdvg = zqnavg - zpnavg - stcavg
        mi = ptnavg - stcavg - zpnavg + (log . fromIntegral . length $ M.keys lkl)
    return (pqdvg,mi,pqdvg/mi)

empiricalDivergenceStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => [(Response k,s)]
    -> Int
    -> Proxy k'
    -> Random r Informations
empiricalDivergenceStatistics zxss n _ = do
    let nzxmp = empiricalTuningCurves $ stimulusResponseMap zxss
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let sidxss = G.convert <$> idxss
        sublkls = subSampleEmpiricalTuningCurves nzxmp <$> sidxss
    (dvgs,mis,nrmdvgs) <- unzip3 <$> mapM estimateInformations sublkls
    let (dvgmu,dvgvr) = estimateMeanVariance dvgs
        (mimu,mivr) = estimateMeanVariance mis
        (nrmdvgmu,nrmdvgvr) = estimateMeanVariance nrmdvgs
        stp = "Step " ++ (show . natValInt $ (Proxy :: Proxy k')) ++ "; dvgmu: " ++ show dvgmu
    return . trace stp $ Informations dvgmu (sqrt dvgvr) mimu (sqrt mivr) nrmdvgmu (sqrt nrmdvgvr)

analyzeDivergence0
    :: forall k s r . (Ord s, Read s, KnownNat k)
    => Int
    -> [([Int], s)]
    -> Proxy k
    -> Random r [Informations]
analyzeDivergence0 nsmps zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0

    (alldvgs0 :: B.Vector k Informations)
        <- B.generatePM' $ empiricalDivergenceStatistics zxss nsmps

    return $ B.toList alldvgs0

analyzeDivergence
    :: (Ord x, Read x)
    => Int
    -> [([Int],x)]
    -> Random r [Informations]
analyzeDivergence nsmps zxs =
    withNat (getPopulationSize zxs) $ analyzeDivergence0 nsmps zxs


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

    let clc = Collection clcstr

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    csvss <- case clcstr of
               "coen-cagli-2015" -> do
                   (zxss :: [[([Int],Int)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeDivergence nsmps) zxss
               "patterson-2013" -> do
                   (zxss :: [[([Int],Double)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeDivergence nsmps) zxss
               _ -> error "Invalid project"

    forM_ (zip csvss dsts) $ \(csvs, Dataset dststr') ->
        BS.writeFile ("projects/" ++ clcstr ++ "/analysis/dvg/" ++ dststr' ++ ".csv")
        $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (dvgOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



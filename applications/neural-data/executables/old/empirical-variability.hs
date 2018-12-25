{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import Goal.Plot
import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV
import qualified Data.Map as M

import Data.Semigroup ((<>))



--- Analysis ---

data Variation = Variation
    { meanResponseVariation :: Double
    , sdResponseVariation :: Double
    , meanTuningCurveVariation :: Double
    , sdTuningCurveVariation :: Double }
    deriving (Generic, Show)

instance CSV.FromNamedRecord Variation
instance CSV.ToNamedRecord Variation
instance CSV.DefaultOrdered Variation
instance NFData Variation


responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)


tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

empiricalVariabilityStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => [(Response k,s)]
    -> Int
    -> Proxy k'
    -> Random r Variation
empiricalVariabilityStatistics zxss n _ = do
    let zxmp = stimulusResponseMap zxss
        nzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let sidxss = G.convert <$> idxss
        subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleEmpiricalTuningCurves nzxmp <$> sidxss
        rvars = sqrt . snd . estimateMeanVariance . map fromIntegral . responseSums <$> subrs
        tcvars = sqrt . snd . estimateMeanVariance . tuningCurveSums <$> subtcss
        (rmu,rvr) = estimateMeanVariance rvars
        (tcmu,tvarr) = estimateMeanVariance tcvars
    return $ Variation rmu (sqrt rvr) tcmu (sqrt tvarr)

analyzeVariability0
    :: forall k s r . (Ord s, Read s, KnownNat k)
    => Int
    -> [([Int], s)]
    -> Proxy k
    -> Random r [Variation]
analyzeVariability0 nsmps zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0

    (allvars0 :: B.Vector k Variation)
        <- B.generatePM' $ empiricalVariabilityStatistics zxss nsmps

    let allvars = B.toList allvars0

    return allvars

analyzeVariability
    :: (Ord x, Read x)
    => Int
    -> [([Int],x)]
    -> Random r [Variation]
analyzeVariability nsmps zxs =
    withNat (getPopulationSize zxs) $ analyzeVariability0 nsmps zxs


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

varOpts :: Parser AnalysisOpts
varOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 10)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts clcstr dststr nsmps) = do

    let clc = Collection clcstr

    let pth = "projects/" ++ clcstr ++ "/analysis/var"

    createDirectoryIfMissing True pth

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    csvss <- case clcstr of
               "coen-cagli-2015" -> do
                   (zxss :: [[([Int],Int)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeVariability nsmps) zxss
               "patterson-2013" -> do
                   (zxss :: [[([Int],Double)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeVariability nsmps) zxss
               's':'y':'n':'t':'h':_ -> do
                   (zxss :: [[([Int],Double)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeVariability nsmps) zxss
               _ -> error "Invalid project"

    forM_ (zip csvss dsts) $ \(csvs, Dataset dststr') ->
        BS.writeFile (pth ++ "/" ++ dststr' ++ ".csv") $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (varOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



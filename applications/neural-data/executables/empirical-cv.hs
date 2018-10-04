{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

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


responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)

tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

empiricalCVStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => [(Response k,s)]
    -> Int
    -> Proxy k'
    -> Random r CoefficientsOfVariation
empiricalCVStatistics zxss n _ = do
    let zxmp = stimulusResponseMap zxss
        nzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let sidxss = G.convert <$> idxss
        subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleEmpiricalTuningCurves nzxmp <$> sidxss
        rcvs = sqrt . snd . estimateMeanVariance . map fromIntegral . responseSums <$> subrs
        tccvs = sqrt . snd . estimateMeanVariance . tuningCurveSums <$> subtcss
        (rmu,rvr) = estimateMeanVariance rcvs
        (tcmu,tcvr) = estimateMeanVariance tccvs
    return $ CoefficientsOfVariation rmu (sqrt rvr) tcmu (sqrt tcvr)

analyzeCoefficientOfVariation0
    :: forall k s r . (Ord s, Read s, KnownNat k)
    => Int
    -> [([Int], s)]
    -> Proxy k
    -> Random r [CoefficientsOfVariation]
analyzeCoefficientOfVariation0 nsmps zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0

    (allcvs0 :: B.Vector k CoefficientsOfVariation)
        <- B.generatePM' $ empiricalCVStatistics zxss nsmps

    let allcvs = B.toList allcvs0

    return allcvs

analyzeCoefficientOfVariation
    :: (Ord x, Read x)
    => Int
    -> [([Int],x)]
    -> Random r [CoefficientsOfVariation]
analyzeCoefficientOfVariation nsmps zxs =
    withNat (getPopulationSize zxs) $ analyzeCoefficientOfVariation0 nsmps zxs


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
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
                   realize $ mapM (analyzeCoefficientOfVariation nsmps) zxss
               "patterson-2013" -> do
                   (zxss :: [[([Int],Double)]]) <- mapM (getNeuralData clc) dsts
                   realize $ mapM (analyzeCoefficientOfVariation nsmps) zxss
               _ -> error "Invalid project"

    forM_ (zip csvss dsts) $ \(csvs, Dataset dststr') ->
        BS.writeFile ("projects/" ++ clcstr ++ "/analysis/cv/" ++ dststr' ++ ".csv")
        $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



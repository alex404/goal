{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import Data.Semigroup ((<>))


--- Globals ---


--- Analysis ---


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
        <- B.generatePM' $ fmap fst . responseStatistics zxss nsmps

    return $ B.toList allcvs0

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
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 1)

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

    forM_ csvss $ \csvs ->
        BS.writeFile (clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv")
        $ CSV.encodeDefaultOrderedByName csvs


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



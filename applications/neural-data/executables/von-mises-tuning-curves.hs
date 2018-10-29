{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Data.Semigroup ((<>))
import qualified Data.List as L

import qualified Goal.Core.Vector.Storable as S

--- Globals ---


ananm :: String
ananm = "tuning-curves"

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 1000

--- CLI ---


--- Functions ---


--- CLI ---

analyzeTuningCurves
    :: forall k . KnownNat k
    => [([Int],Double)]
    -> Proxy k
    -> [[Double]]
analyzeTuningCurves zxs0 _ =
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
        ppc = fitPPC zxs
        (rho0,rprms) = populationCodeRectificationParameters ppc xsmps
        rcrv = rectificationCurve rho0 rprms xsmps
        tcs = listCoordinates . dualTransition <$> ppc >$>* xsmps
        mxtcs = maximum <$> tcs
        (mupr,vrpr) = estimateMeanVariance . map (head . tail . listCoordinates . toSource)
            . S.toList . toRows . snd $ splitAffine ppc
        stdoustr = concat ["mupr: ", show mupr, "; sdpr: ", show $ sqrt vrpr]
    in trace stdoustr $ zipWith (++) (L.transpose [xsmps,sumOfTuningCurves ppc xsmps,rcrv,mxtcs]) tcs

data AnalysisOpts = AnalysisOpts String String Int

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto
        (long "tck" <> help "subsampled population size" <> short 'k' <> value 0)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg knrns0) = do

    dsts <- if dstarg == ""
               then fromJust <$> goalReadDatasetsCSV prjnm expnm
               else return [Dataset dstarg]

    forM_ dsts $ \dst -> do

        zxs <- getNeuralData expnm dst

        let knrns = if knrns0 == 0
                       then getPopulationSize zxs
                       else knrns0

        let tcs = withNat knrns $ analyzeTuningCurves zxs

        goalWriteAnalysis' prjnm expnm ananm (Just dst) tcs


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



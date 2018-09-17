{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import Paths_neural_data
import Data.Semigroup ((<>))


--- CSV ---


data Coefficients = Coefficients
    { responseCV :: Double
    , minimumCV :: Double
    , averageCV :: Double
    , maximumCV :: Double }
    deriving Generic

instance FromNamedRecord Coefficients
instance ToNamedRecord Coefficients
instance DefaultOrdered Coefficients
instance NFData Coefficients


--- Globals ---


--- Functions ---


tuningCurveSums
    :: (Ord s, KnownNat k)
    => M.Map s (Mean # Neurons k)
    -> M.Map s Double
tuningCurveSums mzxmp = S.sum . coordinates <$> mzxmp

responseSums
    :: M.Map s [Response k]
    -> [Int]
responseSums zs = sum <$> concat (M.elems zs)

responseStatistics
    :: forall s k k' r . (Ord s, KnownNat k, KnownNat k', k' <= k)
    => M.Map s [Response k]
    -> Int
    -> Proxy k'
    -> Random r Coefficients
responseStatistics zxmp n _ = do
    let mzxmp = empiricalTuningCurves zxmp
    (idxss :: [B.Vector k' Int]) <- replicateM n $ generateIndices (Proxy :: Proxy k)
    let subrs = subSampleResponses zxmp <$> idxss
        subtcss = subSampleTuningCurves mzxmp . G.convert <$> idxss
        rcvs = estimateCoefficientOfVariation . map fromIntegral . responseSums <$> subrs
        scvs = estimateCoefficientOfVariation . tuningCurveSums <$> subtcss
    return $ Coefficients (average rcvs) (minimum scvs) (average scvs) (maximum scvs)


--- CLI ---


data CVOpts = CVOpts Bool Int

cvOpts :: Parser CVOpts
cvOpts = CVOpts
    <$> switch (long "analyze" <> short 'a' <> help "analyze the data")
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 1)

data AllOpts = AllOpts CVOpts GNUPlotOpts

allOpts :: Parser AllOpts
allOpts = AllOpts <$> cvOpts <*> gnuPlotOpts

--- Main ---


analyzeCoefficientOfVariation
    :: forall x . (Ord x, Read x) => String -> Proxy x -> Int -> Dataset -> IO ()
analyzeCoefficientOfVariation prj _ nsmps dst = do
    (zxs :: [([Int],x)]) <- getNeuralData prj dst
    let k = getPopulationSize zxs
    withNat k $ analyzeCoefficientOfVariation0 nsmps prj dst zxs


analyzeCoefficientOfVariation0
    :: forall k s . (Ord s, Read s, KnownNat k)
    => Int
    -> String
    -> Dataset
    -> [([Int], s)]
    -> Proxy k
    -> IO ()
analyzeCoefficientOfVariation0 nsmps prj nd zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0
        zxmp = stimulusResponseMap zxss

    (allcvs :: B.Vector k Coefficients) <- realize (B.generatePM' $ responseStatistics zxmp nsmps)

    let prj' =  prj ++ "/analysis/cv"

    goalWriteCSV prj' (datasetName nd) $ B.toList allcvs

runOpts :: AllOpts -> IO ()
runOpts (AllOpts (CVOpts abl nsmps) (GNUPlotOpts prj dststr pbl ibl)) = do

    void . goalCreateProject $ prj ++ "/plots/cv"

    dss <- if dststr == ""
              then goalReadCSV prj "datasets"
              else return [Dataset dststr]

    when abl $ do
        let mapper =
                case prj of
                  "patterson-2013" -> analyzeCoefficientOfVariation prj (Proxy :: Proxy Double) nsmps
                  "coen-cagli-2015" -> analyzeCoefficientOfVariation prj (Proxy :: Proxy Int) nsmps
                  _ -> error "Invalid Project"

        mapM_ mapper dss

    let gplts = (\(Dataset ds) -> GNUPlotOpts prj ds pbl ibl) <$> dss
    gpltpth <- getDataFileName "plots/neural-data-cv.gpi"
    mapM_ (runGNUPlotOpts gpltpth) gplts


main :: IO ()
main = runOpts =<< execParser opts
  where
    opts = info (allOpts <**> helper)
      ( fullDesc
     <> progDesc "Analyze and Plot Neural Data"
     <> header "Analyze and Plot Neural Data" )

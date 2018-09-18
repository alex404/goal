{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

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


data CVOpts = CVOpts String String Int

cvOpts :: Parser CVOpts
cvOpts = CVOpts
    <$> strOption
        ( long "collection" <> short 'c' <> help "Which data collection to plot" <> value "")
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 1)

--- Main ---


analyzeCoefficientOfVariation
    :: forall x . (Ord x, Read x) => Proxy x -> Int -> Collection -> Dataset -> IO ()
analyzeCoefficientOfVariation _ nsmps clc dst = do
    (zxs :: [([Int],x)]) <- getNeuralData clc dst
    let k = getPopulationSize zxs
    withNat k $ analyzeCoefficientOfVariation0 zxs nsmps clc dst


analyzeCoefficientOfVariation0
    :: forall k s . (Ord s, Read s, KnownNat k)
    => [([Int], s)]
    -> Int
    -> Collection
    -> Dataset
    -> Proxy k
    -> IO ()
analyzeCoefficientOfVariation0 zxss0 nsmps (Collection clcstr) (Dataset dststr) _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0
        zxmp = stimulusResponseMap zxss

    (allcvs :: B.Vector k Coefficients) <- realize (B.generatePM' $ responseStatistics zxmp nsmps)

    BS.writeFile (clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv")
        . CSV.encodeDefaultOrderedByName $ B.toList allcvs

runOpts :: CVOpts -> IO ()
runOpts (CVOpts clcstr dststr nsmps) = do

    void . goalCreateProject $ clcstr ++ "/plots/cv"

    clcs <- if clcstr == ""
                 then getCollections
                 else return [Collection clcstr]
    dstss <- if dststr == ""
               then mapM getDatasets clcs
               else return [[Dataset dststr]]

    sequence_ $ do
        (clc,dsts) <- zip clcs dstss
        dst <- dsts
        return $ case clc of
          Collection "patterson-2013" ->
              analyzeCoefficientOfVariation (Proxy :: Proxy Double) nsmps clc dst
          Collection "coen-cagli-2015" ->
              analyzeCoefficientOfVariation (Proxy :: Proxy Int) nsmps clc dst
          _ -> error "Invalid Project"

main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts


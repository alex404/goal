{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

import qualified Data.Map as M
import System.Process
import Paths_neural_data


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

nsmps :: Int
nsmps = 100

ptprj :: String
ptprj = "patterson-2013"

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


--- Main ---


analyzeCoefficientOfVariation
    :: forall x . (Ord x, Read x) => Proxy x -> Dataset -> IO ()
analyzeCoefficientOfVariation _ dst = do
    (zxs :: [([Int],x)]) <- getNeuralData ptprj dst
    let k = getPopulationSize zxs
    withNat k $ analyzeCoefficientOfVariation0 ptprj dst zxs


analyzeCoefficientOfVariation0
    :: forall k s . (Ord s, Read s, KnownNat k)
    => String
    -> Dataset
    -> [([Int], s)]
    -> Proxy k
    -> IO ()
analyzeCoefficientOfVariation0 prj nd zxss0 _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0
        zxmp = stimulusResponseMap zxss

    (allcvs :: B.Vector k Coefficients) <- realize (B.generatePM' $ responseStatistics zxmp nsmps)

    let prj' =  prj ++ "/analysis/cv"

    goalWriteCSV prj' (datasetName nd) $ B.toList allcvs

main :: IO ()
main = do

    ptdsts <- goalReadCSV ptprj "datasets"
    mapM_ (analyzeCoefficientOfVariation (Proxy :: Proxy Double)) ptdsts
    flnm <- getDataFileName "plots/plot-neural-cv.gpi"
    void $ spawnProcess "gnuplot" [flnm]

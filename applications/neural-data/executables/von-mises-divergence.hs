{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Generic as G

import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import Data.Semigroup ((<>))


--- CSV ---


data PosteriorDivergence = PosteriorDivergence
    { meanDivergence :: Double
    , sdDivergence :: Double }
    deriving (Generic, Show)

instance FromNamedRecord PosteriorDivergence
instance ToNamedRecord PosteriorDivergence
instance DefaultOrdered PosteriorDivergence
instance NFData PosteriorDivergence

errbnd :: Double
errbnd = 1e-12

stmsmps :: [Double]
stmsmps = tail $ range 0 (2*pi) 1000

--- Analysis ---

approximateConditionalLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
approximateConditionalLogPartitionFunction lkl z =
    potential $ z *<.< snd (splitAffine lkl)

conditionalLogPartitionFunction
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Double
conditionalLogPartitionFunction lkl z =
    let nx = z *<.< snd (splitAffine lkl)
        logupst x = nx <.> sufficientStatistic x - head (sumOfTuningCurves lkl [x])
        mx = maximum $ logupst <$> stmsmps
        upst0 x = exp $ logupst x - mx
     in (+ mx) . log1p . subtract 1 . fst $ integrate errbnd upst0 0 (2*pi)

-- Assumes a uniform prior over stimuli
estimateDivergence
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Random r Double
estimateDivergence lkl = do
    let nzs = dualTransition <$> lkl >$>* stmsmps
        stcavg = average $ sum . listCoordinates <$> nzs
    zs <- mapM samplePoint nzs
    let dffs = average [ approximateConditionalLogPartitionFunction lkl z
               - conditionalLogPartitionFunction lkl z | z <- zs ]
    return $ dffs - stcavg

vonMisesDivergenceStatistics
    :: forall k k' r . (KnownNat k, KnownNat k', k' <= k)
    => [(Response k,Double)]
    -> Int
    -> Proxy k'
    -> Random r PosteriorDivergence
vonMisesDivergenceStatistics zxss nsmps _ = do
    let lkl = fitPPC 200 (-0.2) zxss
    (idxss :: [B.Vector k' Int]) <- replicateM nsmps $ generateIndices (Proxy :: Proxy k)
    let sidxss = G.convert <$> idxss
        subtcss = subSamplePPC lkl <$> sidxss
    dvgs <- mapM estimateDivergence subtcss
    let (dvgmu,dvgvr) = estimateMeanVariance dvgs
    return $ PosteriorDivergence dvgmu (sqrt dvgvr)

analyzeDivergence0
    :: forall k r . KnownNat k
    => Int
    -> [([Int], Double)]
    -> Proxy k
    -> Random r [PosteriorDivergence]
analyzeDivergence0 nsmps zxss0 _ = do

    let zxss :: [(Response k, Double)]
        zxss = strengthenNeuralData zxss0

    (alldvgs0 :: B.Vector k PosteriorDivergence)
        <- B.generatePM' $ vonMisesDivergenceStatistics zxss nsmps

    return $ B.toList alldvgs0

analyzeDivergence
    :: Int
    -> [([Int],Double)]
    -> Random r [PosteriorDivergence]
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



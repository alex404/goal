{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
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

import Data.Semigroup ((<>))
import Data.List


--- Globals ---


stms :: [Double]
stms = init $ range 0 (2*pi) 9

xsmps :: [Double]
xsmps = init $ range 0 (2*pi) 1000


--- CSV ---


data Posteriors = Posteriors
    { stimulusValue :: Double
    , nNeurons :: Int
    , numericalPosterior :: Double
    , approximatePosterior :: Double }
    deriving (Generic, Show)

instance CSV.FromNamedRecord Posteriors
instance CSV.ToNamedRecord Posteriors
instance CSV.DefaultOrdered Posteriors
instance NFData Posteriors


--- Functions ---


--- CLI ---

-- Under the assumption of a flat prior
posteriorPoints0
    :: forall k k' . (KnownNat k, KnownNat k', k' <= k)
    => Mean #> Natural # Neurons k <* VonMises
    -> Response k
    -> Proxy k'
    -> [Posteriors]
posteriorPoints0 ppc z prxk' =
    let k' = natValInt prxk'
        ppc' :: Mean #> Natural # Neurons k' <* VonMises
        ppc' = subPPC ppc
        (z',_ :: B.Vector (k - k') Int) = B.splitAt z
        npstrs = numericalVonMisesPPCPosterior ppc' z' <$> xsmps
        apstrs = approximateVonMisesPPCPosterior ppc' z' <$> xsmps
     in zipWith4 Posteriors xsmps (repeat k') npstrs apstrs

-- Under the assumption of a flat prior
posteriorPoints
    :: forall k r . KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Double
    -> Random r [Posteriors]
posteriorPoints ppc x = do
    idxs <- generateIndices (Proxy :: Proxy k)
    let ppc' :: Mean #> Natural # Neurons k <* VonMises
        ppc' = subSamplePPC ppc $ G.convert idxs
    z <- samplePoint $ ppc' >.>* x
    let pstss :: B.Vector k [Posteriors]
        pstss = B.generateP (posteriorPoints0 ppc' z)
    return . concat $ B.toList pstss

examplePosteriorsFromData
    :: forall k r . KnownNat k
    => [([Int],Double)]
    -> Double
    -> Proxy k
    -> Random r [Posteriors]
examplePosteriorsFromData zxs0 x _ = do
    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0
        ppc = fitPPC zxs
    posteriorPoints ppc x

data AnalysisOpts = AnalysisOpts String String Double

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto
        (long "stimulus" <> help "stimulus value" <> short 's' <> value pi)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts clcstr dststr stm) = do

    let clc = Collection clcstr

    let pth = "projects/" ++ clcstr ++ "/analysis/psts"

    createDirectoryIfMissing True pth

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    forM_ dsts $ \dst@(Dataset dstflnm) -> do

        zxs <- getNeuralData clc dst

        let knrns = getPopulationSize zxs

        psts <- realize (withNat knrns $ examplePosteriorsFromData zxs stm)
        let flnm = concat [pth,"/",dstflnm,".csv"]

        BS.writeFile flnm $ CSV.encodeDefaultOrderedByName psts


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



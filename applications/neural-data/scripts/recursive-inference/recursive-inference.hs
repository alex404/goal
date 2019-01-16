{-# LANGUAGE
    DeriveGeneric,
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    ScopedTypeVariables,
    DataKinds #-}

import NeuralData
import NeuralData.VonMises

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Data.Semigroup ((<>))


--- CSV ---


data RecursiveInformations = RecursiveInformations
    { meanLinearDivergence :: Double
    , meanAffineDivergence :: Double
    , meanDecoderDivergence :: Double
    , meanVonMisesMutualInformation :: Double
    , meanLinearDivergenceRatio :: Double
    , meanAffineDivergenceRatio :: Double
    , meanDecoderDivergenceRatio :: Double }
    deriving (Generic, Show)

instance FromNamedRecord RecursiveInformations
instance ToNamedRecord RecursiveInformations
instance DefaultOrdered RecursiveInformations
instance NFData RecursiveInformations



--- Analysis ---

ananm :: String
ananm = "recursive-inference"

nstms :: Int
nstms = 1000

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) (nstms+1)

informationsFolder
    :: KnownNat k
    => Int
    -> Mean #> Natural # Neurons k <* VonMises
    -> Natural # VonMises
    -> Mean #> Natural # VonMises <* Neurons k
    -> (Double,Double,Double,Double,Double)
    -> Double
    -> Random r (Double,Double,Double,Double,Double)
informationsFolder n lkl rprms dcd (zpn,zqn,zcn,ptn,dcddvg) x = do
    zs <- sample n $ lkl >.>* x
    let (dns,zpn') = numericalRecursiveBayesianInference 1e-12 0 (2*pi) 100 (repeat lkl) zs (const 1)
        zqn' = potential $ rectifiedRecursiveBayesianInference' zero lkl zs zero
        zcn' = potential $ rectifiedRecursiveBayesianInference' rprms lkl zs zero
        ptn' = sum [ sufficientStatistic z <.> (snd (splitAffine lkl) >.>* x) | z <- zs ]
        dcddvg' = linearDecoderDivergence dcd dns $ head zs
    return (zpn + zpn',zqn + zqn',zcn + zcn',ptn + ptn',dcddvg + dcddvg')

-- Assumes a uniform prior over stimuli
estimateInformations
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Mean #> Natural # VonMises <* Neurons k
    -> Int
    -> Random r RecursiveInformations
estimateInformations lkl dcd n = do
    let (rho0,rprms) = regressRectificationParameters lkl xsmps
    (zpnavg0,zqnavg0,zcnavg0,ptnavg0,dcddvg0) <- foldM (informationsFolder n lkl rprms dcd) (0,0,0,0,0) xsmps
    let k' = fromIntegral nstms
        n' = fromIntegral n
        (zpnavg,zqnavg,zcnavg,ptnavg,dcddvg) = (zpnavg0/k',zqnavg0/k',zcnavg0/k',ptnavg0/k',dcddvg0/k')
        pq0dvg = zqnavg - zpnavg - rho0*n'
        pqdvg = zcnavg - zpnavg - rho0*n'
        mi = ptnavg - zpnavg - rho0*n'
    return . traceGiven $ RecursiveInformations pq0dvg pqdvg dcddvg mi (pq0dvg/mi) (pqdvg/mi) (dcddvg/mi)

vonMisesInformationsStatistics
    :: KnownNat k
    => Int
    -> Mean #> Natural # Neurons k <* VonMises
    -> Random r [RecursiveInformations]
vonMisesInformationsStatistics n lkl = do
    dcd <- fitLinearDecoder lkl xsmps
    mapM (estimateInformations lkl dcd) [1..n]

fitAnalyzeInformations
    :: forall k r . KnownNat k
    => Int
    -> [([Int], Double)]
    -> Proxy k
    -> Random r [RecursiveInformations]
fitAnalyzeInformations n zxss0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxss0

    lkl <- fitIPLikelihood zxs

    vonMisesInformationsStatistics n lkl


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int

vminfOpts :: Parser AnalysisOpts
vminfOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "nsteps" <> help "number of inference steps" <> short 'n' <> value 10)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg nsmps) = do

    let expmnt = Experiment prjnm expnm

    dsts <- if null dstarg
               then fromJust <$> goalReadDatasetsCSV expmnt
               else return [dstarg]

    forM_ dsts $ \dst -> do

        (k,zxs :: [([Int], Double)]) <- getNeuralData expnm dst

        let rinfs = case someNatVal k of
                    SomeNat prxk -> fitAnalyzeInformations nsmps zxs prxk

        infs <- realize rinfs

        goalWriteNamedAnalysis True expmnt (Just $ SubExperiment ananm dst) infs


--- Main ---


main :: IO ()
main = do

    let opts = info (vminfOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



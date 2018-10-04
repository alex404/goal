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

import Data.Semigroup ((<>))


--- Globals ---


eps :: Double
eps = -0.2

nepchs :: Int
nepchs = 100

xsmps :: [Double]
xsmps = tail $ range 0 (2*pi) 9

--- CLI ---


--- Functions ---


ppcStimulusDerivatives
    :: KnownNat k => Mean #> Natural # Neurons k <* VonMises
    -> SamplePoint VonMises
    -> S.Vector k Double
ppcStimulusDerivatives ppc x =
    let fxs = coordinates . dualTransition $ ppc >.> mx
        tcs = toRows . snd $ splitAffine ppc
     in S.zipWith zipper fxs tcs
    where mx = sufficientStatistic x
          (cx,sx) = S.toPair $ coordinates mx
          zipper fx (Point cs) =
              let (tht1,tht2) = S.toPair cs
               in fx*(cx * tht2 - sx * tht1)

fisherInformation
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> SamplePoint VonMises
    -> Double
fisherInformation ppc x =
    let fxs2' = S.map square $ ppcStimulusDerivatives ppc x
        fxs = coordinates . dualTransition $ ppc >.>* x
     in S.sum $ S.zipWith (/) fxs2' fxs

averageFisherInformation
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Double
averageFisherInformation ppc =
    average $ fisherInformation ppc <$> tail (range 0 (2*pi) 100)

inferenceStatistics
    :: forall k k' r . (KnownNat k, KnownNat k', k' <= k)
    => Int
    -> [(Response k,Double)]
    -> Proxy k'
    -> Random r ([[Double]],[[Double]],[[Posteriors]])
inferenceStatistics n zxss prxk' = do
    let ppc = fitPPC zxss
        err = 1e-12
    (idxss,mnidxs,mxidxs) <- snd <$> empiricalCVStatistics zxss n prxk'
    let mnppc = subSamplePPC ppc mnidxs
        mxppc = subSamplePPC ppc mxidxs
    let ppcss = subSamplePPC ppc . G.convert <$> idxss
    pstrss <- examplePosteriors err xsmps 200 $ head ppcss
    let mntcs = zipWith (:) xsmps $ listCoordinates . dualTransition <$> mnppc >$>* xsmps
        mxtcs = zipWith (:) xsmps $ listCoordinates . dualTransition <$> mxppc >$>* xsmps
    return (mntcs,mxtcs,pstrss)


--- Analysis ---


inferenceAnalysis0
    :: forall k r . KnownNat k
    => Int
    -> Int
    -> [([Int], Double)]
    -> Proxy k
    -> Random r ([[Double]],[[Double]],[[Posteriors]])
inferenceAnalysis0 nsmps nsbnrns zxs0 _ = do

    let zxs :: [(Response k, Double)]
        zxs = strengthenNeuralData zxs0

    (allcvs :: B.Vector k ([[Double]],[[Double]],[[Posteriors]]))
        <- B.generatePM' (inferenceStatistics nsmps zxs)

    return $  B.unsafeIndex allcvs nsbnrns

inferenceAnalysis
    :: Int
    -> Int
    -> [([Int],Double)]
    -> Random r ([[Double]],[[Double]],[[Posteriors]])
inferenceAnalysis nsmps nsbnrns zxs =
    withNat (getPopulationSize zxs) $ inferenceAnalysis0 nsmps nsbnrns zxs


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Int Int

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument
        ( help "Which data collection to plot" )
    <*> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto
        (long "nsamples" <> help "number of samples to generate" <> short 'n' <> value 10)
    <*> option auto
        (long "tck" <> help "subsampled population size" <> short 'k' <> value 10)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts clcstr dststr nsmps nnrns) = do

    let clc = Collection clcstr

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    zxss <- mapM (getNeuralData clc) dsts

    allcsvss <- realize $ mapM (inferenceAnalysis nsmps nnrns) zxss

    forM_ (zip dsts allcsvss) $ \(Dataset dststr', (mntcs,mxtcs,pstrss)) -> do

        BS.writeFile ("projects/" ++ clcstr ++ "/analysis/tcs/" ++ dststr' ++ "-min.csv")
            $ CSV.encode mntcs

        BS.writeFile ("projects/" ++ clcstr ++ "/analysis/tcs/" ++ dststr' ++ "-max.csv")
            $ CSV.encode mxtcs

        sequence_ $ do
            (k,pstrs) <- zip [0 :: Int ..] pstrss
            let flnm = "projects/" ++ clcstr ++ "/analysis/inf/" ++ dststr' ++ "-stim-" ++ show k ++ ".csv"
            return . BS.writeFile flnm $ CSV.encodeDefaultOrderedByName pstrs


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts



{-# LANGUAGE DeriveGeneric,TupleSections,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B

import qualified Data.Map as M
import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import GHC.Generics


--- CSV ---

data VonMisesCSV = VonMisesCSV
    { averageFisherInformation :: Double
    , minimalCVFisherInformation :: Double
    , maximalCVFisherInformation :: Double }
    deriving Generic

instance FromNamedRecord VonMisesCSV
instance ToNamedRecord VonMisesCSV
instance DefaultOrdered VonMisesCSV
instance NFData VonMisesCSV


--- Globals ---


eps :: Double
eps = -0.2

nepchs :: Int
nepchs = 100

pltsmps :: [Double]
pltsmps = range 0 (2*pi) 500

--- Functions ---


-- Von Mises Statistics --

fitPPC
    :: forall k . KnownNat k
    => M.Map Double [Response k]
    -> ([Double],[Response k],[Mean #> Natural # Neurons k <* VonMises])
fitPPC xzmp =
    let sps :: S.Vector k (Source # VonMises)
        sps = S.map (\mu -> Point $ S.fromTuple (mu,1)) $ S.range 0 (2*pi)
        ppc0 = vonMisesPopulationEncoder True (Left 1) sps
        (xs,zs) = unzip $ concat [ (x,) <$> zs0 | (x,zs0) <- M.toList xzmp ]
        backprop p = joinTangentPair p $ stochasticConditionalCrossEntropyDifferential xs zs p
     in (xs,zs,take nepchs $ vanillaGradientSequence backprop eps defaultAdamPursuit ppc0)

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

-- Under the assumption of a flat prior
--numericalVonMisesPPCPosterior0
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Response k
--    -> Double
--    -> Double
--numericalVonMisesPPCPosterior0 lkl z =
--    let tcs = snd $ splitAffine lkl
--        nxs = z *<.< tcs
--        uldns x = nxs <.> sufficientStatistic x - S.sum (coordinates . dualTransition $ lkl >.>* x)
--        avg = integrate 1e-500 uldns 0 (2*pi) / (2*pi)
--        udns x = exp $ uldns x - avg
--        nrm = trace (show $ udns <$> range 0 (2*pi) 8) $ integrate 1e-5000 udns 0 (2*pi)
--     in (/nrm) . udns
--
---- Under the assumption of a flat prior
--approximateVonMisesPPCPosterior0
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Response k
--    -> Double
--    -> Double
--approximateVonMisesPPCPosterior0 lkl z =
--    let tcs = snd $ splitAffine lkl
--     in density (z *<.< tcs)


--- CLI ---


newtype VMOpts = VMOpts String

vmOpts :: Parser VMOpts
vmOpts = VMOpts
    <$> strOption
        ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")


runOpts :: VMOpts -> IO ()
runOpts (VMOpts dststr) = do

    dsts <- if dststr == ""
               then getDatasets $ Collection "patterson-2013"
               else return [Dataset dststr]

    mapM_ analyzeVonMises dsts


--- Main ---


analyzeVonMisesStatistics
    :: forall x . (Ord x, Read x) => Proxy x -> Int -> Collection -> Dataset -> IO ()
analyzeVonMisesStatistics _ nsmps clc dst = do
    (zxs :: [([Int],x)]) <- getNeuralData clc dst
    let k = getPopulationSize zxs
    withNat k $ analyzeVonMisesStatistics0 zxs nsmps clc dst

analyzeVonMisesStatistics0
    :: forall k s . (Ord s, Read s, KnownNat k)
    => [([Int], s)]
    -> Int
    -> Collection
    -> Dataset
    -> Proxy k
    -> IO ()
analyzeVonMisesStatistics0 zxss0 nsmps (Collection clcstr) (Dataset dststr) _ = do

    let zxss :: [(Response k, s)]
        zxss = strengthenNeuralData zxss0
        zxmp = stimulusResponseMap zxss

    (allcvs :: B.Vector k VonMisesCSV) <- realize (B.generatePM' $ vonMisesStatistics zxmp nsmps)

    BS.writeFile (clcstr ++ "/analysis/cv/" ++ dststr ++ ".csv")
        . CSV.encodeDefaultOrderedByName $ B.toList allcvs

main :: IO ()
main = do

    let opts = info (vmOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the fisher information of neural population data"

    runOpts =<< execParser opts


--- Plot Graveyard ---

--tuningCurveLayout
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Layout Double Double
--tuningCurveLayout lkl = execEC $ do
--
--    goalLayout
--    radiansAbscissa
--
--    layout_x_axis . laxis_title .= "Stimulus"
--    layout_y_axis . laxis_title .= "Rate"
--    --layout_x_axis . laxis_override .= axisGridHide . axisLabelsOverride [(0,"0"),(0.5,"0.5"),(1,"1"),(1.5,"1.5")]
--
--    plot . liftEC $ do
--        plot_lines_style .= solidLine 3 (opaque red)
--        plot_lines_values .= toList (toList <$> tuningCurves pltsmps lkl)
--
--    plot . liftEC $ do
--        plot_lines_style .= solidLine 3 (opaque black)
--        plot_lines_values .= [ zip pltsmps $ sumOfTuningCurves lkl pltsmps ]
--
--
--decodingLayout
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> Response k
--    -> Layout Double Double
--decodingLayout lkl z = execEC $ do
--
--    goalLayout
--    radiansAbscissa
--
--    layout_x_axis . laxis_title .= "Stimulus"
--    layout_y_axis . laxis_title .= "Posterior Density"
--
--    let pstnm = numericalVonMisesPPCPosterior0 lkl z
--        pstht = approximateVonMisesPPCPosterior0 lkl z
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque red)
--        plot_lines_values .= [zip pltsmps $ pstnm <$> pltsmps]
--        plot_lines_title .= "Numerical Posterior"
--
--    plot . liftEC $ do
--
--        plot_lines_style .= solidLine 3 (opaque blue)
--        plot_lines_values .= [zip pltsmps $ pstht <$> pltsmps]
--        plot_lines_title .= "Approximate Posterior"
--
--generateDecodingLayouts
--    :: KnownNat k
--    => Mean #> Natural # Neurons k <* VonMises
--    -> [Double]
--    -> Random r [Layout Double Double]
--generateDecodingLayouts lkl xs = do
--    zs <- mapM samplePoint $ lkl >$>* xs
--    return $ decodingLayout lkl <$> zs
--
--

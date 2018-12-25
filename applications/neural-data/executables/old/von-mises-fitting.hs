{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}
{-# LANGUAGE DeriveGeneric,FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}


--- Imports ---


import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Boxed as B
import Data.Semigroup ((<>))
import qualified Data.List as L


--- Globals ---


ananm :: String
ananm = "von-mises-fitting"

nwghts0 :: (KnownNat n, 1 <= n) => Natural # Categorical Int n
nwghts0 = zero

swghts0 :: (KnownNat n, 1 <= n) => Source # Categorical Int n
swghts0 = transition nwghts0


--- CSV ---

data TrainingStatistics = TrainingStatistics
    { meanTrainingCrossEntropy :: Double
    , sdTrainingCrossEntropy :: Double
    , meanValidationCrossEntropy :: Double
    , sdValidationCrossEntropy :: Double }
    deriving (Generic, Show)

instance FromNamedRecord TrainingStatistics
instance ToNamedRecord TrainingStatistics
instance DefaultOrdered TrainingStatistics
instance NFData TrainingStatistics


--- Functions ---


ppcTraining
    :: forall k r . KnownNat k
    => Double
    -> Int
    -> [(Response k, Double)]
    -> [(Response k, Double)]
    -> Random r [(Double,Double)]
ppcTraining eps nepchs tzxs vzxs = do
    kps <- S.replicateM $ uniformR (0.2,0.6)
    let sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps $ S.range 0 (2*pi)
    gns <- S.replicateM $ uniformR (0.5,2)
    let ppc0 :: Mean #> Natural # Neurons k <* VonMises
        ppc0 = vonMisesPopulationEncoder True (Right $ Point gns) sps
        (tzs,txs) = unzip tzxs
        (vzs,vxs) = unzip vzxs
        backprop ppc = joinTangentPair ppc
            $ stochasticConditionalCrossEntropyDifferential txs tzs ppc
        ppcs = take nepchs $ vanillaGradientSequence backprop (negate eps) defaultAdamPursuit ppc0
        ents ppc = ( stochasticConditionalCrossEntropy txs tzs ppc
                   , stochasticConditionalCrossEntropy vxs vzs ppc )
        tcs = concat . map (('\n':) . show . listCoordinates . toSource) . S.toList . toRows . snd . splitAffine $ last ppcs
        trc2 = "\nPPC TuningCurves:\n" ++ tcs
    return . trace trc2 $ ents <$> ppcs

ppcTrainingStatistics
    :: forall k r . KnownNat k
    => Double
    -> Int
    -> Int
    -> Double
    -> [([Int], Double)]
    -> Proxy k
    -> Random r [TrainingStatistics]
ppcTrainingStatistics eps nepchs nsmps vld zxs0 _ = do
    let zxs1 :: [(Response k, Double)]
        zxs1 = strengthenNeuralData zxs0
    zxs <- shuffleList zxs1
    let spltn = round . (*vld) . fromIntegral $ length zxs
        (vzxs,tzxs) = splitAt spltn zxs
    (scess :: [[(Double,Double)]])
        <- replicateM nsmps $ ppcTraining eps nepchs tzxs vzxs
    return $ entropiesToTrainingStatistics <$> L.transpose scess

mixturePPCTraining
    :: forall n k r . (KnownNat k, KnownNat n, 1 <= n)
    => Double
    -> Int
    -> [(Response k, Double)]
    -> [(Response k, Double)]
    -> Proxy n
    -> Random r [(Double,Double)]
mixturePPCTraining eps nepchs tzxs vzxs _ = do
    kps <- S.replicateM $ uniformR (0.2,0.6)
    let sps = S.zipWith (\kp mu -> Point $ S.doubleton mu kp) kps $ S.range 0 (2*pi)
    gnss <- S.replicateM . fmap Point . S.replicateM $ uniformR (0.5,2)
    let mglm0 :: Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
        mglm0 = vonMisesMixturePopulationEncoder True swghts0 gnss sps
        (tzs,txs) = unzip tzxs
        (vzs,vxs) = unzip vzxs
        backprop mglm = joinTangentPair mglm
            $ mixtureStochasticConditionalCrossEntropyDifferential txs tzs mglm
        mglms = take nepchs $ vanillaGradientSequence backprop (negate eps) defaultAdamPursuit mglm0
        tracer mglm = concat
            [ "Average Absolute mGLM Coordinate: "
            , show . average $ abs <$> listCoordinates mglm
            , "\nInfinite mGLM Coordinate? "
            , show $ any (\c -> isNaN c || isInfinite c) $ listCoordinates mglm ]
        --ents mglm = trace (tracer mglm)
        ents mglm = ( mixtureStochasticConditionalCrossEntropy txs tzs mglm
                    , mixtureStochasticConditionalCrossEntropy vxs vzs mglm )
        tcs = concat . map (('\n':) . show . listCoordinates . toSource) . S.toList . toRows . fst . splitBottomSubLinear $ last mglms
        trc2 = "\nMGLM TuningCurves:\n" ++ tcs
    return . trace trc2 $ ents <$> mglms

entropiesToTrainingStatistics :: [(Double,Double)] -> TrainingStatistics
entropiesToTrainingStatistics tvces =
    let (tces,vces) = unzip tvces
        (tmu,tvr) = estimateMeanVariance tces
        (vmu,vvr) = estimateMeanVariance vces
     in TrainingStatistics tmu (sqrt tvr) vmu (sqrt vvr)

mixturePPCTrainingStatistics
    :: forall n k r . (KnownNat k, KnownNat n, 1 <= n)
    => Double
    -> Int
    -> Int
    -> Double
    -> [([Int], Double)]
    -> Proxy k
    -> Proxy n
    -> Random r [[TrainingStatistics]]
mixturePPCTrainingStatistics eps nepchs nsmps vld zxs0 _ _ = do
    let zxs1 :: [(Response k, Double)]
        zxs1 = strengthenNeuralData zxs0
    zxs <- shuffleList zxs1
    let spltn = round . (*vld) . fromIntegral $ length zxs
        (vzxs,tzxs) = splitAt spltn zxs
    (scesss :: B.Vector n [[(Double,Double)]])
        <- B.generatePM' (replicateM nsmps . mixturePPCTraining eps nepchs tzxs vzxs)
    return [ entropiesToTrainingStatistics <$> L.transpose scess | scess <- B.toList scesss ]


--- CLI ---


data AnalysisOpts = AnalysisOpts String String Double Int Int Int Double

cvOpts :: Parser AnalysisOpts
cvOpts = AnalysisOpts
    <$> strArgument ( help "Which data collection to plot" )
    <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
    <*> option auto (long "epsilon" <> short 'e' <> value 0.05)
    <*> option auto (long "nepochs" <> short 'n' <> value 1000)
    <*> option auto (long "nmixers" <> help "number of mixers; 0 means 'classic' algorithm only" <> short 'm' <> value 5)
    <*> option auto (long "nsmps" <> short 's' <> value 10)
    <*> option auto (long "validation" <> short 'v' <> help "fraction of the dataset to use for validation" <> value 0.2)

runOpts :: AnalysisOpts -> IO ()
runOpts (AnalysisOpts expnm dstarg eps nepchs m nsmps vld) = do

    dsts <- if dstarg == ""
               then fromJust <$> goalReadDatasetsCSV prjnm expnm
               else return [Dataset dstarg]

    if take 4 expnm == "true"

       then error "Fit-testing requires data"

       else forM_ dsts $ \dst -> do

                (k,(zxs :: [([Int], Double)])) <- getNeuralData expnm dst

                ppcsts <- realize $ withNat k (ppcTrainingStatistics eps nepchs nsmps vld zxs)
                goalWriteNamedAnalysis prjnm expnm ananm (Just dst) ppcsts

                when (m > 0) $ do
                    (mglmsts :: [[TrainingStatistics]]) <- realize
                        (withNat1 m (withNat k (mixturePPCTrainingStatistics eps nepchs nsmps vld zxs)))
                    mapM_ (goalAppendNamedAnalysis prjnm expnm ananm (Just dst)) mglmsts


--- Main ---


main :: IO ()
main = do

    let opts = info (cvOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Analyze the coefficient of variation of the total population activity in neural data"

    runOpts =<< execParser opts


--- Graveyard ---


--analyzeTuningCurves
--    :: forall k . KnownNat k
--    => [([Int],Double)]
--    -> Proxy k
--    -> [[Double]]
--analyzeTuningCurves zxs0 _ =
--    let zxs :: [(Response k, Double)]
--        zxs = strengthenNeuralData zxs0
--        ppc = fitPPC zxs
--        (rho0,rprms) = regressRectificationParameters ppc xsmps
--        rcrv = rectificationCurve rho0 rprms xsmps
--        tcs = listCoordinates . dualTransition <$> ppc >$>* xsmps
--        mxtcs = maximum <$> tcs
----        (mupr,vrpr) = estimateMeanVariance . map (head . tail . listCoordinates . toSource)
----            . S.toList . toRows . snd $ splitAffine ppc
--        stdoustr = concat ["mupr: ", show mupr, "; sdpr: ", show $ sqrt vrpr]
--    in trace stdoustr $ zipWith (++) (L.transpose [xsmps,potential <$> ppc >$>* xsmps,rcrv,mxtcs]) tcs



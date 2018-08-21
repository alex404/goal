{-# LANGUAGE FlexibleContexts,TypeFamilies,TypeOperators,ScopedTypeVariables,DataKinds #-}

import KohnLab

import Goal.Core
import Goal.Geometry
import Goal.Probability
import qualified Goal.Core.Vector.Boxed as B
import qualified Goal.Core.Vector.Storable as S
import qualified Data.Vector.Storable as V

import qualified Data.Map as M
import qualified System.Random.MWC.Distributions as MWC
import qualified Data.List as L
import qualified Numeric.AD as AD


--- Globals ---


--- Functions ---


tuningCurveSumStatistics :: [V.Vector Double] -> Double
tuningCurveSumStatistics tcs =
    let stcs = V.sum <$> tcs
        (mu,vr) = estimateMeanVariance stcs
     in sqrt vr / mu

shuffleTuningCurves
    :: [V.Vector Double] -> Random s [V.Vector Double]
shuffleTuningCurves tcs = do
    let tcs' = V.concat  tcs
    tcs'' <- Prob $ MWC.uniformShuffle tcs'
    let nn = V.length $ head tcs
        tcs''' = breakEvery nn $ V.toList tcs''
    return $ V.fromList <$> tcs'''

subSampleTuningCurves
    :: Int
    -> [V.Vector Double]
    -> Random s [V.Vector Double]
subSampleTuningCurves nn' tcs = do
    tcs' <- mapM (Prob . MWC.uniformShuffle) tcs
    return $ V.take nn' <$> tcs'



--- Plots ---

subSampleLayout :: [Int] -> [Double] -> Layout Int Double
subSampleLayout nns tcss = execEC $ do

    goalLayout

    let f :: RealFloat x => [x] -> x
        f = parametricExponentialError (zip (realToFrac <$> nns) (realToFrac <$> tcss))
        ab = AD.gradientDescent f [-1,0.1] !! 10000


    plot . liftEC $ do

        plot_lines_style .= solidLine 5 (opaque black)
        plot_lines_values .= [zip nns tcss]

    plot . liftEC $ do

        plot_lines_style .= solidLine 3 (opaque red)
        plot_lines_values .= [zip nns $ parametricExponential ab . fromIntegral <$> nns ]


parametricExponential :: RealFloat x => [x] -> x -> x
parametricExponential [a,b] x = b + exp (a*x)
parametricExponential _ _ = error "Incorrect number of parameters"

parametricExponentialError :: RealFloat x => [(x,x)] -> [x] -> x
parametricExponentialError xys [a,b] =
    average $ square <$> [log (parametricExponential [a,b] x) - log y | (x,y) <- xys]
parametricExponentialError _ _ = error "Incorrect number of parameters"

--- Main ---

showFoo :: Double -> String
showFoo x = showFFloat (Just 5) x ""

analyzeCoefficientOfVariation
    :: forall nn t1 t2
    . (KnownNat nn, KnownNat t1, KnownNat t2, 1 <= nn, 1 <= t1, 1 <= t2)
    => KohnExperiment nn t1 t2
    -> IO (Double,Double)
analyzeCoefficientOfVariation kxp = do

    let kpdr = kohnProjectPath kxp

    (stmstrm0 :: [(Stimulus,M.Map NeuronID [SpikeTime])]) <- read <$> goalReadFile kpdr "stmstrm0"

    let smps0 = L.groupBy (\x y -> fst x == fst y) $ L.sortOn fst stmstrm0
        smps :: [[B.Vector nn Int]]
        smps = map (fromJust . B.fromList . M.elems . fmap length . snd) <$> smps0
        tcs0 :: [Mean # R nn Poisson]
        tcs0 = sufficientStatisticT <$> smps
        tcs = S.fromSized . coordinates <$> tcs0
        cv = tuningCurveSumStatistics tcs
        nn = natValInt (Proxy :: Proxy nn)
        nns = [1..nn]

    tcsss <- realize . replicateM 100 $ mapM (`subSampleTuningCurves` tcs) nns

    let tcss = average <$> L.transpose [ tuningCurveSumStatistics <$> tcss0 | tcss0 <- tcsss ]


    goalRenderableToPNG kpdr "cv-subsampling" 800 600 . toRenderable $ subSampleLayout nns tcss

    putStrLn ""
    putStrLn $ "Experiment: " ++ experiment kxp

    putStrLn "Coefficient of Variation:"
    putStrLn $ showFoo cv

    stcs <- realize $ shuffleTuningCurves tcs

    let scv = tuningCurveSumStatistics stcs
    putStrLn "Shuffled Coefficient of Variation:"
    putStrLn $ showFoo scv

    return (cv,scv)

main :: IO ()
main = do

    (cv0,scv0) <- analyzeCoefficientOfVariation experiment112l44
    (cv1,scv1) <- analyzeCoefficientOfVariation experiment112l45
    (cv2,scv2) <- analyzeCoefficientOfVariation experiment112r35
    (cv3,scv3) <- analyzeCoefficientOfVariation experiment112r36
    (cv4,scv4) <- analyzeCoefficientOfVariation experiment105r62
    (cv5,scv5) <- analyzeCoefficientOfVariation experiment107l114
    (cv6,scv6) <- analyzeCoefficientOfVariation experiment112l16
    (cv7,scv7) <- analyzeCoefficientOfVariation experiment112r32

    putStrLn ""

    putStrLn "Average Small CV:"
    putStrLn . showFoo $ average [cv0,cv1,cv2,cv3]
    putStrLn "Average Small Shuffled CV:"
    putStrLn . showFoo $ average [scv0,scv1,scv2,scv3]
    putStrLn ""

    putStrLn "Average Large CV:"
    putStrLn . showFoo $ average [cv4,cv5,cv6,cv7]
    putStrLn "Average Large Shuffled CV:"
    putStrLn . showFoo $ average [scv4,scv5,scv6,scv7]
    putStrLn ""

    putStrLn "Average CV:"
    putStrLn . showFoo $ average [cv0,cv1,cv2,cv3,cv4,cv5,cv6,cv7]
    putStrLn "Average Shuffled CV:"
    putStrLn . showFoo $ average [scv0,scv1,scv2,scv3,scv4,scv5,scv6,scv7]
    putStrLn ""

    void $ analyzeCoefficientOfVariation small40Pooled

    putStrLn ""

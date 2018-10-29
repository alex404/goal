{-# LANGUAGE BangPatterns,Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import Data.Char
import qualified Data.Vector as V

--- Globals ---


-- Training --

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 10
trnepchn = 2000
trnbrnn = 100
nepchs = 20
trnln = 100000

-- Testing --

rho0n = 1000
tstepch = 1000


--- Functions ---


receptiveFields :: Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli) -> [Layout Int Int]
receptiveFields hrm = do
    n <- [0..nx-1]
    let mx = Mixture # fromList (Replicated Bernoulli nx) (replicate n 0 ++ [1] ++ replicate (nx - 1 - n) 0)
        mz = dualTransition $ czaff >.> mx
        pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- listCoordinates mz]
    return . execEC $ do
        pixMapLayout
        layout_plots .= [pixMapPlot pxss]
    where (Harmonium (Replicated Bernoulli nx) _) = manifold hrm
          czaff = conditionalObservableDistribution hrm

harmoniumRectificationConstantAndError
    :: Int
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s (Double,Double,Double)
harmoniumRectificationConstantAndError n rx rhrm = do
    let (nx,nz,imtx) = splitHarmonium rhrm
        rhox = rx <-> nx
        f x = potential (nz <+> matrixTranspose imtx >.> sufficientStatistic ym x) - rhox <.> sufficientStatistic ym x
    xs <- replicateM n $ standardGenerate rx
    let rho0 = average $ f <$> xs
        se x = (f x - rho0)^2
        rre = divergence (dualTransition rx) nx
    return (rho0, rre, sqrt . average $ se <$> xs)

approximateRectifiedNegativeLogLikelihood
    :: Double
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> [Bool]
    -> Double
approximateRectifiedNegativeLogLikelihood rho0 rx rhrm z =
    let (nx,nz,imtx) = splitHarmonium rhrm
     in negate $ sufficientStatistic zm z <.> nz + potential (nx <+> imtx >.> sufficientStatistic zm z) - potential rx - rho0

filterFun :: Double -> Bool
filterFun l = not (isInfinite l) && not (isNaN l)

mnistNegativeLogLikelihoodHistogram
    :: V.Vector (Mixture :#: Replicated Bernoulli)
    -> Double
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s (Int,Double,Layout Double Int)
mnistNegativeLogLikelihoodHistogram dgs' rho0 rx hrm = do
    zs <- replicateM tstepch $ sampleRandomMNIST dgs'
    let lqs = approximateRectifiedNegativeLogLikelihood rho0 rx hrm <$> zs
        (lqs',rst) = partition filterFun lqs
        n = length lqs'
        mn99 = round $ 0.01 * fromIntegral n
        mx99 = round $ 0.99 * fromIntegral n
        (lqs'',_) = splitAt mx99 lqs'
        (_',lqs''') = splitAt mn99 lqs''
        mn = head lqs'''
        mx = last lqs'''
        lyt = execEC $ do

            goalLayout
            histogramLayout 20 mn mx

            plot . fmap plotBars . liftEC $ histogramPlot 20 mn mx [lqs''']

    return (length rst, average lqs''', lyt)

testEpoch
    :: String
    -> Int
    -> V.Vector (Mixture :#: Replicated Bernoulli)
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> IO ()
testEpoch ttl k dgs' rx hrm = do

    let dr = "mnist/shallow/" ++ map toLower ttl
        hrmstr = "hrm" ++ show k
    (rho0,rre,rerr) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n rx hrm)
    (oob,avg,lllyt) <- runWithSystemRandom (mnistNegativeLogLikelihoodHistogram dgs' rho0 rx hrm)
    putStrLn $ ttl ++ ":\n"
    putStrLn $ "Rectification Constant: " ++ show rho0
    putStrLn $ "Average Rectification Distance: " ++ show rre
    putStrLn $ "Mean Squared Rectification Error: " ++ show rerr
    putStrLn $ "Average -Log-Likelihood: " ++ show avg
    putStrLn $ "Infinite -Log-Likelihoods: " ++ show oob ++ "\n"
    goalRenderableToSVG dr (hrmstr ++ "-" ++ "negative-log-likelihoods") 1200 600 $ toRenderable lllyt
    sequence_ [goalRenderableToSVG (dr ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ receptiveFields hrm]

harmoniumTests dgs' (k,((_,rx,rhrm),(_,sx,shrm))) = do

    testEpoch "Rectified" k dgs' rx rhrm
    testEpoch "Standard" k dgs' sx shrm
    --testEpoch "Twisted" k dgs' tx thrm
    goalWriteFile "mnist/shallow" "rxrhrmsxshrm" $ show ((rx,rhrm),(sx,shrm))

--- Main ---


main = do

    dgs <- mnistUnsupervisedTrainingData
    --zss <- cycle . breakEvery nbtch <$> runWithSystemRandom (mapM standardGenerate dgs)

    nx0 <- runWithSystemRandom $ initialize w0 ym
    nz0 <- runWithSystemRandom $ initialize w0 zm

    let imtx0 = zero (Tensor ym zm)
        rx0 = nx0
        hrm0 = joinHarmonium nx0 nz0 imtx0

    let mlfun (zs',rx,hrm) = rectifiedHarmoniumRelativeEntropyDifferentials zs' rx hrm

    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ sampleRandomMNIST dgs))

    cdzmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (bulkContrastiveDivergence cdn))
    cdymly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectificationDifferentials cdn nbtch))

    rctfrmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectifierDifferentials nbtch))
    rctmlmly <- runWithSystemRandom (accumulateRandomFunction0 mlfun)

    let rmly = accumulateMealy0 (0,rx0,hrm0) $ proc (zs,(k,!rx,!hrm)) -> do
            drx <- rctfrmly -< (rx,hrm)
            dhrmx <- cdymly -< (rx,hrm)
            dhrmz <- rctmlmly -< (zs,rx,hrm)
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            hrm' <- adamDescent eps bt1 bt2 rg -< (dhrmx <+> dhrmz)
            returnA -< (k+1,rx',trace (show k) hrm')

    let smly = accumulateMealy0 (0,rx0,hrm0) $ proc (zs,(k,!rx,!hrm)) -> do
            drx <- rctfrmly -< (rx,hrm)
            dhrm <- cdzmly -< (zs,hrm)
            hrm' <- adamDescent eps bt1 bt2 rg -< dhrm
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            returnA -< (k+1,rx',trace (show k) hrm')

{-
    let tmly = accumulateMealy0 (rx0,hrm0) $ proc (zs,!(!rx,!hrm)) -> do
            drx <- rctfrmly -< (rx,hrm)
            dhrmz <- rctmlmly -< (zs,rx,hrm)
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            hrm' <- adamDescent eps bt1 bt2 rg -< dhrmz
            returnA -< (rx',hrm')
-}

    dgs' <- mnistUnsupervisedTrainingData

    let rxrhrmsxshrm = head . drop trnln . streamChain $ (rmly &&& smly) <<< dgchn

    harmoniumTests dgs' (100,rxrhrmsxshrm)

{- Graveyard

testStandardEpoch
    :: Int
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> IO ()
testStandardEpoch k sx hrm = do
    let dr = "mnist/standard-harmonium"
        hrmstr = "hrm" ++ show k
    (rho0,rre,rerr) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n sx hrm)
    putStrLn "Standard:\n"
    putStrLn $ "Rectification Constant: " ++ show rho0
    putStrLn $ "Average Rectification Distance: " ++ show rre
    putStrLn $ "Mean Squared Rectification Error: " ++ show rerr
    putStrLn $ "Average Rectification Error: " ++ show rerr ++ "\n"
    sequence_ [goalRenderableToSVG (dr ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ receptiveFields hrm]

-}

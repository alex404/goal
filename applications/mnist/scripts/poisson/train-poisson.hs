{-# LANGUAGE BangPatterns,Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import PMNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Data.Vector as V

--- Globals ---

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 10
trnepchn = 1000
trnbrnn = 0
nepchs = 21

-- Testing --

rho0n = 1000
tstepch = 1000
nmrg = 100

--- Functions ---

likelihoodRenderable anlls = toRenderable . execEC $ do

   goalLayout

   layout_x_axis . laxis_title .= "Epoch"
   layout_y_axis . laxis_title .= "-Log-Likelihood"

   plot . liftEC $ do

       plot_lines_style .= solidLine 3 (opaque black)
       plot_lines_values .= [zip [0..nepchs] anlls]



testEpoch
    :: V.Vector (Mixture :#: Replicated Poisson)
    -> String
    -> (Int,Natural :#: Replicated Bernoulli,Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson))
    -> IO Double
testEpoch dgs' ttl (k,ry,hrm1) = do

    zs <- runWithSystemRandom . replicateM tstepch $ sampleRandomMNIST dgs'
    (rho0,rre,rerr) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n ry hrm1)
    mflds <- runWithSystemRandom $ marginalReceptiveFields nmrg ry hrm1
    let avg = average $ approximateRectifiedNegativeLogLikelihood rho0 ry hrm1 <$> zs
        epch = div k trnepchn
        hrmstr = "hrm" ++ show epch
    goalWriteFile (mnstlyr1dr ++ "/" ++ ttl) ("ryhrm" ++ show epch) $ show (k,ry,hrm1)
    putStrLn $ ttl ++ " epoch " ++ show epch ++ ":\n"
    putStrLn $ "Rectification Constant: " ++ show rho0
    putStrLn $ "Average Rectification Distance: " ++ show rre
    putStrLn $ "Mean Squared Rectification Error: " ++ show rerr
    putStrLn $ "Average -Log-Likelihood: " ++ show avg ++ "\n"
    sequence_ [goalRenderableToSVG (mnstlyr1dr ++ "/" ++ ttl ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ receptiveFields hrm1]
    sequence_ [goalRenderableToSVG (mnstlyr1dr ++ "/" ++ ttl ++ "/marginal-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] mflds]
    return avg

testEpochs dgs' ((rk,ry,rhrm),(sk,sy,shrm)) = do
    rnll <- testEpoch dgs' "rectified" (rk,ry,rhrm)
    snll <- testEpoch dgs' "standard" (sk,sy,shrm)
    return (rnll,snll)

--- Main ---


main = do

    dgs <- mnistUnsupervisedTrainingData
    --zss <- cycle . breakEvery nbtch <$> runWithSystemRandom (mapM standardGenerate dgs)

    ny0 <- runWithSystemRandom $ initialize w0 ym
    nz0 <- runWithSystemRandom $ initialize w0 zm

    let nyz0 = zero (Tensor ym zm)
        ry0 = ny0
        hrm0 = joinHarmonium ny0 nz0 nyz0
        {-
    rfl <- goalReadFile mnstlyr1dr "rectified/ryhrm100"
    sfl <- goalReadFile mnstlyr1dr "standard/ryhrm100"
    rnllsfl <- goalReadFile mnstlyr1dr "rectified/rnlls"

    let (k0,rry0,rnyz0) = read rfl
    let (k0',sry0,snyz0) = read sfl
    let rnlls0 = read rnllsfl
    print $ (k0' :: Int)
    -}
    let k0 = 0
        rnlls0 = []

    let mlfun (zs',ry,hrm1) = rectifiedHarmoniumRelativeEntropyDifferentials zs' ry hrm1

    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ sampleRandomMNIST dgs))

    cdymly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectificationDifferentials cdn nbtch))
    cdzmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (bulkContrastiveDivergence cdn))

    rctfrmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectifierDifferentials nbtch))
    rctmlmly <- runWithSystemRandom (accumulateRandomFunction0 mlfun)

    let rmly = accumulateMealy0 (k0,ry0,hrm0) $ proc (zs,(k,!ry,!hrm1)) -> do
            dry <- rctfrmly -< (ry,hrm1)
            dhrm1y <- cdymly -< (ry,hrm1)
            dhrm1z <- rctmlmly -< (zs,ry,hrm1)
            ry' <- adamDescent eps bt1 bt2 rg -< dry
            hrm1' <- adamDescent eps bt1 bt2 rg -< (dhrm1y <+> dhrm1z)
            returnA -< (k+1,ry',hrm1')

    let smly = accumulateMealy0 (k0,ry0,hrm0) $ proc (zs,(k,!ry,!hrm1)) -> do
            dry <- rctfrmly -< (ry,hrm1)
            dhrm1y <- cdymly -< (ry,hrm1)
            dhrm1z <- cdzmly -< (zs,hrm1)
            ry' <- adamDescent eps bt1 bt2 rg -< dry
            hrm1' <- adamDescent eps bt1 bt2 rg -< (dhrm1y <+> dhrm1z)
            returnA -< (k+1,ry',hrm1')

    dgs' <- mnistUnsupervisedTestData

    let kryrhrmsyshrm1s = take nepchs $ takeEvery trnepchn . drop trnbrnn . streamChain $ (rmly &&& smly) <<< dgchn

    nllss <- mapM (testEpochs dgs') kryrhrmsyshrm1s
    let (rnlls1,_) = unzip nllss
    let rnlls = rnlls0 ++ rnlls1

    goalWriteFile (mnstlyr1dr ++ "/rectified") "rnlls" $ show rnlls
    goalRenderableToSVG (mnstlyr1dr ++ "/rectified") "rnlls" 500 250 $ likelihoodRenderable rnlls

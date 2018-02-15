{-# LANGUAGE BangPatterns,Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

import qualified Data.Vector as V

--- Globals ---


-- Training --

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 50
trnepchn = 10000
trnbrnn = 100
nepchs = 10

-- Testing --

rho0n = 1000
tstepch = 1000
nmrgs = 100

--- Functions ---


testEpoch
    :: V.Vector (Mixture :#: Replicated Bernoulli)
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> (Int,Natural :#: Replicated Bernoulli,Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli))
    -> IO Double
testEpoch dgs' ry hrm1 (k,rx,hrm2) = do

    ys <- runWithSystemRandom . replicateM tstepch $ resampleRandomMNIST dgs' hrm1
    (rho0,rre,rerr) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n rx hrm2)
    let (ny1,nz,nyz) = splitHarmonium hrm1
        (nx,ny2,nxy) = splitHarmonium hrm2
        dhrm = joinDeepHarmonium nx (ny2 <+> ry <-> ny1) nz nxy nyz
    rflds <- runWithSystemRandom $ deepReceptiveFields ry dhrm
    mflds <- runWithSystemRandom $ deepMarginalReceptiveFields nmrgs rx ry dhrm
    let avg = average $ approximateRectifiedNegativeLogLikelihood rho0 rx hrm2 <$> ys
        epch = div k trnepchn
        hrmstr = "hrm" ++ show epch
    putStrLn $ "Epoch " ++ show epch ++ ":\n"
    putStrLn $ "Rectification Constant: " ++ show rho0
    putStrLn $ "Average Rectification Distance: " ++ show rre
    putStrLn $ "Mean Squared Rectification Error: " ++ show rerr
    putStrLn $ "Average -Log-Likelihood: " ++ show avg ++ "\n"
    sequence_ [goalRenderableToSVG (mnstlyr2dr ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ rflds]
    sequence_ [goalRenderableToSVG (mnstlyr2dr ++ "/marginal-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] mflds]
    return avg

--- Main ---


main = do

    hrm1str <- goalReadFile mnstlyr1dr "kryhrm1"
    let (k1,ry,hrm1) = read hrm1str
        (ny1,_,_) = splitHarmonium hrm1
    --zss <- cycle . breakEvery nbtch <$> runWithSystemRandom (mapM standardGenerate dgs)

    nx0 <- runWithSystemRandom $ initialize w0 xm
    print (k1 :: Int)

    let ny0 = ry <+> ny1
        xy0 = zero (Tensor xm ym)
        rx0 = nx0
        hrm20 = joinHarmonium nx0 ny0 xy0

    let mlfun (ys',rx,hrm2) = rectifiedHarmoniumRelativeEntropyDifferentials ys' rx hrm2

    dgs <- mnistUnsupervisedTrainingData
    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ resampleRandomMNIST dgs hrm1))
    cdxmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectificationDifferentials cdn nbtch))
    rctfrmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectifierDifferentials nbtch))
    rctmlmly <- runWithSystemRandom (accumulateRandomFunction0 mlfun)

    let rmly = accumulateMealy0 (0,rx0,hrm20) $ proc (ys,!(k,!rx,!hrm2)) -> do
            drx <- rctfrmly -< (rx,hrm2)
            dhrm2x <- cdxmly -< (rx,hrm2)
            dhrm2y <- rctmlmly -< (ys,rx,hrm2)
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            hrm2' <- adamDescent eps bt1 bt2 rg -< (dhrm2x <+> dhrm2y)
            returnA -< (k+1,rx',hrm2')

    let krxhrm2s = take nepchs $ takeEvery trnepchn . drop trnbrnn . streamChain $ rmly <<< dgchn

    dgs' <- mnistUnsupervisedTrainingData
    nlls <- mapM (testEpoch dgs' ry hrm1) krxhrm2s

    goalWriteFile mnstlyr2dr "krxhrm2" . show $ last krxhrm2s
    goalWriteFile mnstlyr2dr "nlls" $ show nlls

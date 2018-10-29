{-# LANGUAGE BangPatterns,Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation

--- Globals ---


-- Training --

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 10
trnepchn = 10000
trnbrnn = 100
nepchs = 10

-- Testing --

rho0n = 1000
tstepch = 1000
nmrgs = 100


--- Functions ---


testEpoch (k,rx,ry,dhrm) = do

    let (nx,ny,nz,nxy,nyz) = splitDeepHarmonium dhrm
        hrm1 = joinHarmonium ny nz nyz
        hrm2 = joinHarmonium nx ry nxy
    (rho0x,rrex,_) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n rx hrm2)
    (rho0y,rrey,_) <- runWithSystemRandom (harmoniumRectificationConstantAndError rho0n ry hrm1)
    rflds <- runWithSystemRandom $ deepReceptiveFields ry dhrm
    mflds <- runWithSystemRandom $ deepMarginalReceptiveFields nmrgs rx ry dhrm
    let epch = div k trnepchn
        hrmstr = "hrm" ++ show epch
    putStrLn $ "Epoch " ++ show epch ++ ":\n"
    putStrLn $ "X Rectification Constant: " ++ show rho0x
    putStrLn $ "X Average Rectification Distance: " ++ show rrex
    putStrLn $ "Y Rectification Constant: " ++ show rho0y
    putStrLn $ "Y Average Rectification Distance: " ++ show rrey
    sequence_ [goalRenderableToSVG (mnstdpftdr ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ rflds]
    sequence_ [goalRenderableToSVG (mnstdpftdr ++ "/marginal-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] mflds]
    return ()

--- Main ---


main = do

    dgs <- mnistUnsupervisedTrainingData
    --zss <- cycle . breakEvery nbtch <$> runWithSystemRandom (mapM standardGenerate dgs)

    nx0 <- runWithSystemRandom $ initialize w0 xm
    ny0 <- runWithSystemRandom $ initialize w0 ym
    nz0 <- runWithSystemRandom $ initialize w0 zm

    let nxy0 = zero (Tensor xm ym)
        nyz0 = zero (Tensor ym zm)
        rx0 = nx0
        ry0 = ny0
        dhrm0 = joinDeepHarmonium nx0 ny0 nz0 nxy0 nyz0

    let rctmlfun (zs',rx,ry,dhrm) = dualContrastiveDivergence cdn zs' rx ry dhrm
        rctfrfun (rx,ry,dhrm) = strongRectifierDifferentials nbtch rx ry dhrm
        rctxyfun (rx,ry,dhrm) = strongRectificationDifferentials cdn nbtch rx ry dhrm

    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ sampleRandomMNIST dgs))

    rctmlmly <- runWithSystemRandom (accumulateRandomFunction0 rctmlfun)
    rctfrmly <- runWithSystemRandom (accumulateRandomFunction0 rctfrfun)
    rctxymly <- runWithSystemRandom (accumulateRandomFunction0 rctxyfun)

    let rmly = accumulateMealy0 (0,rx0,ry0,dhrm0) $ proc (zs,!(k,!rx,!ry,!dhrm)) -> do
            (drx,dry) <- rctfrmly -< (rx,ry,dhrm)
            ddhrmxy <- rctxymly -< (rx,ry,dhrm)
            ddhrmz <- rctmlmly -< (zs,rx,ry,dhrm)
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            ry' <- adamDescent eps bt1 bt2 rg -< dry
            dhrm' <- adamDescent eps bt1 bt2 rg -< (ddhrmxy <+> ddhrmz)
            returnA -< (k+1,rx',ry',trace (show k) dhrm')

    let krxryhrms = take nepchs $ takeEvery trnepchn . drop trnbrnn . streamChain $ rmly <<< dgchn
    mapM_ testEpoch  krxryhrms
    --harmoniumTests dgs' (100,rxrhrmsxshrm)

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

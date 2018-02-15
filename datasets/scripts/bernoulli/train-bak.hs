{-# LANGUAGE Arrows,TypeOperators,FlexibleContexts #-}

--- Imports ---


import MNIST

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Simulation


--- Globals ---

eps = 0.001
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 10
w0 = Standard # fromList Normal [0,0.0001]
hrmxs = 64
zm = Replicated Bernoulli $ dghght * dgwdth
xm = Replicated Bernoulli hrmxs
hrmm = Harmonium xm zm

--- Functions ---

{-

mnistTest
    :: [(Int, Mixture :#: Replicated Bernoulli)]
    -> Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> [(Int,[Double])]
mnistTest ldgs hrm =
    let aff = conditionalLatentDistribution hrm
        pss = do
            n <- xcats
            let grp = snd <$> filter ((==n) . fst) ldgs
                ps = map average . transpose $ (\p -> density p <$> xcats) <$> (aff >$> grp)
            return ps
     in zip xcats $ transpose pss
-}

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

mnistMaximumLikelihood
    :: [Mixture :#: Replicated Bernoulli]
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s Double
mnistMaximumLikelihood mzs rx hrm = do
    zs <- mapM standardGenerate mzs
    xs <- replicateM (length mzs) $ standardGenerate rx
    let qz = averagePoint $ conditionalObservableDistribution hrm >$>* xs
    return . average $ log . density qz <$> zs


--- Main ---


main = do

    ldgs <- mnistTrainingData
    let (_,dgs) = unzip ldgs

    nx0 <- runWithSystemRandom $ initialize w0 xm
    nz0 <- runWithSystemRandom $ initialize w0 zm

    let imtx0 = zero (Tensor xm zm)
        rx0 = nx0
        hrm0 = joinHarmonium nx0 nz0 imtx0

    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ standardGenerate =<< randomElement dgs))

{-
    cdzchn <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (bulkContrastiveDivergence cdn))

    let cdtrnchn = accumulateMealy0 hrm0 $ proc ((),hrm) -> do
            zs <- dgchn -< ()
            dhrm <- cdzchn -< (zs,hrm)
            adamDescent eps bt1 bt2 rg -< dhrm
            --vanillaGradientDescent eps -< dhrm
            --}

    let mlfun (zs',rx,hrm) = rectifiedHarmoniumRelativeEntropyDifferentials zs' rx hrm

    cdxchn <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectificationDifferentials cdn nbtch))
    mlchn <- runWithSystemRandom (accumulateRandomFunction0 mlfun)
    rctchn <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (rectifierDifferentials nbtch))
    --ldgs' <- mnistTestData
    --let mnistTester = mnistTest ldgs'
    --
    let rcttrnchn = accumulateMealy0 (rx0,hrm0) $ proc ((),(rx,hrm)) -> do
            zs <- dgchn -< ()
            drx <- rctchn -< (rx,hrm)
            dhrmx <- cdxchn -< (rx,hrm)
            dhrmz <- mlchn -< (zs,rx,hrm)
            --dhrmz <- cdzchn -< (zs,hrm)
            rx' <- adamDescent eps bt1 bt2 rg -< drx
            hrm' <- adamDescent eps bt1 bt2 rg -< (dhrmx <+> dhrmz)
            returnA -< (rx',hrm')

    --ldgs' <- mnistTestData
    --let mlTest = mnistMaximumLikelihood (snd <$> ldgs')

    let hrms = snd <$> streamChain rcttrnchn
    let (hrm1:hrms') = drop 100 hrms
    let (hrm2:hrms'') = drop 1000 hrms'
    let (hrm3:_) = drop 10000 hrms''
        {-
    tst1 <- runWithSystemRandom $ mlTest rx1 hrm1
    print tst1
    tst2 <- runWithSystemRandom $ mlTest rx2 hrm2
    print tst2
    tst3 <- runWithSystemRandom $ mlTest rx3 hrm3
    print tst3
    -}

    sequence_ [goalRenderableToSVG "mnist" ("hrm1-" ++ show n) 400 400 $ toRenderable lyt | (n,lyt) <- zip [0..] $ receptiveFields hrm1]
    sequence_ [goalRenderableToSVG "mnist" ("hrm2-" ++ show n) 400 400 $ toRenderable lyt | (n,lyt) <- zip [0..] $ receptiveFields hrm2]
    sequence_ [goalRenderableToSVG "mnist" ("hrm3-" ++ show n) 400 400 $ toRenderable lyt | (n,lyt) <- zip [0..] $ receptiveFields hrm3]

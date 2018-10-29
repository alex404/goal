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


-- Manifolds --

xcat = Categorical [0..9]
xzm = Harmonium xcat zm

-- Training --

eps = 0.01
bt1 = 0.9
bt2 = 0.999
rg = 1e-8
nbtch = 10
cdn = 10
trnepchn = 1000
trnbrnn = 100
nepchs = 20

-- Testing --

mln = 10000

--- Functions ---

bulkCategoricalTrainer
    :: [[Bool]] -- ^ The observations
    -> (Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)) -- ^ The harmonium
    -> RandST s (Differentials :#: Tangent Natural (Harmonium (Categorical [Int]) (Replicated Bernoulli))) -- ^ The resulting stochastic gradient
bulkCategoricalTrainer zs' hrm = do
    xs' <- mapM standardGenerate $ conditionalLatentDistribution hrm >$>* zs'
    xzs <- sampleCategoricalHarmonium (length zs') hrm
    let Harmonium x z = manifold hrm
        (xs,zs) = unzip xzs
        mxs = sufficientStatistic x <$> xs
        mxs' = sufficientStatistic x <$> xs'
        mzs = sufficientStatistic z <$> zs
        mzs' = sufficientStatistic z <$> zs'
    return $ harmoniumGradientCalculator mxs mxs' mzs mzs' hrm


categoricalReceptiveFields  :: Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli) -> [Layout Int Int]
categoricalReceptiveFields hrm = do
    n <- [0..9]
    let mz = dualTransition $ conditionalObservableDistribution hrm >.>* n
        pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- listCoordinates mz]
    return . execEC $ do
        pixMapLayout
        layout_plots .= [pixMapPlot pxss]

categoricalNegativeLogLikelihood
    :: Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> [Bool]
    -> Double
categoricalNegativeLogLikelihood hrm z =
    let (rho0,rx) = categoricalHarmoniumRectificationParameters hrm
        (nx,nz,nxz) = splitHarmonium hrm
     in negate $ sufficientStatistic zm z <.> nz + potential (nx <+> nxz >.> sufficientStatistic zm z) - potential rx - rho0

averageCategoricalNegativeLogLikelihood
    :: Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> V.Vector (Mixture :#: Replicated Bernoulli)
    -> RandST s Double
averageCategoricalNegativeLogLikelihood hrm dgs = do
    dgs' <- replicateM mln $ randomElement' dgs
    blss <- mapM standardGenerate dgs'
    return . average $ categoricalNegativeLogLikelihood hrm <$> blss

categoricalWeightHistogram
    :: Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> Layout Double Double
categoricalWeightHistogram hrm = execEC $ do

    let (nx,nz,nxz) = splitHarmonium hrm

    plot . fmap plotBars . liftEC $ do
        logHistogramPlot0 2 20 [listCoordinates nx, listCoordinates nz, listCoordinates nxz]
        plot_bars_titles .= ["Latent Biases", "Visible Biases", "Interactions"]
        --plot_bars_item_styles .= [(solidFillStyle $ opaque blue, Nothing)]

categoricalRectificationDistance
    :: Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> Double
categoricalRectificationDistance hrm = do
    let rx = snd $ categoricalHarmoniumRectificationParameters hrm
        (nx,_,_) = splitHarmonium hrm
     in divergence (dualTransition rx) nx

testEpoch
    :: String
    -> Int
    -> V.Vector (Mixture :#: Replicated Bernoulli)
    -> Natural :#: Harmonium (Categorical [Int]) (Replicated Bernoulli)
    -> IO ()
testEpoch ttl k dgs' hrm = do

    let rho0 = fst $ categoricalHarmoniumRectificationParameters hrm
        rerr = categoricalRectificationDistance hrm
        dr = "mnist/categorical/" ++ map toLower ttl
        hrmstr = "hrm" ++ show k
    avg <- runWithSystemRandom $ averageCategoricalNegativeLogLikelihood hrm dgs'
    putStrLn $ ttl ++ ":\n"
    putStrLn $ "Rectification Constant: " ++ show rho0
    putStrLn $ "Rectification Distance: " ++ show rerr
    putStrLn $ "Average -Log-Likelihood: " ++ show avg ++ "\n"
    goalRenderableToSVG dr (hrmstr ++ "-" ++ "weights") 1200 600 . toRenderable $ categoricalWeightHistogram hrm
    sequence_ [goalRenderableToSVG (dr ++ "/receptive-fields") (hrmstr ++ "-" ++ show n) 400 400 $ toRenderable lyt
      | (n,lyt) <- zip [0..] $ categoricalReceptiveFields hrm]

harmoniumTests dgs' (k,(rcthrm,cdnhrm)) = do

    putStrLn $ "Epoch " ++ show k ++ ":\n"

    testEpoch "Rectified" k dgs' rcthrm
    testEpoch "Standard" k dgs' cdnhrm
    --goalWriteFile "mnist/categorical" "rxrhrmsxshrm" $ show ((rx,rhrm),(sx,shrm))

--- Main ---


main = do

    dgs <- mnistUnsupervisedTrainingData
    --zss <- cycle . breakEvery nbtch <$> runWithSystemRandom (mapM standardGenerate dgs)

    hrm0 <- runWithSystemRandom $ initialize w0 xzm


    dgchn <- runWithSystemRandom (accumulateRandomFunction0 $ const (replicateM nbtch $ sampleRandomMNIST dgs))
    rmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (bulkCategoricalTrainer))
    cmly <- runWithSystemRandom (accumulateRandomFunction0 $ uncurry (bulkContrastiveDivergence cdn))


    let rctmly = accumulateMealy0 hrm0 $ rmly >>> adamDescent eps bt1 bt2 rg
    let cdnmly = accumulateMealy0 hrm0 $ cmly >>> adamDescent eps bt1 bt2 rg

    dgs' <- mnistUnsupervisedTrainingData

    let hrmss = take nepchs $ takeEvery trnepchn . drop trnbrnn . streamChain $ (rctmly &&& cdnmly) <<< dgchn

    mapM_ (harmoniumTests dgs') $ zip [0..] hrmss

    --harmoniumTests dgs' (100,rxrhrmsxshrm)

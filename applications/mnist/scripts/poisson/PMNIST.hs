{-# LANGUAGE TypeOperators #-}

module PMNIST where

import Goal.Core
import Goal.Geometry
import Goal.Probability

import Data.IDX
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V



--- Harmoniums ---


-- Initialization --

hrmxs = 256
hrmys = 12^2
xm = Replicated Bernoulli hrmxs
ym = Replicated Bernoulli hrmys
zm = Replicated Poisson $ dghght * dgwdth
w0 = Standard # fromList Normal [0,0.001]
scl = 1/255


--- MNIST ---

-- IO --

mnstdr0 = "mnist"
mnstdr = "pmnist"
mnstlyr1dr = mnstdr ++ "/deep/layer1"
mnstlyr2dr = mnstdr ++ "/deep/layer2"
mnstdpftdr = mnstdr ++ "/deep/fine-tune"
trnlblfl = "train-labels-idx1-ubyte"
trnimgfl = "train-images-idx3-ubyte"
tstlblfl = "t10k-labels-idx1-ubyte"
tstimgfl = "t10k-images-idx3-ubyte"

dghght = 28
dgwdth = 28

mnistUnsupervisedTrainingData :: IO (V.Vector (Mixture :#: Replicated Poisson))
mnistUnsupervisedTrainingData = do
    ldgs <- mnistData trnlblfl trnimgfl
    return $ snd <$> ldgs

mnistTrainingData :: IO (V.Vector (Int, Mixture :#: Replicated Poisson))
mnistTrainingData = mnistData trnlblfl trnimgfl

mnistUnsupervisedTestData :: IO (V.Vector (Mixture :#: Replicated Poisson))
mnistUnsupervisedTestData = do
    ldgs <- mnistData tstlblfl tstimgfl
    return $ snd <$> ldgs

mnistTestData :: IO (V.Vector (Int, Mixture :#: Replicated Poisson))
mnistTestData = mnistData tstlblfl tstimgfl

mnistData :: String -> String -> IO (V.Vector (Int, Mixture :#: Replicated Poisson))
mnistData lblfl imgfl = do

    lblpth <- goalDataLocation mnstdr0 lblfl
    imgpth <- goalDataLocation mnstdr0 imgfl
    mlbls <- decodeIDXLabelsFile lblpth
    mimgs <- decodeIDXFile imgpth

    let (lbls,dgs) = unzip . fromJust $ labeledDoubleData (fromJust mlbls) (fromJust mimgs)
        n = dghght * dgwdth

    return . V.fromList $ zip lbls [fromList (Replicated Poisson n) . U.toList $ U.map (*scl) dg | dg <- dgs]

digitToPixMap :: Mixture :#: Replicated Poisson -> [[AlphaColour Double]]
digitToPixMap dg =
    let scl2 = maximum $ listCoordinates dg
     in breakEvery dgwdth [ opaque $ rgb px px px | px <- (/scl2) <$> listCoordinates dg ]

sampleRandomLabelledMNIST :: V.Vector (Int, Mixture :#: Replicated Poisson) -> RandST s (Int,[Int])
sampleRandomLabelledMNIST ldgs = do
    (l,mz) <- randomElement' ldgs
    dg <- standardGenerate mz
    return (l,dg)

sampleHarmoniumLabelledMNIST
    :: Int
    -> V.Vector (Int, Mixture :#: Replicated Poisson)
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson)
    -> RandST s [(Int,[Bool])]
sampleHarmoniumLabelledMNIST nbtch ldgs hrm = do
            lblss <- replicateM nbtch $ sampleRandomLabelledMNIST ldgs
            let (ls,blss) = unzip lblss
            blss' <- mapM standardGenerate $ conditionalLatentDistribution hrm >$>* blss
            return $ zip ls blss'

sampleRandomMNIST :: V.Vector (Mixture :#: Replicated Poisson) -> RandST s [Int]
sampleRandomMNIST dgs = standardGenerate =<< randomElement' dgs

resampleRandomMNIST
    :: V.Vector (Mixture :#: Replicated Poisson)
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson)
    -> RandST s [Bool]
resampleRandomMNIST dgs hrm = do
    dg <-  randomElement' dgs
    bls <- standardGenerate dg
    standardGenerate $ conditionalLatentDistribution hrm >.>* bls

marginalReceptiveFields
    :: Int
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson)
    -> RandST s [Layout Int Int]
marginalReceptiveFields nbtch ry hrm = do
    ys <- replicateM nbtch $ standardGenerate ry
    let mzs = conditionalObservableDistribution hrm >$>* ys
    return $ do
        mz <- mzs
        let pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- (/maximum (listCoordinates mz)) <$> listCoordinates mz]
        return . execEC $ do
            pixMapLayout
            layout_plots .= [pixMapPlot pxss]

receptiveFields :: Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson) -> [Layout Int Int]
receptiveFields hrm = do
    n <- [0..nx-1]
    let mx = Mixture # fromList (Replicated Bernoulli nx) (replicate n 0 ++ [1] ++ replicate (nx - 1 - n) 0)
        mz = dualTransition $ czaff >.> mx
        pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- (/maximum (listCoordinates mz)) <$> listCoordinates mz]
    return . execEC $ do
        pixMapLayout
        layout_plots .= [pixMapPlot pxss]
    where (Harmonium (Replicated Bernoulli nx) _) = manifold hrm
          czaff = conditionalObservableDistribution hrm
{-

deepMarginalReceptiveFields
    :: Int
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: DeepHarmonium (Replicated Bernoulli) (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s [Layout Int Int]
deepMarginalReceptiveFields nbtch rx ry dhrm = do
    xs <- replicateM nbtch $ standardGenerate rx
    let (_,_,nz,nxy,nyz) = splitDeepHarmonium dhrm
    ys <- sequence $ standardGenerate <$> joinAffine ry (matrixTranspose nxy) >$>* xs
    let mzs = dualTransition <$> joinAffine nz (matrixTranspose nyz) >$>* ys
    return $ do
        mz <- mzs
        let pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- listCoordinates mz]
        return . execEC $ do
            pixMapLayout
            layout_plots .= [pixMapPlot pxss]

deepReceptiveFields
    :: Natural :#: Replicated Bernoulli
    -> Natural :#: DeepHarmonium (Replicated Bernoulli) (Replicated Bernoulli) (Replicated Bernoulli)
    -> RandST s [Layout Int Int]
deepReceptiveFields ry dhrm = do
    let (DeepHarmonium (Replicated Bernoulli nx) _ _) = manifold dhrm
        mxs = [ Mixture # fromList (Replicated Bernoulli nx) (replicate n 0 ++ [1] ++ replicate (nx - 1 - n) 0) | n <- [0..nx-1] ]
        (_,_,nz,nxy,nyz) = splitDeepHarmonium dhrm
    ys <- sequence $ standardGenerate <$> joinAffine ry (matrixTranspose nxy) >$> mxs
    let mzs = dualTransition <$> joinAffine nz (matrixTranspose nyz) >$>* ys
    return $ do
        mz <- mzs
        let pxss = breakEvery dgwdth [opaque $ rgb px px px | px <- listCoordinates mz]
        return . execEC $ do
            pixMapLayout
            layout_plots .= [pixMapPlot pxss]

-}
harmoniumRectificationConstantAndError
    :: Int
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson)
    -> RandST s (Double,Double,Double)
harmoniumRectificationConstantAndError n rx rhrm = do
    let (nx,nz,imtx) = splitHarmonium rhrm
        rhox = rx <-> nx
        f x = potential (nz <+> matrixTranspose imtx >.> sufficientStatistic xm x) - rhox <.> sufficientStatistic xm x
    xs <- replicateM n $ standardGenerate rx
    let rho0 = average $ f <$> xs
        se x = (f x - rho0)^2
        rre = divergence (dualTransition rx) nx
    return (rho0, rre, sqrt . average $ se <$> xs)

approximateRectifiedNegativeLogLikelihood
    :: Double
    -> Natural :#: Replicated Bernoulli
    -> Natural :#: Harmonium (Replicated Bernoulli) (Replicated Poisson)
    -> [Int]
    -> Double
approximateRectifiedNegativeLogLikelihood rho0 rx rhrm z =
    let (nx,nz,imtx) = splitHarmonium rhrm
     in negate $ sufficientStatistic zm z <.> nz + potential (nx <+> imtx >.> sufficientStatistic zm z) - potential rx - rho0



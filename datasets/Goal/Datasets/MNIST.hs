{-# LANGUAGE TypeOperators,DataKinds #-}
module Goal.Datasets.MNIST where

import Goal.Core

import qualified Goal.Core.Vector.Boxed as B
import qualified Data.Vector.Generic as G

import Data.IDX



--- Harmoniums ---


-- Initialization --

type MNISTHeight = 28
type MNISTWidth = 28
type MNISTSize = MNISTHeight * MNISTWidth
--type Digit m = Replicated (Height * Width) m


--- MNIST ---

mnstdr,trnlblfl,trnimgfl,tstlblfl,tstimgfl :: String
mnstdr = "mnist"
trnlblfl = "train-labels-idx1-ubyte"
trnimgfl = "train-images-idx3-ubyte"
tstlblfl = "t10k-labels-idx1-ubyte"
tstimgfl = "t10k-images-idx3-ubyte"

-- IO --

mnistData :: String -> String -> IO [(B.Vector MNISTSize Double, Int)]
mnistData lblfl imgfl = do

    lblpth <- goalDatasetLocation mnstdr lblfl
    imgpth <- goalDatasetLocation mnstdr imgfl
    mlbls <- decodeIDXLabelsFile lblpth
    mimgs <- decodeIDXFile imgpth

    let (lbls,dgs) = unzip . fromJust $ labeledIntData (fromJust mlbls) (fromJust mimgs)

    return $ zip (fmap ((/255) . fromIntegral) . fromJust . B.toSized . G.convert <$> dgs) lbls

mnistTrainingData :: IO [(B.Vector MNISTSize Double, Int)]
mnistTrainingData = mnistData trnlblfl trnimgfl

mnistTestData :: IO [(B.Vector MNISTSize Double, Int)]
mnistTestData = mnistData tstlblfl tstimgfl

{-
mnistUnsupervisedTrainingData :: IO (V.Vector (Mean # Replicated n Bernoulli))
mnistUnsupervisedTrainingData = do
    ldgs <- mnistData trnlblfl trnimgfl
    return $ snd <$> ldgs

mnistUnsupervisedTestData :: IO (V.Vector (Mean # Replicated n Bernoulli))
mnistUnsupervisedTestData = do
    ldgs <- mnistData tstlblfl tstimgfl
    return $ snd <$> ldgs

digitToPixMap :: Mean # Replicated n Bernoulli -> [[AlphaColour Double]]
digitToPixMap dg = breakEvery dgwdth [ opaque $ rgb px px px | px <- listCoordinates dg ]

sampleRandomLabelledMNIST :: V.Vector (Int, Mean # Replicated n Bernoulli) -> Random s (Int,[Bool])
sampleRandomLabelledMNIST ldgs = do
    (l,mz) <- randomElement' ldgs
    dg <- standardGenerate mz
    return (l,dg)

sampleRandomMNIST :: V.Vector (Mean # Replicated n Bernoulli) -> Random s [Bool]
sampleRandomMNIST dgs = standardGenerate =<< randomElement' dgs

-}

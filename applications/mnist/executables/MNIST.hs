{-# LANGUAGE TypeOperators,DataKinds #-}

module MNIST where

import Goal.Core

import qualified Goal.Core.Vector.Boxed as B
import qualified Data.Vector.Generic as G

import Data.IDX

import Paths_mnist



--- Harmoniums ---


-- Initialization --

type Height = 28
type Width = 28
type Length = Height * Width
--type Digit m = Replicated (Height * Width) m


---  ---

mnstdr,trnlblfl,trnimgfl,tstlblfl,tstimgfl :: String
mnstdr = "mnist"
trnlblfl = "train-labels-idx1-ubyte"
trnimgfl = "train-images-idx3-ubyte"
tstlblfl = "t10k-labels-idx1-ubyte"
tstimgfl = "t10k-images-idx3-ubyte"

-- IO --

mnistData :: String -> String -> IO [(B.Vector Length Double, Int)]
{-# INLINE mnistData #-}
mnistData lblfl imgfl = do

    lblpth <- getDataFileName lblfl
    imgpth <- getDataFileName imgfl
    mlbls <- decodeIDXLabelsFile lblpth
    mimgs <- decodeIDXFile imgpth

    let (lbls,dgs) = unzip . fromJust $ labeledIntData (fromJust mlbls) (fromJust mimgs)
        dgs' = fmap ((/255) . fromIntegral) . fromJust . B.toSized . G.convert <$> dgs

    return $ zip dgs' lbls

mnistTrainingData :: IO [(B.Vector Length Double, Int)]
{-# INLINE mnistTrainingData #-}
mnistTrainingData = mnistData trnlblfl trnimgfl

mnistTestData :: IO [(B.Vector Length Double, Int)]
{-# INLINE mnistTestData #-}
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

sampleRandomLabelled :: V.Vector (Int, Mean # Replicated n Bernoulli) -> Random s (Int,[Bool])
sampleRandomLabelled ldgs = do
    (l,mz) <- randomElement' ldgs
    dg <- standardGenerate mz
    return (l,dg)

sampleRandom :: V.Vector (Mean # Replicated n Bernoulli) -> Random s [Bool]
sampleRandom dgs = standardGenerate =<< randomElement' dgs

-}

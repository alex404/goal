{-# LANGUAGE ScopedTypeVariables #-}

--- Imports ---


-- Goal --

import NeuralData

import Goal.Core

-- Other --

import qualified Data.List as L


--- Experiments ---


amirExperiment :: FilePath
amirExperiment = "amir-acute"

expmnt :: Experiment
expmnt = Experiment prjnm amirExperiment

dsts :: [String]
dsts = ["bias2to1_acute", "bias2to1_acute2", "bias4to1_acute"]

--- IO ---


parseCSV :: String -> IO ()
parseCSV sbexpnm = do
    ((stms:rws) :: [[Int]]) <- goalImport expmnt sbexpnm
    let cls = L.transpose rws
        k :: Int
        k = length rws
        zxs :: [([Int], Double)]
        zxs = filter ((< 180) . snd) . zip cls $ fromIntegral <$> stms
        (zs,xs) = unzip zxs
    goalWriteDataset expmnt sbexpnm $ show (k,zip zs [pi*x/90 | x <- xs ])

main :: IO ()
main = do

    mapM_ parseCSV dsts
    goalWriteDatasetsCSV expmnt dsts


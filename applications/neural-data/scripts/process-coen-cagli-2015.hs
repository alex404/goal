{-# LANGUAGE ScopedTypeVariables #-}

--- Imports ---


-- Goal --

import NeuralData

import Goal.Core


--- Experiments ---


coenCagli :: String
coenCagli = "coen-cagli-2015"

expmnt :: Experiment
expmnt = Experiment prjnm coenCagli

dsts :: [String]
dsts = ("session" ++) . show <$> [1..10 :: Int]

--- IO ---


parseCSV :: String -> IO ()
parseCSV sbexpnm = do
    Just (rws :: [[Double]]) <- goalImport expmnt sbexpnm
    let stms = [ pi*(x+90)/90 | x <- head <$> rws ]
        rspns = map round . tail <$> rws
        k :: Int
        k = length $ head rspns
        zxs :: [([Int], Double)]
        zxs = zip rspns stms
    putStrLn $ concat ["Parsing ", sbexpnm, "\n", "Number of Neurons: ", show k, "\n", "Number of Trials: ", show $ length rspns]
    goalWriteDataset expmnt sbexpnm $ show (k,zxs)

main :: IO ()
main = do

    mapM_ parseCSV dsts
    goalWriteDatasetsCSV expmnt dsts


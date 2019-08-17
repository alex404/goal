{-# LANGUAGE ScopedTypeVariables #-}

--- Imports ---


-- Goal --

import NeuralData

import Goal.Core


--- Experiments ---


coenCagli :: String
coenCagli = "coen-cagli-2015"

dsts :: [String]
dsts = ("session" ++) . show <$> [1..10 :: Int]

--- IO ---


parseCSV :: String -> IO ()
parseCSV dst = do
    Right (rws :: [[Double]]) <- goalImport $ loadPath (coenCagli ++ "/import") dst
    let stms = [ pi*(x+90)/90 | x <- head <$> rws ]
        rspns = map round . tail <$> rws
        k :: Int
        k = length $ head rspns
        zxs :: [([Int], Double)]
        zxs = zip rspns stms
    putStrLn $ concat [ "Parsing ", dst, "\n"
                      , "Number of Neurons: ", show k, "\n"
                      , "Number of Trials: ", show $ length rspns ]
    let ldpth = loadPath coenCagli dst
    createDirectoryIfMissing True ldpth
    writeFile (ldpth ++ "/dataset.dat") $ show (k,zxs)

main :: IO ()
main = mapM_ parseCSV dsts


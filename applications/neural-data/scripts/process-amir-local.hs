--- Imports ---


-- Goal --

import NeuralData

import Goal.Core

-- Other --

import qualified Data.Vector as V
import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV

import Data.List


--- Experiments ---


coenCagliPath :: FilePath
coenCagliPath = "projects/coen-cagli-2015"


--- IO ---


parseCSV :: String -> [([Int], Int)]
parseCSV csvstr =
    let (stms:rspns) = [ read $ concat ["[", ln, "]"] | ln <- lines csvstr ]
     in zip (transpose rspns) stms

main = return ()
--poolData :: [[(V.Vector Int, Int)]] -> [(V.Vector Int, Int)]
--poolData = foldr1 (zipWith zipper) . map (sortOn snd)
--    where zipper(z1,x1) (z2,x2)
--            | x1 == x2 = (z1 V.++ z2,x1)
--            | otherwise = error "mismatched stimuli"
--
--sessionNames :: [String]
--sessionNames = [ "session" ++ show k | k <- [1 :: Int .. 10]]
--
--sessions :: [String]
--sessions = [ coenCagliPath ++ "/raw-data/" ++ ssn ++ ".csv" | ssn <- sessionNames ]
--
--main :: IO ()
--main = do
--
--    createDirectoryIfMissing True (coenCagliPath ++ "/data")
--
--    zxss <- mapM (fmap parseCSV . readFile) sessions
--    let plzxs = poolData zxss
--
--    sequence_ $ do
--        (ssnnm,zxs) <- zip sessionNames zxss
--        let zxfl = coenCagliPath ++ "/data/" ++ ssnnm ++ ".dat"
--        return $ do
--            putStrLn $ "Coen-Cagli "  ++ ssnnm ++ " Number of Neurons:"
--            print . V.length . fst $ head zxs
--            writeFile zxfl . show $ zxs
--
--    putStrLn "Coen-Cagli Total Number of Neurons:"
--    print . V.length . fst $ head plzxs
--    writeFile (coenCagliPath ++ "/data/pooled.dat") $ show plzxs
--
--    let dsts = Dataset <$> "pooled" : sessionNames
--
--    BS.writeFile (coenCagliPath ++ "/datasets.csv") $ CSV.encodeDefaultOrderedByName dsts

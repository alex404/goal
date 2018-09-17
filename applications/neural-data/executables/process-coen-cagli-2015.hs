--- Imports ---


-- Goal --

import Goal.Core

-- Other --

import qualified Data.Vector as V

import Data.List
import Paths_neural_data


--- Experiments ---

ccprj :: String
ccprj = "coen-cagli-2015"

--- IO ---


parseCSV :: String -> [(V.Vector Int, Int)]
parseCSV csvstr =
    let lns = lines csvstr
     in do
         ln <- lns
         let zx0 = V.fromList . read $ "[" ++ ln ++ "]"
         return (V.tail zx0, V.head zx0)

poolData :: [[(V.Vector Int, Int)]] -> [(V.Vector Int, Int)]
poolData = foldr1 (zipWith zipper) . map (sortOn snd)
    where zipper(z1,x1) (z2,x2)
            | x1 == x2 = (z1 V.++ z2,x1)
            | otherwise = error "mismatched stimuli"

sessions0 :: [String]
sessions0 = [ "session" ++ show k | k <- [1 :: Int .. 10]]

sessions :: [String]
sessions = (++ ".csv") <$> sessions0

main :: IO ()
main = do

    csvpths <- mapM (getDataFileName . ("coen-cagli-2015/" ++)) sessions

    zxss <- mapM (fmap parseCSV . readFile) csvpths
    let plzxs = poolData zxss

    sequence_ $ do
        (ssn,zxs) <- zip sessions0 zxss
        let zxttl = ssn ++ ".dat"
        return $ do
            putStrLn $ "Coen-Cagli "  ++ ssn ++ " Number of Neurons:"
            print . V.length . fst $ head zxs
            goalWriteFile (ccprj ++ "/data") zxttl . show $ zxs

    putStrLn "Coen-Cagli Total Number of Neurons:"
    print . V.length . fst $ head plzxs
    goalWriteFile (ccprj ++ "/data") "pooled.dat" $ show plzxs

    let dst = Dataset <$> "pooled" : sessions0

    goalWriteCSV ccprj "datasets" dst

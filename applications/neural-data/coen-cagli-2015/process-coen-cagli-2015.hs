--- Imports ---


-- Goal --

import Goal.Core

-- Other --

import qualified Data.Vector as V

import Data.List


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

sessionNames :: [String]
sessionNames = [ "session" ++ show k | k <- [1 :: Int .. 10]]

sessions :: [String]
sessions = [ "raw-data/" ++ ssn ++ ".csv" | ssn <- sessionNames ]

main :: IO ()
main = do

    zxss <- mapM (fmap parseCSV . readFile) sessions
    let plzxs = poolData zxss

    sequence_ $ do
        (ssnnm,zxs) <- zip sessionNames zxss
        let zxfl = "data/" ++ ssnnm ++ ".dat"
        return $ do
            putStrLn $ "Coen-Cagli "  ++ ssnnm ++ " Number of Neurons:"
            print . V.length . fst $ head zxs
            writeFile zxfl . show $ zxs

    putStrLn "Coen-Cagli Total Number of Neurons:"
    print . V.length . fst $ head plzxs
    writeFile "data/pooled.dat" $ show plzxs

    let dst = Dataset <$> "pooled" : sessionNames

    goalWriteCSV ccprj "datasets" dst

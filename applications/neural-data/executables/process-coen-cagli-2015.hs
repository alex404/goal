--- Imports ---


-- Goal --

import Goal.Core

-- Other --

import qualified Data.Map as M
import qualified Data.Vector as V
import qualified System.Directory as D

import Data.List
import Data.Char

import System.IO


--- Experiments ---


--- IO ---

newtype DummyVector v = Vector v deriving Show

hReadFile :: String -> IO String
hReadFile pth = do
    hnd <- openFile pth ReadMode
    str <- hGetContents hnd
    str `deepseq` hClose hnd
    return str

parseCSV :: String -> [V.Vector Int]
parseCSV csv0 =
    let lns = lines csv0
     in do
         str <- lns
         let ns = read $ "[" ++ str ++ "]"
         return $ V.fromList ns

parseFilename :: String -> (Int,Int)
parseFilename flnm =
    let (flstr,flnm') = span isDigit $ dropWhile isAlpha flnm
     in (read flstr, read $ dropWhile isAlpha flnm')

groupData :: [(Int,Int)] -> [[V.Vector Int]] -> [[(V.Vector Int, Int)]]
groupData grpxs zss =
    let grpzxs = [ (grp,zip zs $ repeat x) | ((grp,x),zs) <- zip grpxs zss ]
     in M.elems $ M.fromListWith (++) grpzxs

showData :: [(V.Vector Int, Int)] -> String
showData zxs = show [[(Vector z,x) | (z,x) <- zxs]]

poolData :: [[(V.Vector Int, Int)]] -> [(V.Vector Int, Int)]
poolData = foldr1 (zipWith zipper) . map (sortOn snd)
    where zipper(z1,x1) (z2,x2)
            | x1 == x2 = (z1 V.++ z2,x1)
            | otherwise = error "mismatched stimuli"


main :: IO ()
main = do

    gdpth <- goalRawDataDirectory
    let ccstr = "coen-cagli-2015"
        csvdr = gdpth ++ "/" ++ ccstr ++ "/" ++ "csvs"

    flnms <- D.listDirectory csvdr

    zss <- mapM (fmap parseCSV . hReadFile) [ csvdr ++ "/" ++ flnm | flnm <- flnms ]
    let grpxs = parseFilename <$> flnms
        grpzxss = groupData grpxs zss
        grpzxs = poolData grpzxss

    putStrLn "Coen-Cagli Total Number of Neurons:"
    print . V.length . fst $ head grpzxs
    goalWriteFile (ccstr ++ "/data") "zxs-pooled" $ showData grpzxs

    sequence_ $ do
        (grp,zxs) <- zip [(1 :: Int)..10] grpzxss
        let zxttl = "zxs" ++ show grp
        return $ do
            putStrLn $ "Coen-Cagli Group "  ++ show grp ++ " Number of Neurons:"
            print . V.length . fst $ head zxs
            goalWriteFile (ccstr ++ "/data") zxttl . showData $ zxs

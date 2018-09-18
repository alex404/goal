{-# LANGUAGE DeriveGeneric #-}

-- | This module exports a set of generic numerical and list manipulation functions, as well as a
-- set of Goal-specific functions for file and directory manipulation. These functions use the XDG
-- directory specification to save files in appropriate directories.
module Goal.Core.Util
    ( -- * List Manipulation
      takeEvery
    , breakEvery
    -- * Low-Level
    , traceGiven
    -- * Numeric
    , roundSD
    , toPi
    , integrate
    , logistic
    , logit
    , square
    -- ** List Numerics
    , average
    , range
    , discretizeFunction
    -- * Goal IO
    -- ** Project Management
    , goalProjectDirectory
    , goalCreateProject
    , goalRemoveProject
    -- ** File Management
    , goalFilePath
    , goalDoesFileExist
    , goalWriteFile
    , goalReadFile
    , goalDeleteFile
    -- *** CSV Management
    , goalReadCSV
    , goalWriteCSV
    -- ** Plot Rendering
    , goalRenderableToPNG
    , goalRenderableToSVG
    , goalRenderableToPDF
    -- ** Csvs
    , Collection (Collection)
    , Dataset (Dataset)
    , getCollections
    , getDatasets
    ) where


--- Imports ---


import Control.Monad
import Debug.Trace
import System.Directory
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo
import GHC.Generics

-- Qualified --

import qualified Numeric.Integration.TanhSinh as I
import qualified Data.Csv as CSV
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V

--- General Functions ---


-- | Takes every nth element, starting with the head of the list.
takeEvery :: Int -> [x] -> [x]
{-# INLINE takeEvery #-}
takeEvery m = map snd . filter (\(x,_) -> mod x m == 0) . zip [0..]

-- | Break the list up into lists of length n.
breakEvery :: Int -> [x] -> [[x]]
{-# INLINE breakEvery #-}
breakEvery _ [] = []
breakEvery n xs = take n xs : breakEvery n (drop n xs)

-- | Runs traceShow on the given element.
traceGiven :: Show a => a -> a
traceGiven a = traceShow a a


--- Data management ---


newtype Collection = Collection { collection :: String } deriving (Read,Show,Generic)

newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance CSV.FromNamedRecord Collection
instance CSV.ToNamedRecord Collection
instance CSV.DefaultOrdered Collection

instance CSV.FromNamedRecord Dataset
instance CSV.ToNamedRecord Dataset
instance CSV.DefaultOrdered Dataset

getCollections :: IO [Collection]
getCollections = do

    bstrm <- BS.readFile "collections.csv"
    let Right (_,as) = CSV.decodeByName bstrm

    return $ V.toList as

getDatasets :: Collection -> IO [Dataset]
getDatasets (Collection clc) = do

    bstrm <- BS.readFile $ clc ++ "/datasets.csv"
    let Right (_,as) = CSV.decodeByName bstrm

    return $ V.toList as



--- Numeric ---

-- | Numerically integrates a 1-d function over an interval.
integrate
    :: Double -- ^ Error Tolerance
    -> (Double -> Double) -- ^ Function
    -> Double -- ^ Interval beginning
    -> Double -- ^ Interval end
    -> Double -- ^ Integral
integrate err f mn mx =
    I.result . I.absolute err $ I.parTrap f mn mx

-- | Rounds the number to the specified significant digit.
roundSD :: RealFloat x => Int -> x -> x
{-# INLINE roundSD #-}
roundSD n x =
    let n' :: Int
        n' = round $ 10^n * x
     in fromIntegral n'/10^n

-- | Modulo's a real value to be in [0,2pi]
toPi :: RealFloat x => x -> x
{-# INLINE toPi #-}
toPi x =
    let xpi = x / (2*pi)
        f = xpi - fromIntegral (floor xpi :: Int)
     in 2 * pi * f

-- | A standard sigmoid function.
logistic :: Floating x => x -> x
{-# INLINE logistic #-}
logistic x = 1 / (1 + exp(negate x))

-- | The inverse of the logistic.
logit :: Floating x => x -> x
{-# INLINE logit #-}
logit x = log $ x / (1 - x)

-- | The square of a number (for avoiding endless default values).
square :: Floating x => x -> x
{-# INLINE square #-}
square x = x^(2::Int)

-- Lists --

-- | Average value of a 'Traversable' of 'Fractional's.
average :: (Foldable f, Fractional x) => f x -> x
{-# INLINE average #-}
average = uncurry (/) . foldr (\e (s,c) -> (e+s,c+1)) (0,0)

-- | Returns n numbers which uniformly partitions the interval [mn,mx].
range
    :: RealFloat x => x -> x -> Int -> [x]
{-# INLINE range #-}
range _ _ 0 = []
range mn mx 1 = [(mn + mx) / 2]
range mn mx n =
    [ x * mx + (1 - x) * mn | x <- (/ (fromIntegral n - 1)) . fromIntegral <$> [0 .. n-1] ]

-- | Takes range information in the form of a minimum, maximum, and sample count,
-- a function to sample, and returns a list of pairs (x,f(x)) over the specified
-- range.
discretizeFunction :: Double -> Double -> Int -> (Double -> Double) -> [(Double,Double)]
{-# INLINE discretizeFunction #-}
discretizeFunction mn mx n f =
    let rng = range mn mx n
    in zip rng $ f <$> rng


--- Goal directory management ---


-- | Returns the xdg-based directory where projects are stored in Goal.
goalProjectDirectory :: IO FilePath
goalProjectDirectory = getXdgDirectory XdgData "goal/projects"

-- | Returns the path to a file given its name and the name of a Goal project.
goalFilePath
    :: String -- ^ Goal Project
    -> String -- ^ File name
    -> IO FilePath -- ^ Path
{-# INLINE goalFilePath #-}
goalFilePath sbdr flnm = do
    gldr <- goalProjectDirectory
    return $ gldr ++ "/" ++ sbdr ++ "/" ++ flnm

-- | Writes a file with the given filename.
goalWriteFile
    :: String -- ^ Goal Project
    -> String -- ^ File name
    -> String -- ^ File contents
    -> IO ()
{-# INLINE goalWriteFile #-}
goalWriteFile sbdr flnm txt = do
    sbpth <- goalCreateProject sbdr
    let fpth' =  sbpth ++ "/" ++ flnm
    writeFile fpth' txt

-- | Read a file in the given Goal project with the given file name.
goalReadFile
    :: String -- ^ Goal project name
    -> String -- ^ File name
    -> IO String -- ^ File Contents
{-# INLINE goalReadFile #-}
goalReadFile sbdr flnm = do
    fpth <- goalFilePath sbdr flnm
    readFile fpth



--- CSV ---


goalReadCSV ::  CSV.FromNamedRecord a => String -> String -> IO [a]
goalReadCSV prj flnm = do

    csvpth <- goalFilePath prj (flnm ++ ".csv")
    bstrm <- BS.readFile csvpth

    let Right (_,as) = CSV.decodeByName bstrm

    return $ V.toList as

goalWriteCSV ::  (CSV.DefaultOrdered a, CSV.ToNamedRecord a) => String -> String -> [a] -> IO ()
goalWriteCSV prj flnm as = do
    sbpth <- goalCreateProject prj
    let fpth = sbpth ++ "/" ++ flnm ++ ".csv"
    BS.writeFile fpth $ CSV.encodeDefaultOrderedByName as

-- | Checks the existence of a file in the given project and with the given name.
goalDoesFileExist :: String -> String -> IO Bool
goalDoesFileExist sbdr flnm = goalFilePath sbdr flnm >>= doesFileExist

-- | Deletes a file in the given project and with the given name.
goalDeleteFile :: String -> String -> IO ()
goalDeleteFile sbdr flnm = goalFilePath sbdr flnm >>= removeFile

-- | Removes a project with the given name.
goalRemoveProject :: String -> IO ()
goalRemoveProject sbdr = do
    gldr <- goalProjectDirectory
    let sbpth = gldr ++ "/" ++ sbdr
    removeDirectoryRecursive sbpth

-- | Creates a project directory with the given name and returns its absolute path.
goalCreateProject :: String -> IO FilePath
goalCreateProject sbdr = do
    gldr <- goalProjectDirectory
    let sbpth = gldr ++ "/" ++ sbdr
    createDirectoryIfMissing True sbpth
    return sbpth

-- | Given the project name and file name (without the extension), saves
-- the given renderable in the project directory as a PNG.
goalRenderableToPNG
    :: String -- ^ Project name
    -> String -- ^ File name
    -> Int -- ^ Pixel width
    -> Int -- ^ Pixel height
    -> Renderable a -- ^ The Renderable
    -> IO ()
goalRenderableToPNG sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateProject sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".png"
    void $ renderableToFile (FileOptions (xn,yn) PNG) fpth' rnbl

-- | Given the project name and file name (without the extension), saves
-- the given renderable in the project directory as a SVG.
goalRenderableToSVG
    :: String -- ^ Project name
    -> String -- ^ File name
    -> Int -- ^ Pixel width
    -> Int -- ^ Pixel height
    -> Renderable a -- ^ The Renderable
    -> IO ()
goalRenderableToSVG sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateProject sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".svg"
    void $ renderableToFile (FileOptions (xn,yn) SVG) fpth' rnbl

-- | Given the project name and file name (without the extension), saves
-- the given renderable in the project directory as a PDF.
goalRenderableToPDF
    :: String -- ^ Project name
    -> String -- ^ File name
    -> Int -- ^ Pixel width
    -> Int -- ^ Pixel height
    -> Renderable a -- ^ The Renderable
    -> IO ()
goalRenderableToPDF sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateProject sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".pdf"
    void $ renderableToFile (FileOptions (xn,yn) PDF) fpth' rnbl

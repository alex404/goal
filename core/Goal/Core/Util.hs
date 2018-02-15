-- | This module exports a set of generic numerical and list manipulation functions, as well as a
-- set of Goal-specific functions for file and directory manipulation. These functions use the XDG
-- directory specification to save files in appropriate directories.
module Goal.Core.Util
    ( -- * Lists
      takeEvery
    , breakEvery
    -- * Low-Level
    , traceGiven
    -- * Numeric
    , roundSD
    , toPi
    , integrate
    -- ** Functions
    , logistic
    , logit
    -- ** Lists
    , average
    , range
    , discretizeFunction
    -- * Goal IO
    -- ** File Management
    , goalDoesFileExist
    , goalWriteFile
    , goalReadFile
    , goalDeleteFile
    , goalRenderableToPNG
    , goalRenderableToSVG
    , goalRenderableToPDF
    -- ** Directory Management
    , goalProjectDirectory
    , goalDatasetDirectory
    , goalCreateSubdirectory
    , goalRemoveSubdirectory
    , goalProjectLocation
    , goalDatasetLocation
    ) where


--- Imports ---


import Control.Monad
import Debug.Trace
import System.Directory
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo

-- Qualified --

import qualified Numeric.Integration.TanhSinh as I

--- General Functions ---


-- | Takes every nth element, starting with the head of the list.
takeEvery :: Int -> [x] -> [x]
takeEvery m = map snd . filter (\(x,_) -> mod x m == 0) . zip [0..]

-- | Break the list up into lists of length n.
breakEvery :: Int -> [x] -> [[x]]
breakEvery _ [] = []
breakEvery n xs = take n xs : breakEvery n (drop n xs)

-- | Runs traceShow on the given element.
traceGiven :: Show a => a -> a
traceGiven a = traceShow a a


--- Numeric ---

integrate :: Double -> (Double -> Double) -> Double -> Double -> Double
integrate err f mn mx =
    I.result . I.absolute err $ I.parTrap f mn mx

-- | Rounds the number to the specified significant digit.
roundSD :: RealFloat x => Int -> x -> x
{-# INLINE roundSD #-}
roundSD n x =
    let n' :: Int
        n' = round $ 10^n * x
     in fromIntegral n'/10^n

-- | Modulo's a real value to be in [-pi,pi]
toPi :: RealFloat x => x -> x
{-# INLINE toPi #-}
toPi x =
    let xpi = x / pi
        n :: Int
        n = floor xpi
        f = xpi - fromIntegral n
    in if even n then pi * f else -(pi * (1 - f))

-- | A standard sigmoid function.
logistic :: Floating x => x -> x
{-# INLINE logistic #-}
logistic x = 1 / (1 + exp(negate x))

-- | The inverse of the logistic.
logit :: Floating x => x -> x
{-# INLINE logit #-}
logit x = log $ x / (1 - x)

-- Lists --

-- | Average value of a 'Traversable' of 'Fractional's.
average :: (Foldable f, Fractional x) => f x -> x
{-# INLINE average #-}
average = uncurry (/) . foldr (\e (s,c) -> (e+s,c+1)) (0,0)

-- | Returns n  numbers from mn to mx.
range :: RealFloat x => x -> x -> Int -> [x]
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


goalProjectDirectory :: IO FilePath
goalProjectDirectory = getXdgDirectory XdgData "goal/projects"

goalDatasetDirectory :: IO FilePath
goalDatasetDirectory = getXdgDirectory XdgData "goal/datasets"

goalProjectLocation :: String -> String -> IO FilePath
goalProjectLocation sbdr flnm = do
    gldr <- goalProjectDirectory
    return $ gldr ++ "/" ++ sbdr ++ "/" ++ flnm

goalDatasetLocation :: String -> String -> IO FilePath
goalDatasetLocation sbdr flnm = do
    gldr <- goalDatasetDirectory
    return $ gldr ++ "/" ++ sbdr ++ "/" ++ flnm


-- | Writes a file to the goal data directory. The first two arguments are the
-- subdirectory and filename, and the third is the string to write.
goalWriteFile :: String -> String -> String -> IO ()
goalWriteFile sbdr flnm txt = do
    sbpth <- goalCreateSubdirectory sbdr
    let fpth' =  sbpth ++ "/" ++ flnm
    writeFile fpth' txt

-- | Read a file from the goal data directory. The first two arguments are the
-- subdirectory and filename.
goalReadFile :: String -> String -> IO String
goalReadFile sbdr flnm = goalProjectLocation sbdr flnm >>= readFile

goalDoesFileExist :: String -> String -> IO Bool
goalDoesFileExist sbdr flnm = goalProjectLocation sbdr flnm >>= doesFileExist

-- | Read a file from the goal data directory. The first two arguments are the
-- subdirectory and filename.
goalDeleteFile :: String -> String -> IO ()
goalDeleteFile sbdr flnm = goalProjectLocation sbdr flnm >>= removeFile

goalRemoveSubdirectory :: String -> IO ()
goalRemoveSubdirectory sbdr = do
    gldr <- goalProjectDirectory
    let sbpth = gldr ++ "/" ++ sbdr
    removeDirectoryRecursive sbpth

goalCreateSubdirectory :: String -> IO String
goalCreateSubdirectory sbdr = do
    gldr <- goalProjectDirectory
    let sbpth = gldr ++ "/" ++ sbdr
    createDirectoryIfMissing True sbpth
    return sbpth

-- | Given the desired subdirectory and filename (without the extension), saves
-- the given renderable in the goal data directory as a PNG.
goalRenderableToPNG :: String -> String -> Int -> Int -> Renderable a -> IO ()
goalRenderableToPNG sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateSubdirectory sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".png"
    void $ renderableToFile (FileOptions (xn,yn) PNG) fpth' rnbl

-- | Given the desired subdirectory and filename (without the extension), saves
-- the given renderable in the goal data directory as a SVG.
goalRenderableToSVG :: String -> String -> Int -> Int -> Renderable a -> IO ()
goalRenderableToSVG sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateSubdirectory sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".svg"
    void $ renderableToFile (FileOptions (xn,yn) SVG) fpth' rnbl

-- | Given the desired subdirectory and filename (without the extension), saves
-- the given renderable in the goal data directory as a PDF.
goalRenderableToPDF :: String -> String -> Int -> Int -> Renderable a -> IO ()
goalRenderableToPDF sbdr flnm xn yn rnbl = do
    sbpth <- goalCreateSubdirectory sbdr
    let fpth' =  sbpth ++ "/" ++ flnm ++ ".pdf"
    void $ renderableToFile (FileOptions (xn,yn) PDF) fpth' rnbl

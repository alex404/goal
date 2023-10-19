-- | Helper functions for exporting data as JSON and running Python scripts.
module Goal.Core.Project (
    exportJSON,
    runPythonScriptWithArg,
    toJSON,
    (.=),
    findCabalProjectRoot,
    benchFilePath,
    resultsFilePath,
) where

--- Imports ---

import Data.Aeson (Key, ToJSON, Value, encode, object, (.=))
import Data.ByteString.Lazy qualified as BL
import System.Directory (createDirectoryIfMissing, doesFileExist, getCurrentDirectory)
import System.FilePath (takeDirectory, (</>))
import System.Process (callCommand)

--- Globals ---

dataDir :: FilePath
dataDir = "data-dir"

benchDir :: FilePath
benchDir = "bench-reports"

resultsDir :: FilePath
resultsDir = "results"

plotDir :: FilePath
plotDir = "examples" </> "plots"

--- Functions ---

-- Export data as JSON
exportJSON :: (ToJSON a) => FilePath -> a -> IO ()
exportJSON pth dat = BL.writeFile pth (encode dat)

-- Give JSON creation a less generic name
toJSON :: [(Key, Value)] -> Value
toJSON = object

-- Run a Python script
runPythonScriptWithArg :: FilePath -> FilePath -> IO ()
runPythonScriptWithArg pynm jsnnm = do
    cblrt <- findCabalProjectRoot
    putStrLn $ "Running Python Script: " ++ pynm
    let pypth = cblrt </> plotDir </> pynm
    callCommand $ "python " ++ pypth ++ " " ++ jsnnm

-- Find the project root directory
findCabalProjectRoot :: IO FilePath
findCabalProjectRoot = do
    currentDir <- getCurrentDirectory
    findCabalProjectRoot' currentDir

findCabalProjectRoot' :: FilePath -> IO FilePath
findCabalProjectRoot' dir = do
    let projectFile = dir </> "cabal.project"
    exists <- doesFileExist projectFile
    if exists
        then return dir
        else
            let parentDir = takeDirectory dir
             in if parentDir == dir -- We are at the root
                    then error "cabal.project file not found in any parent directory."
                    else findCabalProjectRoot' parentDir

-- Function to create and return the path for bench reports
benchFilePath :: FilePath -> IO FilePath
benchFilePath flnm = do
    cblrt <- findCabalProjectRoot
    let benchPath = cblrt </> dataDir </> benchDir
    createDirectoryIfMissing True benchPath
    return $ benchPath </> flnm

-- Function to create and return the path for results
resultsFilePath :: FilePath -> IO FilePath
resultsFilePath flnm = do
    cblrt <- findCabalProjectRoot
    let resultsPath = cblrt </> dataDir </> resultsDir
    createDirectoryIfMissing True resultsPath
    return $ resultsPath </> flnm

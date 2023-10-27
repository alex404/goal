{- | Helper functions for processing data in the workspace directory, e.g. exporting
data as JSON.
-}
module Goal.Core.Workspace (
    exportJSON,
    toJSON,
    (.=),
    findWorkspace,
    benchFilePath,
    resultsFilePath,
) where

--- Imports ---

import Data.Aeson (Key, ToJSON, Value, encode, object, (.=))
import Data.ByteString.Lazy qualified as BL
import System.Directory (createDirectoryIfMissing, doesFileExist, getCurrentDirectory)
import System.FilePath (takeDirectory, (</>))

--- Globals ---

dataDir :: FilePath
dataDir = "workspace"

benchDir :: FilePath
benchDir = "bench-reports"

resultsDir :: FilePath
resultsDir = "results"

--- Functions ---

-- Export data as JSON
exportJSON :: (ToJSON a) => FilePath -> a -> IO ()
exportJSON pth dat = BL.writeFile pth (encode dat)

-- Give JSON creation a less generic name
toJSON :: [(Key, Value)] -> Value
toJSON = object

-- Find the project root directory
findWorkspace :: IO FilePath
findWorkspace = do
    currentDir <- getCurrentDirectory
    prjctrt <- findCabalProjectRoot currentDir
    return $ prjctrt </> dataDir

findCabalProjectRoot :: FilePath -> IO FilePath
findCabalProjectRoot dir = do
    let projectFile = dir </> "cabal.project"
    exists <- doesFileExist projectFile
    if exists
        then return dir
        else
            let parentDir = takeDirectory dir
             in if parentDir == dir -- We are at the root
                    then error "cabal.project file not found in any parent directory."
                    else findCabalProjectRoot parentDir

-- Function to create and return the path for bench reports
benchFilePath :: FilePath -> IO FilePath
benchFilePath flnm = do
    wrkspc <- findWorkspace
    let benchPath = wrkspc </> benchDir
    createDirectoryIfMissing True benchPath
    return $ benchPath </> flnm

-- Function to create and return the path for results
resultsFilePath :: FilePath -> IO FilePath
resultsFilePath flnm = do
    wrkspc <- findWorkspace
    let resultsPath = wrkspc </> resultsDir
    createDirectoryIfMissing True resultsPath
    return $ resultsPath </> flnm

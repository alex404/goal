{- | Helper functions for processing data in the workspace directory, e.g. exporting
data as JSON.
-}
module Goal.Core.Workspace (
    exportJSON,
    toJSON,
    importJSON,
    (.=),
    findWorkspace,
    benchFilePath,
    resultsFilePath,
    dataFilePath,
) where

--- Imports ---

import Data.Aeson (FromJSON, Key, ToJSON, Value, eitherDecode, encode, object, (.=))
import Data.ByteString.Lazy qualified as BL
import System.Directory (createDirectoryIfMissing, doesFileExist, getCurrentDirectory)
import System.FilePath (takeDirectory, (</>))

--- Globals ---

workspaceDir :: FilePath
workspaceDir = "workspace"

benchDir :: FilePath
benchDir = "bench-reports"

resultsDir :: FilePath
resultsDir = "results"

dataDir :: FilePath
dataDir = "data-files"

--- Functions ---

-- Export data as JSON
exportJSON :: (ToJSON a) => FilePath -> a -> IO ()
exportJSON pth dat = BL.writeFile pth (encode dat)

importJSON :: (FromJSON a) => FilePath -> IO a
importJSON pth = do
    dat <- BL.readFile pth
    return $ either error id (eitherDecode dat)

-- Give JSON creation a less generic name
toJSON :: [(Key, Value)] -> Value
toJSON = object

-- Find the project root directory
findWorkspace :: IO FilePath
findWorkspace = do
    currentDir <- getCurrentDirectory
    prjctrt <- findCabalProjectRoot currentDir
    return $ prjctrt </> workspaceDir

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

-- Function to create and return the path for results
subdirFilePath :: FilePath -> FilePath -> IO FilePath
subdirFilePath dr flnm = do
    wrkspc <- findWorkspace
    let pth = wrkspc </> dr
    createDirectoryIfMissing True pth
    return $ pth </> flnm

-- Function to create and return the path for bench reports
benchFilePath :: FilePath -> IO FilePath
benchFilePath = subdirFilePath benchDir

-- Function to create and return the path for results
resultsFilePath :: FilePath -> IO FilePath
resultsFilePath = subdirFilePath resultsDir

-- Function to create and return the path for results
dataFilePath :: FilePath -> IO FilePath
dataFilePath = subdirFilePath dataDir

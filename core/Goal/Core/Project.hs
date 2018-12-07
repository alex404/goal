{-# LANGUAGE DeriveGeneric #-}

-- | A set of Goal-specific functions for file, directory, and csv manipulation.
-- These functions use the XDG directory specification to save files in
-- appropriate directories.
module Goal.Core.Project
    (
    -- * Project Management
      goalProjectDirectory
    , goalRawDataDirectory
    , goalProjectPath
    , goalExperimentPath
    -- * File Management
    , goalReadFile
    , goalWriteFile
    -- * Dataset Management
    , Dataset (Dataset)
    , goalWriteDatasetsCSV
    , goalReadDatasetsCSV
    , goalWriteDataset
    , goalReadDataset
    -- * Analysis
    , goalWriteAnalysis
    , goalWriteNamedAnalysis
    , goalAppendAnalysis
    , goalAppendNamedAnalysis
    -- * Criterion
    , goalCriterionMain
    ) where


--- Imports ---


import System.Directory
import GHC.Generics
import Data.String
import Control.Monad

-- Qualified --

import qualified Data.Vector as V
import qualified Data.ByteString.Lazy as BS
import qualified Data.Csv as CSV
import qualified Criterion.Main as C
import qualified Criterion.Types as C

--- Goal Projects ---


-- | Returns the xdg-based directory where projects are stored in Goal.
goalProjectDirectory :: IO FilePath
goalProjectDirectory = getXdgDirectory XdgData "goal/projects"

-- | Creates a project directory with the given name and returns its absolute path.
goalProjectPath :: String -> IO FilePath
goalProjectPath prnm = do
    prdr <- goalProjectDirectory
    return $ prdr ++ "/" ++ prnm

-- | Creates a project directory with the given name and returns its absolute path.
goalExperimentPath :: String -> String -> IO FilePath
goalExperimentPath prnm expnm = do
    prpth <- goalProjectPath prnm
    return $ prpth ++ "/" ++ expnm

-- | Creates a project directory with the given name and returns its absolute path.
goalRawDataDirectory :: IO FilePath
goalRawDataDirectory =
    getXdgDirectory XdgData "goal/raw-data"

-- | Read a file in the given Goal project with the given file name.
goalReadFile
    :: String -- ^ Goal project name
    -> String -- ^ Goal experiment name
    -> String -- ^ File name
    -> IO String -- ^ File Contents
{-# INLINE goalReadFile #-}
goalReadFile prnm expnm flnm = do
    fpth <- goalExperimentPath prnm expnm
    readFile $ fpth ++ "/" ++ flnm

-- | Writes a file with the given filename.
goalWriteFile
    :: String -- ^ Goal Name
    -> String -- ^ Experiment Name
    -> String -- ^ File name
    -> String -- ^ File contents
    -> IO ()
{-# INLINE goalWriteFile #-}
goalWriteFile prnm expnm flnm txt = do
    pth <- goalExperimentPath prnm expnm
    createDirectoryIfMissing True pth
    let pth' =  pth ++ "/" ++ flnm
    writeFile pth' txt


--- Data management ---


-- | Newtype for working with multiple datasets in a particular project.
newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance CSV.FromNamedRecord Dataset
instance CSV.ToNamedRecord Dataset
instance CSV.DefaultOrdered Dataset

-- | Write the list of datasets.
goalWriteDatasetsCSV :: String -> String -> [Dataset] -> IO ()
goalWriteDatasetsCSV prjnm expnm dsts = do

    pth <- goalExperimentPath prjnm expnm
    createDirectoryIfMissing True pth

    BS.writeFile (pth ++ "/datasets.csv") $ CSV.encodeDefaultOrderedByName dsts

-- | Read the list of datasets (if it exists).
goalReadDatasetsCSV :: String -> String -> IO (Maybe [Dataset])
goalReadDatasetsCSV prjnm expnm = do

    pth <- goalExperimentPath prjnm expnm
    let fpth = pth ++ "/datasets.csv"
    bl <- doesFileExist fpth
    if bl
       then do
           bstrm <- BS.readFile fpth
           let Right (_,as) = CSV.decodeByName bstrm
           return . Just $ V.toList as
       else return Nothing

-- | Run criterion with some defaults for goal. In particular save the results
-- in the benchmarks project.
goalCriterionMain :: String -> [C.Benchmark] -> IO ()
goalCriterionMain expnm bmrks = do

    exppth <- goalExperimentPath "benchmarks" expnm
    createDirectoryIfMissing True exppth

    let rptpth = exppth ++ "/" ++ "report.html"

    C.defaultMainWith (C.defaultConfig { C.reportFile = Just rptpth}) bmrks

-- | Write a single dataset to a file.
goalWriteDataset
    :: String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> Dataset -- ^ Dataset name
    -> String -- ^ File contents
    -> IO ()
goalWriteDataset prjnm expnm (Dataset dstnm) fl = do

    exppth <- goalExperimentPath prjnm expnm

    let pth = exppth ++ "/data"
    createDirectoryIfMissing True pth

    let flpth = pth ++ "/" ++ dstnm ++ ".dat"
    writeFile flpth fl

-- | Read a single dataset from a file.
goalReadDataset
    :: String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> Dataset -- ^ Dataset name
    -> IO String
goalReadDataset prjnm expnm (Dataset dstnm) = do

    exppth <- goalExperimentPath prjnm expnm
    let flpth = exppth ++ "/data/" ++ dstnm ++ ".dat"

    readFile flpth

analysisFilePath :: Bool -> String -> String -> String -> Maybe Dataset -> IO String
analysisFilePath cbl prjnm expnm ananm mdst = do

    exppth <- goalExperimentPath prjnm expnm

    case mdst of
      Just (Dataset dstnm) -> do
          let flpth0 = exppth ++ "/analysis/" ++ ananm
          when cbl $ createDirectoryIfMissing True flpth0
          return $ concat [flpth0,"/",dstnm,".csv"]
      Nothing -> do
          let flpth0 = exppth ++ "/analysis"
          when cbl $ createDirectoryIfMissing True flpth0
          return $ concat [flpth0,"/",ananm,".csv"]

-- | Write the results of an analysis (in the form of a CSV) to the project
-- directory, using the goal organization structure. If there are multiple
-- datasets, analyses are named by dataset in a folder named after the analysis,
-- and otherwise the csv takes the name of the analysis itself.
goalWriteAnalysis
    :: CSV.ToField x
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [[x]] -- ^ CSVs
    -> IO ()
goalWriteAnalysis prjnm expnm ananm mdst csvs = do

    flpth <- analysisFilePath True prjnm expnm ananm mdst

    BS.writeFile flpth $ CSV.encode csvs

-- | Write an analysis CSV based on named records.
goalWriteNamedAnalysis
    :: (CSV.DefaultOrdered csv, CSV.ToNamedRecord csv)
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [csv] -- ^ CSVs
    -> IO ()
goalWriteNamedAnalysis prjnm expnm ananm mdst csvs = do

    flpth <- analysisFilePath True prjnm expnm ananm mdst

    BS.writeFile flpth $ CSV.encodeDefaultOrderedByName csvs

-- | Append an analysis CSV to an existing analysis CSV. Analysis blocks are
-- seperated by two lines, allowing gnuplot to distinguish them as different CSV
-- blocks.
goalAppendAnalysis
    :: CSV.ToField x
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [[x]] -- ^ CSVs
    -> IO ()
goalAppendAnalysis prjnm expnm ananm mdst csvs = do

    flpth <- analysisFilePath False prjnm expnm ananm mdst

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encode csvs

-- | Append an analysis CSV to an existing analysis CSV with named columns.
goalAppendNamedAnalysis
    :: (CSV.DefaultOrdered csv, CSV.ToNamedRecord csv)
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [csv] -- ^ CSVs
    -> IO ()
goalAppendNamedAnalysis prjnm expnm ananm mdst csvs = do

    flpth <- analysisFilePath False prjnm expnm ananm mdst

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encodeDefaultOrderedByName csvs

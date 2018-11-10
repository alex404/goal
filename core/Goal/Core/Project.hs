{-# LANGUAGE DeriveGeneric #-}

-- | This module exports a set of generic numerical and list manipulation functions, as well as a
-- set of Goal-specific functions for file and directory manipulation. These functions use the XDG
-- directory specification to save files in appropriate directories.
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
goalRawDataDirectory = do
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

-- | Checks the existence of a file in the given project and with the given name.
--goalDoesFileExist :: String -> String -> String -> IO Bool
--goalDoesFileExist sbdr expnm flnm = do
--    pth <- goalExperimentPath sbdr expnm
--    doesFileExist $ pth ++ "/" ++ flnm


--- Data management ---


newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance CSV.FromNamedRecord Dataset
instance CSV.ToNamedRecord Dataset
instance CSV.DefaultOrdered Dataset

goalWriteDatasetsCSV :: String -> String -> [Dataset] -> IO ()
goalWriteDatasetsCSV prjnm expnm dsts = do

    pth <- goalExperimentPath prjnm expnm
    createDirectoryIfMissing True pth

    BS.writeFile (pth ++ "/datasets.csv") $ CSV.encodeDefaultOrderedByName dsts


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

goalCriterionMain :: String -> [C.Benchmark] -> IO ()
goalCriterionMain expnm bmrks = do

    exppth <- goalExperimentPath "benchmarks" expnm
    createDirectoryIfMissing True exppth

    let rptpth = exppth ++ "/" ++ "report.html"

    C.defaultMainWith (C.defaultConfig { C.reportFile = Just rptpth}) bmrks

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

goalReadDataset
    :: String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> Dataset -- ^ Dataset name
    -> IO String
goalReadDataset prjnm expnm (Dataset dstnm) = do

    exppth <- goalExperimentPath prjnm expnm
    let flpth = exppth ++ "/data/" ++ dstnm ++ ".dat"

    readFile flpth

goalWriteAnalysis
    :: CSV.ToField x
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [[x]] -- ^ CSVs
    -> IO ()
goalWriteAnalysis prjnm expnm ananm mdst csvs = do

    exppth <- goalExperimentPath prjnm expnm

    flpth <- case mdst of
               Just (Dataset dstnm) -> do
                   let flpth0 = exppth ++ "/analysis/" ++ ananm
                   createDirectoryIfMissing True flpth0
                   return $ concat [flpth0,"/",dstnm,".csv"]
               Nothing -> do
                   let flpth0 = exppth ++ "/analysis"
                   createDirectoryIfMissing True flpth0
                   return $ concat [flpth0,"/",ananm,".csv"]

    BS.writeFile flpth $ CSV.encode csvs

goalWriteNamedAnalysis
    :: (CSV.DefaultOrdered csv, CSV.ToNamedRecord csv)
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [csv] -- ^ CSVs
    -> IO ()
goalWriteNamedAnalysis prjnm expnm ananm mdst csvs = do

    exppth <- goalExperimentPath prjnm expnm

    flpth <- case mdst of
               Just (Dataset dstnm) -> do
                   let flpth0 = exppth ++ "/analysis/" ++ ananm
                   createDirectoryIfMissing True flpth0
                   return $ concat [flpth0,"/",dstnm,".csv"]
               Nothing -> do
                   let flpth0 = exppth ++ "/analysis"
                   createDirectoryIfMissing True flpth0
                   return $ concat [flpth0,"/",ananm,".csv"]

    BS.writeFile flpth $ CSV.encodeDefaultOrderedByName csvs

goalAppendAnalysis
    :: CSV.ToField x
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [[x]] -- ^ CSVs
    -> IO ()
goalAppendAnalysis prjnm expnm ananm mdst csvs = do

    exppth <- goalExperimentPath prjnm expnm

    flpth <- case mdst of
               Just (Dataset dstnm) -> do
                   let flpth0 = exppth ++ "/analysis/" ++ ananm
                   return $ concat [flpth0,"/",dstnm,".csv"]
               Nothing -> do
                   let flpth0 = exppth ++ "/analysis"
                   return $ concat [flpth0,"/",ananm,".csv"]

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encode csvs


goalAppendNamedAnalysis
    :: (CSV.DefaultOrdered csv, CSV.ToNamedRecord csv)
    => String -- ^ Project Name
    -> String -- ^ Experiment Name
    -> String -- ^ Analysis Name
    -> Maybe Dataset -- ^ Maybe dataset name
    -> [csv] -- ^ CSVs
    -> IO ()
goalAppendNamedAnalysis prjnm expnm ananm mdst csvs = do

    exppth <- goalExperimentPath prjnm expnm

    flpth <- case mdst of
               Just (Dataset dstnm) -> do
                   let flpth0 = exppth ++ "/analysis/" ++ ananm
                   return $ concat [flpth0,"/",dstnm,".csv"]
               Nothing -> do
                   let flpth0 = exppth ++ "/analysis"
                   return $ concat [flpth0,"/",ananm,".csv"]

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encodeDefaultOrderedByName csvs



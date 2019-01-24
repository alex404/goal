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
    , Experiment (Experiment, projectName, experimentName)
    , SubExperiment (SubExperiment, analysisName, datasetName)
    , goalWriteDatasetsCSV
    , goalReadDatasetsCSV
    , goalWriteDataset
    , goalReadDataset
    -- * Analysis
    , goalWriteAnalysis
    , goalWriteNamedAnalysis
    -- * Plotting
    , GnuplotOptions ( GnuplotOptions, maybeOutputDirectory
                     , whetherGeekie, whetherPNG, whetherLatex, whetherInteractive )
    , defaultGnuplotOptions
    , gnuplotCommandString
    , runGnuplot
    -- * Criterion
    , goalCriterionMain
    ) where


--- Imports ---


-- Unqualified --

import System.Directory
import GHC.Generics
import Data.String
import Control.Monad
import System.Process
import Paths_goal_core

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
goalExperimentPath :: Experiment -> IO FilePath
goalExperimentPath (Experiment prnm expnm) = do
    prpth <- goalProjectPath prnm
    return $ prpth ++ "/" ++ expnm

-- | Creates a project directory with the given name and returns its absolute path.
goalRawDataDirectory :: IO FilePath
goalRawDataDirectory =
    getXdgDirectory XdgData "goal/raw-data"

-- | Read a file in the given Goal project with the given file name.
goalReadFile
    :: Experiment
    -> String -- ^ File name
    -> IO String -- ^ File Contents
{-# INLINE goalReadFile #-}
goalReadFile expmnt flnm = do
    fpth <- goalExperimentPath expmnt
    readFile $ fpth ++ "/" ++ flnm

-- | Writes a file with the given filename.
goalWriteFile
    :: Experiment
    -> String -- ^ File name
    -> String -- ^ File contents
    -> IO ()
{-# INLINE goalWriteFile #-}
goalWriteFile expmnt flnm txt = do
    pth <- goalExperimentPath expmnt
    createDirectoryIfMissing True pth
    let pth' =  pth ++ "/" ++ flnm
    writeFile pth' txt


--- Analysis management ---


-- | Newtype for working with multiple datasets in a particular project.
newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance CSV.FromNamedRecord Dataset
instance CSV.ToNamedRecord Dataset
instance CSV.DefaultOrdered Dataset

data Experiment = Experiment
    { projectName :: String
    , experimentName :: String }

data SubExperiment = SubExperiment
    { analysisName :: String
    , datasetName :: String }

-- | Write the list of datasets.
goalWriteDatasetsCSV :: Experiment -> [String] -> IO ()
goalWriteDatasetsCSV expmnt dsts = do

    pth <- goalExperimentPath expmnt
    createDirectoryIfMissing True pth

    BS.writeFile (pth ++ "/datasets.csv") $ CSV.encodeDefaultOrderedByName $ Dataset <$> dsts

-- | Read the list of datasets (if it exists).
goalReadDatasetsCSV :: Experiment -> IO (Maybe [String])
goalReadDatasetsCSV expmnt = do

    pth <- goalExperimentPath expmnt
    let fpth = pth ++ "/datasets.csv"
    bl <- doesFileExist fpth
    if bl
       then do
           bstrm <- BS.readFile fpth
           let Right (_,as) = CSV.decodeByName bstrm
           return . Just $ dataset <$> V.toList as
       else return Nothing

-- | Run criterion with some defaults for goal. In particular save the results
-- in the benchmarks project.
goalCriterionMain :: String -> [C.Benchmark] -> IO ()
goalCriterionMain expnm bmrks = do

    exppth <- goalExperimentPath (Experiment "benchmarks" expnm)
    createDirectoryIfMissing True exppth

    let rptpth = exppth ++ "/" ++ "report.html"

    C.defaultMainWith (C.defaultConfig { C.reportFile = Just rptpth}) bmrks

-- | Write a single dataset to a file.
goalWriteDataset
    :: Experiment -- ^ Experiment
    -> String -- ^ Dataset name
    -> String -- ^ File contents
    -> IO ()
goalWriteDataset expmnt dstnm fl = do

    exppth <- goalExperimentPath expmnt

    let pth = exppth ++ "/data"
    createDirectoryIfMissing True pth

    let flpth = pth ++ "/" ++ dstnm ++ ".dat"
    writeFile flpth fl

-- | Read a single dataset from a file.
goalReadDataset
    :: Experiment -- ^ Experiment
    -> String -- ^ Dataset name
    -> IO String
goalReadDataset expmnt dstnm = do

    exppth <- goalExperimentPath expmnt
    let flpth = exppth ++ "/data/" ++ dstnm ++ ".dat"

    readFile flpth

analysisFilePath
    :: Bool -- ^ Create directories if missing
    -> Experiment
    -> Maybe SubExperiment
    -> IO String
analysisFilePath cbl expmnt (Just (SubExperiment ananm dstnm)) = do

    exppth <- goalExperimentPath expmnt
    let flpth0 = exppth ++ "/analysis/" ++ ananm

    when cbl $ createDirectoryIfMissing True flpth0
    return $ concat [flpth0,"/",dstnm,".csv"]

analysisFilePath cbl expmnt Nothing = do

    exppth <- goalExperimentPath expmnt

    when cbl $ createDirectoryIfMissing True exppth
    return . concat $ exppth ++ "/analysis.csv"

-- | Write the results of an analysis (in the form of a CSV) to the project
-- directory, using the goal organization structure. If there are multiple
-- datasets, analyses are named by dataset in a folder named after the analysis,
-- and otherwise the csv takes the name of the analysis itself.
--
-- The bool determines whether target files should be overwritten or appended
-- to. When appended, analysis blocks are seperated by two lines, allowing
-- gnuplot to distinguish them as different CSV blocks.
goalWriteAnalysis
    :: CSV.ToField x
    => Bool -- ^ Overwrite
    -> Experiment
    -> Maybe SubExperiment
    -> [[x]] -- ^ CSVs
    -> IO ()
goalWriteAnalysis True expmnt msbexp csvs = do

    flpth <- analysisFilePath True expmnt msbexp

    BS.writeFile flpth $ CSV.encode csvs

goalWriteAnalysis False expmnt msbexp csvs = do

    flpth <- analysisFilePath False expmnt msbexp

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encode csvs

-- | Write an analysis CSV based on named records.
goalWriteNamedAnalysis
    :: (CSV.DefaultOrdered csv, CSV.ToNamedRecord csv)
    => Bool
    -> Experiment
    -> Maybe SubExperiment
    -> [csv] -- ^ CSVs
    -> IO ()
goalWriteNamedAnalysis True expmnt msbexp csvs = do

    flpth <- analysisFilePath True expmnt msbexp

    BS.writeFile flpth $ CSV.encodeDefaultOrderedByName csvs

goalWriteNamedAnalysis False expmnt msbexp csvs = do

    flpth <- analysisFilePath False expmnt msbexp

    BS.appendFile flpth . BS.append (fromString "\r\n\r\n") $ CSV.encodeDefaultOrderedByName csvs


--- Plotting ---


gnuplotCommandString
    :: String -- ^ preload.gpi path
    -> String -- ^ .csv pth
    -> String -- ^ output file path
    -> Bool -- ^ Whether to Plot
    -> Bool -- ^ Whether to Latex
    -> Bool -- ^ Whether to interact
    -> String -- ^ .gpi path
    -> String -- ^ gnuplot command string
gnuplotCommandString pldpth anapth pltpthnm pbl lbl ibl gpipth =
     concat [ "gnuplot "
            , " -e \"preload='", pldpth
            , "'; csv='", anapth
            , "'; output_file_path='", pltpthnm
            , "'; do_plot=", if pbl then "1" else "0"
            , "; do_latex=", if lbl then "1" else "0"
            , "; do_interact=", if ibl then "1" else "0"
            ,"\" ", gpipth, if ibl then " -" else "" ]

data GnuplotOptions = GnuplotOptions
    { maybeOutputDirectory :: Maybe FilePath
    , whetherGeekie :: Bool
    , whetherPNG :: Bool
    , whetherLatex :: Bool
    , whetherInteractive :: Bool }

defaultGnuplotOptions :: GnuplotOptions
defaultGnuplotOptions = GnuplotOptions
    { maybeOutputDirectory = Nothing
    , whetherGeekie = False
    , whetherPNG = True
    , whetherLatex = False
    , whetherInteractive = False }

--runGNUPlotOpts (GNUPlotOpts prjnm expnm gpinm dstnm odr cbl rbl gbl pbl lbl ibl)

runGnuplot
    :: Experiment
    -> Maybe SubExperiment
    -> GnuplotOptions
    -> FilePath
    -> IO ()
runGnuplot expmnt msbexp (GnuplotOptions modr gbl pbl lbl ibl) gpipth = do

    pldpth <- getDataFileName "executables/preload.gpi"

    anapth <- analysisFilePath False expmnt msbexp

    exppth <- goalExperimentPath expmnt

    let pltnm = reverse . takeWhile (/= '/') . drop 4 $ reverse gpipth

    let pltpthnm =
            case msbexp of
              Just (SubExperiment ananm dstnm) ->
                   case modr of
                     Just odr -> concat [odr, "/", concat [dstnm,"-",pltnm]]
                     Nothing -> concat [exppth,"/plots/",ananm,"/",concat [dstnm,"-",pltnm]]
              Nothing ->
                  case modr of
                    Just odr -> concat [odr, "/", pltnm]
                    Nothing -> concat [exppth,"/",pltnm]

    when (null modr && not (null msbexp)) $
        createDirectoryIfMissing True (reverse . tail . dropWhile (/= '/') $ reverse pltpthnm)

    when gbl $ void . spawnCommand $ "geeqie " ++ exppth

    when (pbl || lbl || ibl) .
          callCommand $ gnuplotCommandString pldpth anapth pltpthnm pbl lbl ibl gpipth

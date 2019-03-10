{-# LANGUAGE OverloadedStrings,DeriveGeneric,FlexibleContexts #-}

-- | A set of Goal-specific functions for file, directory, and csv manipulation.
-- These functions use the XDG directory specification to save files in
-- appropriate directories.
module Goal.Core.Project
    (
    -- * Experiments
      Experiment (Experiment, projectName, experimentName)
    , Analysis (Analysis, analysisName, datasetName)
    -- * Import/Export
    , goalImport
    , goalExport
    , goalExportLines
    , goalExportNamed
    -- * Datasets
    , goalWriteDataset
    , goalReadDataset
    , goalWriteDatasetsCSV
    , goalReadDatasetsCSV
    -- * Plotting
    , GnuplotOptions ( GnuplotOptions, maybeOutputDirectory
                     , whetherGeekie, whetherPNG, whetherLatex, whetherInteractive, whetherAnimate )
    , defaultGnuplotOptions
    , gnuplotCommandString
    , runGnuplot
    -- * Criterion
    , goalCriterionMain
    -- * Util
    , goalCSVParser
    , goalCSVNamer
    , goalCSVOrder
    , deCamelCase
    ) where


--- Imports ---


-- Unqualified --

import System.Directory
import Control.Monad
import System.Process
import Paths_goal_core
import Data.Csv
import Data.Char
import GHC.Generics

-- Qualified --

import qualified Data.Vector as V
import qualified Data.ByteString.Lazy as BS
import qualified Data.ByteString as BSI
import qualified Criterion.Main as C
import qualified Criterion.Types as C


--- Experiments ---


data Experiment = Experiment
    { projectName :: String
    , experimentName :: String }

data Analysis = Analysis
    { analysisName :: String
    , datasetName :: String }


--- Import/Export ---


-- | Load a CSV file from the import directory of the given experiment.
goalImport
    :: FromRecord r
    => Experiment
    -> String
    -> IO (Maybe [r]) -- ^ CSVs
goalImport expmnt flnm = do

    expdr <- goalExperimentDirectory expmnt

    let flpth = concat [expdr,"/import/",flnm,".csv"]

    fbl <- doesFileExist flpth
    if fbl
       then do
           bstrm <- BS.readFile flpth
           let Right as = decode NoHeader bstrm
           return . Just $ V.toList as
        else return Nothing

-- | Write the results of an export (in the form of a CSV) to the project
-- directory, using the goal organization structure. If there are multiple
-- datasets, analyses are named by dataset in a folder named after the export,
-- and otherwise the csv takes the name of the export itself.
--
-- The bool determines whether target files should be overwritten or appended
-- to. When appended, export blocks are seperated by two lines, allowing
-- gnuplot to distinguish them as different CSV blocks.

goalExport
    :: ToRecord r
    => Bool -- ^ Overwrite
    -> Experiment
    -> Maybe Analysis
    -> [r] -- ^ CSVs
    -> IO ()
goalExport wbl expmnt msbexp csvs = do

    flpth <- exportFilePath True expmnt msbexp

    if wbl
       then BS.writeFile flpth $ encode csvs
       else BS.appendFile flpth . BS.append "\r\n\r\n" $ encode csvs

goalExportLines
    :: ToRecord r
    => Bool -- ^ Overwrite
    -> Experiment
    -> Maybe Analysis
    -> [[r]] -- ^ CSVs
    -> IO ()
goalExportLines wbl expmnt msbexp csvss = do

    flpth <- exportFilePath True expmnt msbexp

    let encs = BS.concat $ BS.tail . BS.tail . BS.append "\r\n" . encode <$> csvss

    if wbl
       then BS.writeFile flpth encs
       else BS.appendFile flpth $ BS.append "\r\n\r\n" encs

goalExportNamed
    :: (ToNamedRecord r, DefaultOrdered r)
    => Bool -- ^ Overwrite
    -> Experiment
    -> Maybe Analysis
    -> [r] -- ^ CSVs
    -> IO ()
goalExportNamed wbl expmnt msbexp csvs = do

    flpth <- exportFilePath True expmnt msbexp

    if wbl
       then BS.writeFile flpth $ encodeDefaultOrderedByName csvs
       else BS.appendFile flpth . BS.append "\r\n\r\n" $ encodeDefaultOrderedByName csvs


--- Dataset Management ---


-- | Newtype for working with multiple datasets in a particular project.
newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance FromNamedRecord Dataset where
    parseNamedRecord = genericParseNamedRecord deCamelCaseCSV

instance ToNamedRecord Dataset where
    toNamedRecord (Dataset dst) = namedRecord [ "Dataset" .= dst]

instance DefaultOrdered Dataset where
    headerOrder _ = header ["Dataset"]

-- | Write the list of datasets.
goalWriteDatasetsCSV :: Experiment -> [String] -> IO ()
goalWriteDatasetsCSV expmnt dsts = do

    pth <- goalExperimentDirectory expmnt
    createDirectoryIfMissing True pth

    BS.writeFile (pth ++ "/datasets.csv") . encodeDefaultOrderedByName $ Dataset <$> dsts

-- | Read the list of datasets (if it exists).
goalReadDatasetsCSV :: Experiment -> IO [String]
goalReadDatasetsCSV expmnt = do

    expdr <- goalExperimentDirectory expmnt
    let flpth = expdr ++ "/datasets.csv"

    fbl <- doesFileExist flpth

    if fbl
       then do
           bstrm <- BS.readFile flpth
           let Right (_,as) = decodeByName bstrm
           return $ dataset <$> V.toList as
       else error $ concat [ "datasets.csv does not exist. Has project "
                           , projectName expmnt, " been created?" ]

-- | Write a single dataset to a file.
goalWriteDataset
    :: Experiment -- ^ Experiment
    -> String -- ^ Dataset name
    -> String -- ^ File contents
    -> IO ()
goalWriteDataset expmnt dstnm fl = do

    expdr <- goalExperimentDirectory expmnt

    let pth = expdr ++ "/data"
    createDirectoryIfMissing True pth

    let flpth = pth ++ "/" ++ dstnm ++ ".dat"
    writeFile flpth fl

-- | Read a single dataset from a file.
goalReadDataset
    :: Experiment -- ^ Experiment
    -> String -- ^ Dataset name
    -> IO (Maybe String)
goalReadDataset expmnt dstnm = do

    expdr <- goalExperimentDirectory expmnt
    let flpth = expdr ++ "/data/" ++ dstnm ++ ".dat"

    fbl <- doesFileExist flpth

    if fbl then Just <$> readFile flpth
           else return Nothing


--- Criterion ---


-- | Run criterion with some defaults for goal. In particular save the results
-- in the benchmarks project.
goalCriterionMain :: String -> [C.Benchmark] -> IO ()
goalCriterionMain expnm bmrks = do

    expdr <- goalExperimentDirectory (Experiment "benchmarks" expnm)
    createDirectoryIfMissing True expdr

    let rptpth = expdr ++ "/" ++ "report.html"

    C.defaultMainWith (C.defaultConfig { C.reportFile = Just rptpth}) bmrks


--- Plotting ---


data GnuplotOptions = GnuplotOptions
    { maybeOutputDirectory :: Maybe FilePath
    , whetherGeekie :: Bool
    , whetherPNG :: Bool
    , whetherLatex :: Bool
    , whetherInteractive :: Bool
    , whetherAnimate :: Bool }

defaultGnuplotOptions :: GnuplotOptions
defaultGnuplotOptions = GnuplotOptions
    { maybeOutputDirectory = Nothing
    , whetherGeekie = False
    , whetherPNG = True
    , whetherLatex = False
    , whetherInteractive = False
    , whetherAnimate = False }

-- | Run gnuplot based on the given arguments.
runGnuplot
    :: Experiment
    -> Maybe Analysis
    -> GnuplotOptions
    -> FilePath
    -> IO ()
runGnuplot expmnt msbexp (GnuplotOptions modr gbl pbl lbl ibl abl) gpipth = do

    pldpth <- getDataFileName "executables/preload.gpi"

    anapth <- exportFilePath False expmnt msbexp

    expdr <- goalExperimentDirectory expmnt

    let pltnm = reverse . takeWhile (/= '/') . drop 4 $ reverse gpipth

    let pltpthnm =
            case msbexp of
              Just (Analysis ananm dstnm) ->
                   case modr of
                     Just odr -> concat [odr, "/", concat [dstnm,"-",pltnm]]
                     Nothing -> concat [expdr,"/plots/",ananm,"/",concat [dstnm,"-",pltnm]]
              Nothing ->
                  case modr of
                    Just odr -> concat [odr, "/", pltnm]
                    Nothing -> concat [expdr,"/",pltnm]

    when (null modr && not (null msbexp)) $
        createDirectoryIfMissing True (reverse . tail . dropWhile (/= '/') $ reverse pltpthnm)

    when gbl $ void . spawnCommand $ "geeqie " ++ expdr

    when (pbl || lbl || ibl || abl) .
          callCommand $ gnuplotCommandString pldpth anapth pltpthnm pbl lbl ibl abl gpipth


--- Util ---


deCamelCaseLoop :: String -> String
deCamelCaseLoop "" = ""
deCamelCaseLoop (c:wrds) =
    let (wrd,wrds') = span isLower wrds
     in (c:wrd) ++ ' ' : deCamelCaseLoop wrds'

deCamelCase :: String -> String
deCamelCase (c:wrds) = init $ deCamelCaseLoop (toUpper c : wrds)
deCamelCase "" = error "How is deCamelCase being run on an empty string?"

deCamelCaseCSV :: Options
deCamelCaseCSV = defaultOptions { fieldLabelModifier = deCamelCase }

goalCSVParser :: (Generic a, GFromNamedRecord (Rep a)) => NamedRecord -> Parser a
goalCSVParser = genericParseNamedRecord deCamelCaseCSV

goalCSVNamer
    :: (Generic a, GToRecord (Rep a) (BSI.ByteString, BSI.ByteString)) => a -> NamedRecord
goalCSVNamer = genericToNamedRecord deCamelCaseCSV

goalCSVOrder :: (Generic a, GToNamedRecordHeader (Rep a)) => a -> Header
goalCSVOrder = genericHeaderOrder deCamelCaseCSV

--- Internal ---


-- | Returns the xdg-based directory where projects are stored in Goal.
goalExperimentDirectory :: Experiment -> IO FilePath
goalExperimentDirectory (Experiment prnm expnm) = do
    xdgdr <- getXdgDirectory XdgData "goal/projects"
    return $ concat [xdgdr, "/",prnm,"/",expnm]

exportFilePath
    :: Bool -- ^ Create directories if missing
    -> Experiment
    -> Maybe Analysis
    -> IO String
exportFilePath cbl expmnt (Just (Analysis ananm dstnm)) = do

    expdr <- goalExperimentDirectory expmnt
    let flpth0 = expdr ++ "/export/" ++ ananm

    when cbl $ createDirectoryIfMissing True flpth0
    return $ concat [flpth0,"/",dstnm,".csv"]

exportFilePath cbl expmnt Nothing = do

    expdr <- goalExperimentDirectory expmnt

    when cbl $ createDirectoryIfMissing True expdr
    return $ expdr ++ "/export.csv"

gnuplotCommandString
    :: String -- ^ preload.gpi path
    -> String -- ^ .csv pth
    -> String -- ^ output file path
    -> Bool -- ^ Whether to Plot
    -> Bool -- ^ Whether to Latex
    -> Bool -- ^ Whether to interact
    -> Bool -- ^ Whether to animate
    -> String -- ^ .gpi path
    -> String -- ^ gnuplot command string
gnuplotCommandString pldpth anapth pltpthnm pbl lbl ibl abl gpipth =
     concat [ "gnuplot "
            , " -e \"preload='", pldpth
            , "'; csv='", anapth
            , "'; output_file_path='", pltpthnm
            , "'; do_plot=", if pbl then "1" else "0"
            , "; do_latex=", if lbl then "1" else "0"
            , "; do_animate=", if abl then "1" else "0"
            , "; do_interact=", if ibl then "1" else "0"
            ,"\" ", gpipth, if ibl then " -" else "" ]



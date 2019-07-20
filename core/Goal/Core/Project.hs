{-# LANGUAGE OverloadedStrings #-}
-- | This module provides functions for incorporating Goal into a
-- data-processing project. In particular, this module provides tools for
-- managing CSV files, and connecting them with gnuplot scripts for plotting.
module Goal.Core.Project
    (
    -- * CSV
      goalImport
    , goalExport
    , goalExportLines
    , goalExportNamed
    , goalExportNamedLines
    -- ** CSV Instances
    , goalCSVParser
    , goalCSVNamer
    , goalCSVOrder
    , deCamelCase
    -- * Util
    , runGnuplot
    , criterionMainWithReport
    ) where


--- Imports ---


-- Unqualified --

import System.Process
import System.Directory
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


--- Import/Export ---

runGnuplot
    :: FilePath
    -> FilePath
    -> IO ()
runGnuplot ldpth gpipth =
    callCommand $ concat [ "gnuplot ", " -e \"load_path='", ldpth, "'\" ",gpipth,".gpi" ]


-- | Load a CSV file from the import directory of the given experiment. The
-- @.csv@ extension is automatically added.
goalImport
    :: FromRecord r
    => FilePath
    -> IO (Either String [r]) -- ^ CSVs
goalImport flpth = do
    bstrm <- decode NoHeader <$> BS.readFile (flpth ++ ".csv")
    case bstrm of
      Right as -> return . Right $ V.toList as
      Left str -> return $ Left str

filePather :: FilePath -> FilePath -> IO FilePath
filePather ldpth flnm = do
    createDirectoryIfMissing True ldpth
    return $ concat [ldpth,"/",flnm,".csv"]

goalExport
    :: ToRecord r
    => FilePath
    -> FilePath
    -> [r] -- ^ CSVs
    -> IO ()
goalExport ldpth flnm csvs = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth $ encode csvs

goalExportLines
    :: ToRecord r
    => FilePath
    -> FilePath
    -> [[r]] -- ^ CSVss
    -> IO ()
goalExportLines ldpth flnm csvss = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth . BS.concat $ BS.tail . BS.tail . BS.append "\r\n" . encode <$> csvss

goalExportNamed
    :: (ToNamedRecord r, DefaultOrdered r)
    => FilePath
    -> FilePath
    -> [r] -- ^ CSVs
    -> IO ()
goalExportNamed ldpth flnm csvs = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth $ encodeDefaultOrderedByName csvs

goalExportNamedLines
    :: (ToNamedRecord r, DefaultOrdered r)
    => FilePath
    -> FilePath
    -> [[r]] -- ^ CSVss
    -> IO ()
goalExportNamedLines ldpth flnm csvss = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth . BS.concat $ BS.append "\r\n" . encodeDefaultOrderedByName <$> csvss


--- Criterion ---


-- | Run criterion with some defaults for goal. In particular save the results
-- in the benchmarks project.
criterionMainWithReport :: String -> [C.Benchmark] -> IO ()
criterionMainWithReport rprtnm =
    C.defaultMainWith (C.defaultConfig { C.reportFile = Just $ rprtnm ++ ".html"})

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



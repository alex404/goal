{-# LANGUAGE OverloadedStrings #-}
-- | This module provides functions for incorporating Goal into a
-- data-processing project. In particular, this module provides tools for
-- managing CSV files, and connecting them with gnuplot scripts for plotting.
-- CSV management is powered by @cassava@.
module Goal.Core.Project
    (
    -- * CSV
      goalImport
    , goalImportNamed
    , goalExport
    , goalExportLines
    , goalExportNamed
    , goalExportNamedLines
    -- ** CSV Instances
    , goalCSVParser
    , goalCSVNamer
    , goalCSVOrder
    -- * Util
    , runGnuplot
    , runGnuplotWithVariables
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


--- Experiments ---


--- Import/Export ---


-- | Runs @gnuplot@ on the given @.gpi@, passing it a @load_path@ variable to
-- help it find Goal-generated csvs.
runGnuplot
    :: FilePath -- ^ Gnuplot loadpath
    -> String -- ^ Gnuplot script
    -> IO ()
{-# INLINE runGnuplot #-}
runGnuplot ldpth gpipth =
    callCommand $ concat [ "gnuplot ", " -e \"load_path='", ldpth, "'\" ",gpipth,".gpi" ]

-- | Runs @gnuplot@ on the given @.gpi@, passing it a @load_path@ variable to
-- help it find Goal-generated csvs, and a list of variables.
runGnuplotWithVariables
    :: FilePath -- ^ Gnuplot loadpath
    -> String -- ^ Gnuplot script
    -> [(String,String)] -- ^ Arguments
    -> IO ()
{-# INLINE runGnuplotWithVariables #-}
runGnuplotWithVariables ldpth gpipth args =
    callCommand . concat $ [ "gnuplot ", " -e \"load_path='", ldpth, "'" ]
                         ++ (mapArgs <$> args) ++ [ "\" ",gpipth,".gpi" ]
        where mapArgs (nm,val) = concat ["; ",nm,"='",val,"'"]


-- | Load the given CSV file. The @.csv@ extension is automatically added.
goalImport
    :: FromRecord r
    => FilePath
    -> IO (Either String [r]) -- ^ CSVs
{-# INLINE goalImport #-}
goalImport flpth = do
    bstrm <- decode NoHeader <$> BS.readFile (flpth ++ ".csv")
    case bstrm of
      Right as -> return . Right $ V.toList as
      Left str -> return $ Left str

-- | Load the given CSV file with headers. The @.csv@ extension is automatically added.
goalImportNamed
    :: FromNamedRecord r
    => FilePath
    -> IO (Either String [r]) -- ^ CSVs
{-# INLINE goalImportNamed #-}
goalImportNamed flpth = do
    bstrm <- decodeByName <$> BS.readFile (flpth ++ ".csv")
    case bstrm of
      Right as -> return . Right . V.toList $ snd as
      Left str -> return $ Left str

filePather :: FilePath -> FilePath -> IO FilePath
{-# INLINE filePather #-}
filePather ldpth flnm = do
    createDirectoryIfMissing True ldpth
    return $ concat [ldpth,"/",flnm,".csv"]

-- | Export the given CSVs to a file in the given directory. The @.csv@
-- extension is automatically added to the file name.
goalExport
    :: ToRecord r
    => FilePath -- load_path
    -> String -- File Name
    -> [r] -- ^ CSVs
    -> IO ()
{-# INLINE goalExport #-}
goalExport ldpth flnm csvs = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth $ encode csvs

-- | Export the given list of CSVs to a file in the given directory, seperating
-- each set of CSVs by a single line. This causes gnuplot to the read CSV as a
-- collection of line segments. The @.csv@ extension is automatically added to
-- the file name.
goalExportLines
    :: ToRecord r
    => FilePath
    -> FilePath
    -> [[r]] -- ^ CSVss
    -> IO ()
{-# INLINE goalExportLines #-}
goalExportLines ldpth flnm csvss = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth . BS.concat $ BS.tail . BS.tail . BS.append "\r\n" . encode <$> csvss

-- | Export the named CSVs to a file in the given directory, adding a header to
-- the @.csv@ file.
goalExportNamed
    :: (ToNamedRecord r, DefaultOrdered r)
    => FilePath
    -> FilePath
    -> [r] -- ^ CSVs
    -> IO ()
{-# INLINE goalExportNamed #-}
goalExportNamed ldpth flnm csvs = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth $ encodeDefaultOrderedByName csvs

-- | Export the given list of named CSVs to a file, breaking it into a set of
-- line segments (with headers).
goalExportNamedLines
    :: (ToNamedRecord r, DefaultOrdered r)
    => FilePath
    -> FilePath
    -> [[r]] -- ^ CSVss
    -> IO ()
{-# INLINE goalExportNamedLines #-}
goalExportNamedLines ldpth flnm csvss = do
    flpth <- filePather ldpth flnm
    BS.writeFile flpth . BS.concat $ BS.append "\r\n" . encodeDefaultOrderedByName <$> csvss


--- Util ---


deCamelCaseLoop :: String -> String
{-# INLINE deCamelCaseLoop #-}
deCamelCaseLoop "" = ""
deCamelCaseLoop (c:wrds) =
    let (wrd,wrds') = span isLower wrds
     in (c:wrd) ++ ' ' : deCamelCaseLoop wrds'

deCamelCase :: String -> String
{-# INLINE deCamelCase #-}
deCamelCase (c:wrds) = init $ deCamelCaseLoop (toUpper c : wrds)
deCamelCase "" = error "How is deCamelCase being run on an empty string?"

deCamelCaseCSV :: Options
{-# INLINE deCamelCaseCSV #-}
deCamelCaseCSV = defaultOptions { fieldLabelModifier = deCamelCase }

-- | A generic @.csv@ parser which reorganizes a header name in camel case into
-- "human readable" text. Useful for instantiating 'FromNamedRecord'.
goalCSVParser :: (Generic a, GFromNamedRecord (Rep a)) => NamedRecord -> Parser a
{-# INLINE goalCSVParser #-}
goalCSVParser = genericParseNamedRecord deCamelCaseCSV

-- | A generic @.csv@ namer which reorganizes a header name in camel case into
-- "human readable" text. Useful for instantiating 'ToNamedRecord'.
goalCSVNamer
    :: (Generic a, GToRecord (Rep a) (BSI.ByteString, BSI.ByteString)) => a -> NamedRecord
{-# INLINE goalCSVNamer #-}
goalCSVNamer = genericToNamedRecord deCamelCaseCSV

-- | A generic @.csv@ order which reorganizes a header name in camel case into
-- "human readable" text. Useful for instantiating 'DefaultOrdered'.
goalCSVOrder :: (Generic a, GToNamedRecordHeader (Rep a)) => a -> Header
{-# INLINE goalCSVOrder #-}
goalCSVOrder = genericHeaderOrder deCamelCaseCSV



{-# LANGUAGE DeriveGeneric #-}

module Goal.Plot
    ( -- * Csvs
      Collection (Collection)
    , Dataset (Dataset)
    , getCollections
    , getDatasets
    -- * Reexports
    , module Options.Applicative
    ) where


--- Imports ---


import GHC.Generics
import Options.Applicative

-- Qualified --

import qualified Data.Csv as CSV
import qualified Data.ByteString.Lazy as BS
import qualified Data.Vector as V


--- Data management ---


newtype Collection = Collection { collection :: String } deriving (Read,Show,Generic)

newtype Dataset = Dataset { dataset :: String } deriving (Read,Show,Generic)

instance CSV.FromNamedRecord Collection
instance CSV.ToNamedRecord Collection
instance CSV.DefaultOrdered Collection

instance CSV.FromNamedRecord Dataset
instance CSV.ToNamedRecord Dataset
instance CSV.DefaultOrdered Dataset

getCollections :: IO [Collection]
getCollections = do

    bstrm <- BS.readFile "collections.csv"
    let Right (_,as) = CSV.decodeByName bstrm

    return $ V.toList as

getDatasets :: Collection -> IO [Dataset]
getDatasets (Collection clc) = do

    bstrm <- BS.readFile $ "projects/" ++ clc ++ "/datasets.csv"
    let Right (_,as) = CSV.decodeByName bstrm

    return $ V.toList as


----- CSV ---
--
--
--goalReadCSV ::  CSV.FromNamedRecord a => String -> String -> IO [a]
--goalReadCSV prj flnm = do
--
--    csvpth <- goalFilePath prj (flnm ++ ".csv")
--    bstrm <- BS.readFile csvpth
--
--    let Right (_,as) = CSV.decodeByName bstrm
--
--    return $ V.toList as
--
--goalWriteCSV ::  (CSV.DefaultOrdered a, CSV.ToNamedRecord a) => String -> String -> [a] -> IO ()
--goalWriteCSV prj flnm as = do
--    sbpth <- goalCreateProject prj
--    let fpth = sbpth ++ "/" ++ flnm ++ ".csv"
--    BS.writeFile fpth $ CSV.encodeDefaultOrderedByName as
--
--

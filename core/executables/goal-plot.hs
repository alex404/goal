--- Imports ---

import Goal.Core

import System.Process

--- Opt Parse ---


--goalWriteDataset ::  String -> Dataset -> String -> IO ()
--goalWriteDataset prj (Dataset flnm) dat = do
--    sbpth <- goalCreateProject (prj ++ "/data")
--    let fpth = sbpth ++ "/" ++ flnm ++ ".dat"
--    writeFile fpth dat
--
--goalReadDataset ::  String -> Dataset -> IO String
--goalReadDataset prj (Dataset flnm) = do
--    sbpth <- goalCreateProject (prj ++ "/data")
--    let fpth = sbpth ++ "/" ++ flnm ++ ".dat"
--    readFile fpth

data GNUPlotOpts = GNUPlotOpts String String String Bool

gnuPlotOpts :: Parser GNUPlotOpts
gnuPlotOpts = GNUPlotOpts
      <$> strArgument (help "GNUPlot script")
      <*> strOption
            ( long "collection" <> short 'c' <> help "Which data collection to plot" <> value "")
      <*> strOption
            ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
      <*> switch ( long "interactive" <> short 'i' <> help "Start interactive session" )

gnuPlotCommand :: String -> Collection -> Dataset -> Bool -> String
gnuPlotCommand gpi (Collection clc) (Dataset dst) ibl =
    concat [ "gnuplot ", " -e \"collection='", clc, "'; dataset='", dst
           , "'; interact=" , if ibl then "1" else "0","\" ", gpi, if ibl then " -" else ""]


runGNUPlotOpts :: GNUPlotOpts -> IO ()
runGNUPlotOpts (GNUPlotOpts gpipth clcstr dststr ibl) = do

    when (ibl && clcstr == "" && dststr /= "") . void $ fail "Unspecified project for dataset"
    when (ibl && dststr == "") . void $ fail "Dataset must be specified for interactive mode"

    when ibl $ do
        void . callCommand $ gnuPlotCommand gpipth (Collection clcstr) (Dataset dststr) ibl
        fail ""

    clcs <- if clcstr == ""
                 then getCollections
                 else return [Collection clcstr]
    dstss <- if dststr == ""
               then mapM getDatasets clcs
               else return [[Dataset dststr]]

    sequence_ $ do
        (clc,dsts) <- zip clcs dstss
        dst <- dsts
        return . spawnCommand $ gnuPlotCommand gpipth clc dst ibl




--- Main ---


main :: IO ()
main = runGNUPlotOpts =<< execParser opts
    where opts = info (gnuPlotOpts <**> helper)
              (fullDesc <> progDesc "Run gnuplot on collections of datasets")

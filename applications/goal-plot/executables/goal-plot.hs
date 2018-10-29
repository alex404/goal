--- Imports ---

import Goal.Plot
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

data GNUPlotOpts = GNUPlotOpts String String String String Bool Bool Bool

gnuPlotOpts :: Parser GNUPlotOpts
gnuPlotOpts = GNUPlotOpts
      <$> strArgument (help "GNUPlot script")
      <*> strArgument ( help "Which data collection to plot")
      <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
      <*> strOption ( long "odir" <> short 'o' <> help "Force this output directory" <> value "")
      <*> switch ( long "interactive" <> short 'p' <> help "Generate plots" )
      <*> switch ( long "latex" <> short 'l' <> help "Produce latex output")
      <*> switch ( long "interactive" <> short 'i' <> help "Start interactive session" )

gnuPlotCommand :: String -> Collection -> Dataset -> String -> Bool -> Bool -> Bool -> String
gnuPlotCommand gpi (Collection clc) (Dataset dst) ostr pbl lbl ibl =
     concat [ "gnuplot "
            , " -e \"collection='", clc
            , "'; dataset='", dst
            , if ostr == "" then "" else "'; odir='" ++ ostr
            , "'; plot=", if pbl then "1" else "0"
            , "; latex=", if lbl then "1" else "0"
            , "; interact=", if ibl then "1" else "0"
            ,"\" ", gpi, if ibl then " -" else ""]


runGNUPlotOpts :: GNUPlotOpts -> IO ()
runGNUPlotOpts (GNUPlotOpts gpipth clcstr dststr ostr pbl lbl ibl) = do

    let clc = Collection clcstr

    dsts <- if dststr == ""
               then getDatasets clc
               else return [Dataset dststr]

    sequence_ $ do
        dst <- dsts
        return . callCommand $ gnuPlotCommand gpipth clc dst ostr pbl lbl ibl




--- Main ---


main :: IO ()
main = runGNUPlotOpts =<< execParser opts
    where opts = info (gnuPlotOpts <**> helper)
              (fullDesc <> progDesc "Run gnuplot on collections of datasets")

--- Imports ---

import Goal.Core
import System.Process
import Paths_goal_core

--- Opt Parse ---


data GNUPlotOpts = GNUPlotOpts String String String String String Bool Bool Bool Bool Bool Bool

gnuPlotOpts :: Parser GNUPlotOpts
gnuPlotOpts = GNUPlotOpts
      <$> strArgument ( help "Which project to work with")
      <*> strArgument ( help "Which experiment to work with")
      <*> strArgument (help "GNUPlot script name" <> value "")
      <*> strOption ( long "dataset" <> short 'd' <> help "Which dataset to plot" <> value "")
      <*> strOption ( long "odir" <> short 'o' <> help "Force this output directory" <> value "")
      <*> switch ( long "copy" <> short 'c' <> help "Copy local gpi to project dir (and then run comands)")
      <*> switch ( long "read" <> short 'r' <> help "Read local gpi")
      <*> switch ( long "geeqie" <> short 'g' <> help "Launch geeqie in the experiment's plots directory" )
      <*> switch ( long "plot" <> short 'p' <> help "Generate plots" )
      <*> switch ( long "latex" <> short 'l' <> help "Produce latex output")
      <*> switch ( long "interactive" <> short 'i' <> help "Start interactive session" )

gnuPlotCommand :: String -> String -> Maybe String -> String -> Bool -> Bool -> Bool -> String -> String
gnuPlotCommand gpipth pth mdstnm odr pbl lbl ibl pldpth =
    let dststr = maybe "" ("'; dataset='" ++) mdstnm
     in concat [ "gnuplot "
               , " -e \"path='", pth
               , dststr
               , "'; preload='", pldpth
               , if odr == "" then "" else "'; odir='" ++ odr
               , "'; plot=", if pbl then "1" else "0"
               , "; latex=", if lbl then "1" else "0"
               , "; interact=", if ibl then "1" else "0"
               ,"\" ", gpipth, if ibl then " -" else "" ]

runGNUPlotOpts :: GNUPlotOpts -> IO ()
runGNUPlotOpts (GNUPlotOpts prjnm expnm gpinm dstnm odr cbl rbl gbl pbl lbl ibl) = do

    mdsts <- if dstnm == ""
               then goalReadDatasetsCSV prjnm expnm
               else return $ Just [Dataset dstnm]

    prjpth <- goalProjectPath prjnm
    exppth <- goalExperimentPath prjnm expnm
    pldpth <- getDataFileName "preload.gpi"

    let gpifl = gpinm ++ ".gpi"
        gpipth = prjpth ++ "/gpis"
        gpiflpth = gpipth ++ "/" ++ gpifl

    when gbl $ void . spawnCommand $ concat ["geeqie ", exppth, "/plots"]

    when cbl $ do
        createDirectoryIfMissing True gpipth
        copyFile gpifl gpiflpth

    let gpiflpth' = if rbl then gpifl else gpiflpth

    let dsts = maybe [Nothing] (\dsts0 -> [Just dst | Dataset dst <- dsts0]) mdsts

    when (pbl || lbl || ibl) . sequence_ $ do
        dst <- dsts
        return . callCommand $ gnuPlotCommand gpiflpth' exppth dst odr pbl lbl ibl pldpth



--- Main ---


main :: IO ()
main = runGNUPlotOpts =<< execParser opts
    where opts = info (gnuPlotOpts <**> helper)
              (fullDesc <> progDesc "Run gnuplot on collections of datasets")

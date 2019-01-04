--- Imports ---

import Goal.Core


--- Opt Parse ---


data GnuplotArgs = GnuplotArgs String String String String String String Bool Bool Bool Bool

gnuplotArgs :: Parser GnuplotArgs
gnuplotArgs = GnuplotArgs
      <$> strArgument ( help "Project Name")
      <*> strArgument ( help "Experiment Name")
      <*> strArgument (help "gnuplot .gpi script" <> value "")
      <*> strOption ( long "analysis" <> short 'a' <> help "(maybe) Analysis Name" <> value "" )
      <*> strOption ( long "dataset" <> short 'd' <> help "(maybe) Dataset Name" <> value "" )
      <*> strOption ( long "odir" <> short 'o' <> help "Force this output directory" <> value "" )
      <*> switch ( long "geeqie" <> short 'g' <> help "Launch geeqie in the experiment's root directory" )
      <*> switch ( long "plot" <> short 'p' <> help "Generate png output" )
      <*> switch ( long "latex" <> short 'l' <> help "Produce latex output" )
      <*> switch ( long "interactive" <> short 'i' <> help "Start interactive session" )

runGnuplotArgs :: GnuplotArgs -> IO ()
runGnuplotArgs (GnuplotArgs prjnm expnm _ _ _ _ True _ _ _) = do

    let expmnt = Experiment prjnm expnm

    runGnuplot expmnt Nothing (GnuplotOptions Nothing True False False False) ""

runGnuplotArgs (GnuplotArgs prjnm expnm gpipth [] _ odr False pbl lbl ibl) = do

    let expmnt = Experiment prjnm expnm

        modr = if null odr then Nothing else Just odr

    runGnuplot expmnt Nothing (GnuplotOptions modr False pbl lbl ibl) gpipth

runGnuplotArgs (GnuplotArgs prjnm expnm gpipth ananm dstnm odr False pbl lbl ibl) = do

    let expmnt = Experiment prjnm expnm

        modr = if null odr then Nothing else Just odr

    mdstnms <- if null dstnm
               then goalReadDatasetsCSV expmnt
               else return $ Just [dstnm]

    when (pbl || lbl || ibl) $ do

        case mdstnms of
          Just dstnms -> sequence_ $ do
                  dstnm' <- dstnms
                  let msbexp = Just $ SubExperiment ananm dstnm'
                  return $ runGnuplot expmnt msbexp (GnuplotOptions modr False pbl lbl ibl) gpipth
          Nothing -> error "No dataset name given and cannot find datasets.csv"



--- Main ---


main :: IO ()
main = runGnuplotArgs =<< execParser opts
    where opts = info (gnuplotArgs <**> helper)
              (fullDesc <> progDesc "Run gnuplot on collections of datasets")

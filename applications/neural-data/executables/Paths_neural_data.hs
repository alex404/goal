module Paths_neural_data where

getDataFileName :: FilePath -> IO FilePath
getDataFileName fp = return $ "../data/" ++ fp

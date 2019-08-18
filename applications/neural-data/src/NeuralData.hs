{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    FlexibleContexts,
    TypeOperators
    #-}

module NeuralData
    ( -- * Variables
      mnx
    , mxx
    , loadRoot
    , loadPath
    -- * CLI
    , ExperimentOpts (ExperimentOpts)
    , experimentOpts
    ) where


--- Imports ---


-- Goal --

import Goal.Core

--- Variables ---

loadRoot :: String
loadRoot = "/home/alex404/development/goal/applications/neural-data/data"

loadPath :: String -> String -> String
loadPath expmnt dst = concat [loadRoot, "/", expmnt, "/", dst]

mnx,mxx :: Double
mnx = 0
mxx = 2*pi


--- Inference ---


data ExperimentOpts = ExperimentOpts String String

experimentOpts :: Parser ExperimentOpts
experimentOpts = ExperimentOpts
    <$> strArgument
        ( help "Which data collection to analyze"
        <> metavar "EXPERIMENT" )
    <*> strArgument
        ( help "Which dataset to plot (if no argument is given than all datasets in the project are analyzed)"
        <> metavar "DATASET" )

{-# LANGUAGE
    GADTs,
    ScopedTypeVariables,
    DataKinds,
    FlexibleContexts,
    TypeOperators
    #-}

module NeuralData
    ( -- * Variables
      prjnm
    , mnx
    , mxx
    -- * IO
    , getNeuralData
    , strengthenNeuralData
    -- * General
    , meanSDInliers
    , dataCovariances
    -- * Subsampling
    , generateIndices
    , subSampleResponses
    -- * CLI
    , ExperimentOpts (ExperimentOpts)
    , experimentOpts
    , readDatasets
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S
import qualified Goal.Core.Vector.Generic as G

-- Qualified --

import qualified Data.Map as M
import qualified Data.List as L


--- Types ---


--- Variables ---

prjnm :: String
prjnm = "neural-data"

mnx,mxx :: Double
mnx = 0
mxx = 2*pi


--- Inference ---


getNeuralData :: Read s => String -> String -> IO (NatNumber,[([Int], s)])
getNeuralData expnm dst = read . fromJust <$> goalReadDataset (Experiment prjnm expnm) dst

strengthenNeuralData :: (KnownNat k, Read s) => [([Int], s)] -> [(Response k, s)]
strengthenNeuralData xss =
    let (ks,ss) = unzip xss
     in zip (fromJust . S.fromList <$> ks) ss


--- Analysis ---


dataCovariances :: (Ord x, KnownNat k) => [(Response k,x)] -> [(Source # MultivariateNormal k,x)]
{-# INLINE dataCovariances #-}
dataCovariances zxs =
    let zxss = L.groupBy (on (==) snd) $ L.sortOn snd zxs
     in [ (mle $ G.convert . G.map realToFrac <$> zs, head xs) | (zs,xs) <- unzip <$> zxss ]


meanSDInliers :: [Double] -> (Double,Double)
meanSDInliers xs =
    let (mu,vr) = estimateMeanVariance xs
        xs' = filter (\x -> square (x-mu) < 4*vr) xs
        (mu',vr') = estimateMeanVariance xs'
     in (mu',sqrt vr')

subSampleResponses
    :: (Ord s, KnownNat k, KnownNat m)
    => M.Map s [Response (k+m)]
    -> Response k
    -> M.Map s [Response k]
subSampleResponses zxmp idxs =
     map (`S.backpermute` idxs) <$> zxmp

generateIndices
    :: forall k m r v . (KnownNat k, KnownNat m, G.VectorClass v Int)
    => Proxy (k + m)
    -> Random r (G.Vector v k Int)
generateIndices _ = do
    let idxs :: G.Vector v (k + m) Int
        idxs = G.generate finiteInt
    subsampleVector idxs

data ExperimentOpts = ExperimentOpts String String

experimentOpts :: Parser ExperimentOpts
experimentOpts = ExperimentOpts
    <$> strArgument
        ( help "Which data collection to analyze"
        <> metavar "EXPERIMENT" )
    <*> strOption
        ( short 'd'
        <> long "dataset"
        <> help "Which dataset to plot (if no argument is given than all datasets in the project are analyzed)"
        <> value "" )

readDatasets :: ExperimentOpts -> IO [String]
readDatasets (ExperimentOpts expnm dstarg) =
    if null dstarg
       then goalReadDatasetsCSV (Experiment prjnm expnm)
       else return [dstarg]

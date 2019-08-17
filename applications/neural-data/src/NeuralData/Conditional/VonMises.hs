{-# LANGUAGE TypeApplications #-}
module NeuralData.Conditional.VonMises
    ( -- * IO
    parametersPath
    , datasetPath
    , readMixtureLikelihood
    , writeMixtureLikelihood
    , strengthenMixtureLikelihood
    , writeDataset
    , readDataset
    ) where


--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S

import NeuralData


--- IO ---


parametersPath :: String -> String -> FilePath
parametersPath expmnt dst = loadPath expmnt dst ++ "/parameters.dat"

datasetPath :: String -> String -> FilePath
datasetPath expmnt dst = loadPath expmnt dst ++ "/dataset.dat"

writeDataset
    :: forall k . KnownNat k
    => String -- ^ Experiment
    -> String -- ^ Dataset
    -> Sample (Neurons k, VonMises)  -- ^ Mixture Likelihood
    -> IO ()
writeDataset expmnt dst zxs = do
    createDirectoryIfMissing True $ loadPath expmnt dst
    let flpth = datasetPath expmnt dst
        (zs,xs) = unzip zxs
        zs' = S.toList <$> zs
    writeFile flpth . show $ (natValInt (Proxy @ k), zip zs' xs)

readDataset
    :: String -- ^ Experiment
    -> String -- ^ Dataset
    -> IO (NatNumber,[([Int],Double)])
readDataset expmnt dst =
    let flpth = datasetPath expmnt dst
     in read <$> readFile flpth

writeMixtureLikelihood
    :: forall k n . (KnownNat k, KnownNat n)
    => String -- ^ Experiment
    -> String -- ^ Dataset
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises -- ^ Mixture Likelihood
    -> IO ()
writeMixtureLikelihood expmnt dst lkl = do
    createDirectoryIfMissing True $ loadPath expmnt dst
    let k = natValInt (Proxy @ k)
        n = natValInt (Proxy @ n)
        flpth = parametersPath expmnt dst
    writeFile flpth $ show (k,n,listCoordinates lkl)

readMixtureLikelihood
    :: String -- ^ Experiment
    -> String -- ^ Dataset
    -> IO (NatNumber,NatNumber,[Double])
readMixtureLikelihood expmnt dst =
    let flpth = parametersPath expmnt dst
     in read <$> readFile flpth

strengthenMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => [Double]
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
strengthenMixtureLikelihood xs = Point . fromJust $ S.fromList xs

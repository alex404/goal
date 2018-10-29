{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}
{-# LANGUAGE FlexibleContexts,GADTs,ScopedTypeVariables,DataKinds,TypeOperators #-}

--- Imports ---


-- Goal --

import NeuralData

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---


-- General --

mnx,mxx :: Double
mnx = 0
mxx = 2*pi

nstms :: Int
nstms = 8

nsmps :: Int
nsmps = 100

stms :: [Double]
stms = tail $ range mnx mxx (nstms + 1)

gn0 :: Double
gn0 = 10

-- Convolutional --

cmus :: KnownNat k => S.Vector k Double
cmus = S.init $ S.range mnx mxx

ckp :: Double
ckp = 1

csps :: KnownNat k => S.Vector k (Source # VonMises)
csps = S.map (Point . flip S.doubleton ckp) cmus

clkl :: KnownNat k => Mean #> Natural # R k Poisson <* VonMises
clkl = vonMisesPopulationEncoder False (Left gn0) csps

-- Random --

rmus :: KnownNat k => Random r (S.Vector k Double)
rmus = S.replicateM $ uniformR (mnx,mxx)

rkps :: KnownNat k => Random r (S.Vector k Double)
rkps = S.replicateM $ uniformR (0.5,1.5)

rsps :: KnownNat k => Random r (S.Vector k (Source # VonMises))
rsps = S.zipWith (\x y -> Point $ S.doubleton x y) <$> rmus <*> rkps

rlklr :: KnownNat k => Random r (Mean #> Natural # R k Poisson <* VonMises)
rlklr = vonMisesPopulationEncoder True (Left gn0) <$> rsps

normalizeLikelihood
    :: KnownNat k
    => Mean #> Natural # Neurons k <* VonMises
    -> Mean #> Natural # Neurons k <* VonMises
normalizeLikelihood lkl0 =
    let (nz,nzx) = splitAffine lkl0
        bnd = 0.0001
        eps = -0.005
        xsmps = range mnx mxx 100
        cauchify = last . take 10000 . cauchySequence euclideanDistance bnd
        rho0 = average $ potential <$> lkl0 >$>* xsmps
        diff = populationCodeRectificationDifferential rho0 zero xsmps nzx
        nz' = cauchify $ vanillaGradientSequence diff eps defaultAdamPursuit nz
     in joinAffine nz' nzx

combineStimuli :: [[Response k]] -> [([Int],Double)]
combineStimuli zss =
    concat $ zipWith (\zs x -> zip (toList <$> zs) $ repeat x) zss stms

-- IO --

syntheticPath :: Int -> String
syntheticPath k = "synthetic-" ++ show k ++ "-neurons"


--- Main ---


newtype SyntheticOpts = SyntheticOpts Int

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> option auto (long "kneurons" <> help "number of neurons" <> short 'k' <> value 50)


synthesizeData :: forall k . KnownNat k => Proxy k -> IO ()
synthesizeData prxk = do

    let k = natValInt prxk

    rlkl <- realize rlklr

    let nrlkl = normalizeLikelihood rlkl

    (czss :: [[Response k]]) <- realize (mapM (sample nsmps) $ clkl >$>* stms)
    (rzss :: [[Response k]]) <- realize (mapM (sample nsmps) $ rlkl >$>* stms)
    (nrzss :: [[Response k]]) <- realize (mapM (sample nsmps) $ nrlkl >$>* stms)


    let czxs,rzxs,nrzxs :: [([Int], Double)]
        czxs = combineStimuli czss
        rzxs = combineStimuli rzss
        nrzxs = combineStimuli nrzss

    let dsts@[cnvdst,rnddst,nrmdst] = Dataset <$> ["convolutional","random","random-normalized"]

    goalWriteDataset prjnm (syntheticPath k) cnvdst $ show czxs
    goalWriteDataset prjnm (syntheticPath k) rnddst $ show rzxs
    goalWriteDataset prjnm (syntheticPath k) nrmdst $ show nrzxs

    goalWriteDatasetsCSV prjnm (syntheticPath k) dsts

runOpts :: SyntheticOpts -> IO ()
runOpts (SyntheticOpts k) = withNat k synthesizeData

--- Main ---


main :: IO ()
main = do

    let opts = info (syntheticOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr = "Generate synthetic data"

    runOpts =<< execParser opts

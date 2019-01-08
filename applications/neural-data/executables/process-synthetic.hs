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

nstms :: Int
nstms = 8

nsmps :: Int
nsmps = 100

stms :: [Double]
stms = tail $ range mnx mxx (nstms + 1)

wghts :: KnownNat n => Natural # Categorical Int n
wghts = zero

fromConditionalOneMixture
    :: Mean #> Natural # MixtureGLM (Neurons k) Int 0 VonMises
    -> Mean #> Natural # Neurons k <* VonMises
fromConditionalOneMixture = breakPoint



-- Convolutional --


cmus :: KnownNat k => S.Vector k Double
cmus = S.init $ S.range mnx mxx

ckp :: Double
ckp = 1

csps :: KnownNat k => S.Vector k (Source # VonMises)
csps = S.map (Point . flip S.doubleton ckp) cmus

cgnss :: forall n k . (KnownNat n, KnownNat k) => S.Vector n (Source # Neurons k)
cgnss = S.generateP generator
    where generator :: (KnownNat j) => Proxy j -> (Source # Neurons k)
          generator prxj = Point . S.replicate $ 10 + 5 * fromIntegral (natValInt prxj)

clkl
    :: (KnownNat k, KnownNat n)
    => Proxy n -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
clkl _ = vonMisesMixturePopulationEncoder True wghts (S.map toNatural cgnss) (S.map toNatural csps)

-- Random --

rmus :: KnownNat k => Random r (S.Vector k Double)
rmus = S.replicateM $ uniformR (mnx,mxx)

rkps :: KnownNat k => Random r (S.Vector k Double)
rkps = S.replicateM $ uniformR (0.5,1.5)

rsps :: KnownNat k => Random r (S.Vector k (Source # VonMises))
rsps = S.zipWith (\x y -> Point $ S.doubleton x y) <$> rmus <*> rkps

rgnss :: (KnownNat n, KnownNat k) => Random r (S.Vector n (Source # Neurons k))
rgnss = S.replicateM . fmap Point . S.replicateM $ uniformR (10,20)

rlklr
    :: (KnownNat k, KnownNat n)
    => Random r (Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises)
rlklr = do
    gnss <- rgnss
    vonMisesMixturePopulationEncoder True wghts (S.map toNatural gnss) . S.map toNatural <$> rsps

normalizeMixtureLikelihood
    :: (KnownNat k, KnownNat n)
    => Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
    -> Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
normalizeMixtureLikelihood lkl0 =
    let (nzk,nzx) = splitBottomSubLinear lkl0
        bnd = 0.0001
        eps = -0.001
        xsmps = range mnx mxx 100
        cauchify = last . take 10000 . cauchySequence euclideanDistance bnd
        rho0 = average $ potential <$> lkl0 >$>* xsmps
        diff = conditionalHarmoniumRectificationDifferential rho0 zero xsmps nzx
        nzk' = cauchify $ vanillaGradientSequence diff eps defaultAdamPursuit nzk
     in joinBottomSubLinear nzk' nzx

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

syntheticExperiment :: Int -> Experiment
syntheticExperiment k = Experiment prjnm $ "synthetic-" ++ show k ++ "k"

trueSyntheticExperiment :: Int -> Experiment
trueSyntheticExperiment k = Experiment prjnm $ "true-" ++ (experimentName $ syntheticExperiment k)

syntheticMixtureExperiment :: Int -> Int -> Experiment
syntheticMixtureExperiment k n = Experiment prjnm $ "synthetic-" ++ show k ++ "k-" ++ show n ++ "n"

trueSyntheticMixtureExperiment :: Int -> Int -> Experiment
trueSyntheticMixtureExperiment k n = Experiment prjnm $ "true-" ++ (experimentName $ syntheticMixtureExperiment k n)

--- Main ---


synthesizeData :: forall k . KnownNat k => Proxy k -> Double -> Double -> Double -> Double -> IO ()
synthesizeData prxk lgmu lgsd ltmu ltsd = do

    rlkl0 <- realize rlklr
    let rlkl = fromConditionalOneMixture rlkl0

    let nrlkl :: Mean #> Natural # Neurons k <* VonMises
        nrlkl = normalizeLikelihood rlkl

    let clkln = fromConditionalOneMixture $ clkl Proxy

    (czss :: [[Response k]]) <- realize (mapM (sample nsmps) $ clkln >$>* stms)
    (rzss :: [[Response k]]) <- realize (mapM (sample nsmps) $ rlkl >$>* stms)
    (nrzss :: [[Response k]]) <- realize (mapM (sample nsmps) $ nrlkl >$>* stms)


    let czxs,rzxs,nrzxs :: [([Int], Double)]
        czxs = combineStimuli czss
        rzxs = combineStimuli rzss
        nrzxs = combineStimuli nrzss

    let dsts@[cnvdst,rnddst,nrmdst] = ["convolutional","random","random-normalized"]

    let k = natValInt prxk

    goalWriteDataset (syntheticExperiment k) cnvdst $ show (k,czxs)
    goalWriteDataset (syntheticExperiment k) rnddst $ show (k,rzxs)
    goalWriteDataset (syntheticExperiment k) nrmdst $ show (k,nrzxs)

    goalWriteDatasetsCSV (syntheticExperiment k) dsts

    goalWriteDataset (trueSyntheticExperiment k) cnvdst $ show (k,listCoordinates clkln)
    goalWriteDataset (trueSyntheticExperiment k) rnddst $ show (k,listCoordinates rlkl)
    goalWriteDataset (trueSyntheticExperiment k) nrmdst $ show (k,listCoordinates nrlkl)

    goalWriteDatasetsCSV (trueSyntheticExperiment k) dsts

synthesizeMixtureData :: forall k n . (KnownNat k, KnownNat n) => Proxy k -> Proxy n -> IO ()
synthesizeMixtureData prxk prxn = do

    (rlkl :: Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises) <- realize rlklr

    let nrlkl :: Mean #> Natural # MixtureGLM (Neurons k) Int n VonMises
        nrlkl = normalizeMixtureLikelihood rlkl

    (rzss :: [[Response k]]) <- realize (mapM (fmap (map hHead) . sample nsmps) $ rlkl >$>* stms)
    (nrzss :: [[Response k]]) <- realize (mapM (fmap (map hHead) . sample nsmps) $ nrlkl >$>* stms)


    let rzxs,nrzxs :: [([Int], Double)]
        rzxs = combineStimuli rzss
        nrzxs = combineStimuli nrzss

    let dsts@[rnddst,nrmdst] = ["random","random-normalized"]

    let k = natValInt prxk
        n = natValInt prxn

    goalWriteDataset (syntheticMixtureExperiment k n) rnddst $ show (k,rzxs)
    goalWriteDataset (syntheticMixtureExperiment k n) nrmdst $ show (k,nrzxs)

    goalWriteDatasetsCSV (syntheticMixtureExperiment k n) dsts

    goalWriteDataset (trueSyntheticMixtureExperiment k n) rnddst $ show (k,n,listCoordinates rlkl)
    goalWriteDataset (trueSyntheticMixtureExperiment k n) nrmdst $ show (k,n,listCoordinates nrlkl)

    goalWriteDatasetsCSV (trueSyntheticMixtureExperiment k n) dsts


--- Opt Parse ---


data SyntheticOpts = SyntheticOpts NatNumber NatNumber Double Double Double Double

syntheticOpts :: Parser SyntheticOpts
syntheticOpts = SyntheticOpts
    <$> option auto
        ( short 'k'
        <> long "k-neurons"
        <> help "Number of neurons in the model population."
        <> showDefault
        <> value 20 )
    <*> option auto
        ( short 'm'
        <> long "n-mixers"
        <> help "Number of mixers to use to model neural correlations."
        <> showDefault
        <> value 0 )
    <*> option auto
        ( short 'g'
        <> long "mean-log-gain"
        <> help "The mean of the log-gains."
        <> showDefault
        <> value 2 )
    <*> option auto
        ( short 'G'
        <> long "sd-log-gain"
        <> help "The standard deviation of the log-gains."
        <> showDefault
        <> value 1 )
    <*> option auto
        ( short 't'
        <> long "mean-log-tuning"
        <> help "The mean of the log tuning-widths."
        <> showDefault
        <> value 0 )
    <*> option auto
        ( short 'T'
        <> long "sd-log-tuning"
        <> help "The standard deviation of the log tuning-width."
        <> showDefault
        <> value 0.5 )



runOpts :: SyntheticOpts -> IO ()
runOpts (SyntheticOpts k n lgmu lgsd ltmu ltsd)
  | n == 0 = case someNatVal k of
              SomeNat prxk -> synthesizeData prxk lgmu lgsd ltmu ltsd
  | otherwise = case someNatVal (n-1) of
                  SomeNat prxn -> case someNatVal k of
                    SomeNat prxk -> synthesizeMixtureData prxk prxn

--- Main ---


main :: IO ()
main = do

    let opts = info (syntheticOpts <**> helper) (fullDesc <> progDesc prgstr)
        prgstr =
            "Generate synthetic data from a model population of\
            \stimulus-dependent, Poisson neurons. Model parameters can be\
            \specified with command line arguments."

    runOpts =<< execParser opts

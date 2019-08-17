{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise #-}

{-# LANGUAGE
    FlexibleContexts,
    TypeFamilies,
    TypeOperators,
    TypeApplications,
    ScopedTypeVariables,
    DataKinds
    #-}


--- Imports ---


import NeuralData
import NeuralData.Conditional.VonMises

import Goal.Core
import Goal.Geometry
import Goal.Probability

import qualified Goal.Core.Vector.Storable as S


--- Globals ---


nbns :: Int
nbns = 10

xsmps :: [Double]
xsmps = init $ range mnx mxx 101


--- Functions ---


fisherInformation
    :: (KnownNat n, KnownNat k)
    => Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
    -> Double
    -> Double
fisherInformation mlkl x =
    let cvr = snd . splitMultivariateNormal . mixturePopulationCovariance $ mlkl >.>* x
        nzx = snd $ splitConditionalDeepHarmonium mlkl
        dsx = Point $ S.fromTuple (-sin x, cos x)
        dmu = coordinates $ nzx >.> dsx
     in S.dotProduct dmu $ S.matrixVectorMultiply cvr dmu

averageLogFisherInformation
    :: (KnownNat n, KnownNat k)
    => Sample VonMises
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
    -> Double
averageLogFisherInformation xs mlkl =
    let avg = average [ 0.5 * log (2*pi*exp 1 / fisherInformation mlkl x) | x <- xs ]
     in log (2*pi) - avg

numericalLogPosterior
    :: (KnownNat n, KnownNat k)
    => Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
    -> SamplePoint (Neurons k)
    -> SamplePoint VonMises
    -> Double
numericalLogPosterior mlkl z =
    let logupst x = logMixtureDensity (mlkl >.>* x) z
        logprt = logIntegralExp 1e-6 logupst 0 (2*pi) xsmps
     in \x -> logupst x - logprt

numericalPosteriorDivergence
    :: (KnownNat n, KnownNat k)
    => Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
    -> SamplePoint (Neurons k)
    -> Double
numericalPosteriorDivergence mlkl z =
    let logdns = numericalLogPosterior mlkl z
        dns = exp . logdns
        klf x = dns x * logdns x
     in (+ log (2*pi)) . fst $ integrate 1e-6 klf 0 (2*pi)

numericalMutualInformation
    :: (KnownNat n, KnownNat k)
    => Sample VonMises
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
    -> Random r Double
numericalMutualInformation xs mlkl = do
    zs <- mapM (fmap hHead . samplePoint) (mlkl >$>* xs)
    return . average $ numericalPosteriorDivergence mlkl <$> zs

subsampleMixtureLikelihood
    :: (KnownNat k, KnownNat n, KnownNat m)
    => Natural #> ConditionalMixture (Neurons (k+m)) n Tensor VonMises
    -> S.Vector k Int
    -> Natural #> ConditionalMixture (Neurons k) n Tensor VonMises
subsampleMixtureLikelihood mlkl idxs =
    let (dhrm,tns) = splitConditionalDeepHarmonium mlkl
        tns' = fromMatrix . S.fromRows . flip S.backpermute idxs . S.toRows $ toMatrix tns
        (cmpnts,wghts) = splitMixture dhrm
        cmpnts' = S.map (Point . flip S.backpermute idxs . coordinates) cmpnts
     in joinConditionalDeepHarmonium (joinMixture cmpnts' wghts) tns'

subsampleInformation
    :: forall k m r . (KnownNat k, KnownNat m)
    => Natural #> ConditionalMixture (Neurons k) m Tensor VonMises -- ^ Likelihood
    -> NatNumber -- ^ Number of subsamples
    -> NatNumber -- ^ Subsample size
    -> Random r [(Double,Double)]
subsampleInformation mlkl npop k' = do
    let k = natVal (Proxy @ k)
        l = k - k'
    case someNatVal k' of
        SomeNat (Proxy :: Proxy k') -> case someNatVal l of
            SomeNat (Proxy :: Proxy l) -> do
                let prxk = Proxy @ k
                    prxkl = Proxy @ (k' + l)
                case sameNat prxk prxkl of
                  Just Refl -> do
                      idxss <- replicateM (fromIntegral npop) $ generateIndices prxkl
                      let mlkls' :: [Natural #> ConditionalMixture (Neurons k') m Tensor VonMises]
                          mlkls' = subsampleMixtureLikelihood mlkl <$> idxss
                          avgfshs = averageLogFisherInformation xsmps <$> mlkls'
                      --avginfs <- mapM (numericalMutualInformation xsmps) mlkls'
                          avginfs = [ average $ fisherInformation mlkl' <$> xsmps | mlkl' <- mlkls' ]
                      return $ zip avgfshs avginfs
                  Nothing -> error "What?"


informationStats :: NatNumber -> [Double] -> (NatNumber,Double,Double,Double)
informationStats k' infs = (k', minimum infs, average infs, maximum infs)

expInformationStats :: NatNumber -> [Double] -> (NatNumber,Double,Double,Double)
expInformationStats k' infs0 =
    let infs = exp <$> infs0
     in (k', minimum infs, average infs, maximum infs)

--- CLI ---


data OrderOpts = OrderOpts NatNumber NatNumber

validationOpts :: Parser OrderOpts
validationOpts = OrderOpts
    <$> option auto
        ( short 'n'
        <> long "n-population"
        <> help "Number of populations to generate per size."
        <> showDefault
        <> value 10 )
    <*> option auto
        ( short 's'
        <> long "count"
        <> help "Number of population sizes to test."
        <> showDefault
        <> value 10 )


data AllOpts = AllOpts ExperimentOpts OrderOpts

allOpts :: Parser AllOpts
allOpts = AllOpts <$> experimentOpts <*> validationOpts

runOpts :: AllOpts -> IO ()
runOpts ( AllOpts (ExperimentOpts expmnt dst0) (OrderOpts npop nstp) ) = do

    let dst = dst0 ++ "/standard-hybrid-em"
        ldpth = loadPath expmnt dst

    (k,m,cs) <- readMixtureLikelihood expmnt dst

    putStrLn "\nNumber of Neurons:"
    print k

    putStrLn "\nNumber of Mixture Components:"
    print m

    let stps = (div k nstp *) <$> [1..nstp]

    case someNatVal k of
        SomeNat (Proxy :: Proxy k) -> case someNatVal m
            of SomeNat (Proxy :: Proxy m) -> do

                let trulkl :: Natural #> ConditionalMixture (Neurons k) m Tensor VonMises
                    trulkl = strengthenMixtureLikelihood cs

                fimistts <- forM stps $ \k' -> do

                    putStrLn "\nSubSample Size:"
                    print k'

                    fimis <- realize $ subsampleInformation trulkl npop k'

                    let (fis,mis) = unzip fimis

                    putStrLn "\nAverage Log-Fisher Informations:"
                    print fis

                    putStrLn "\nMutual Informations:"
                    print mis

                    return (expInformationStats k' fis, informationStats k' mis)

                let (fistts,mistts) = unzip fimistts

                goalExport ldpth "fisher-informations" fistts
                goalExport ldpth "mutual-informations" mistts

                runGnuplot ldpth "information-order"

--- Main ---


main :: IO ()
main = do

    let prgstr = "Measuring the order of correlations in large populations."
        hdrstr = prgstr
        opts = info (allOpts <**> helper) (fullDesc <> progDesc prgstr <> header hdrstr)
    runOpts =<< execParser opts

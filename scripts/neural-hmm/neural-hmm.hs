#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import Control.Concurrent.Async

import qualified Goal.Core.Vector.Storable as S
import qualified Data.Map as M
import qualified Data.List as L

-- Learning

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)

type Neurons n = Replicated n Poisson
type CoMNeurons n = Replicated n CoMPoisson
type Response n = S.Vector n Int

strengthenObservations :: KnownNat n => [[[Int]]] -> [[Response n]]
strengthenObservations nsss =
    [ fromJust . S.fromList <$> nss | nss <- nsss ]

processData :: Int -> IO ()
processData ssn = do

    let flnm = "data/gratings-session-" ++ show ssn

    putStrLn $ "\nProcessing Data: " ++ flnm

    dmpstr <- readFile $ flnm ++ ".map"
    kmpstr <- readFile $ flnm ++ ".keys"

    let dmp :: M.Map Int [[[Int]]]
        dmp = read dmpstr

    let kmp :: M.Map Int (Double,Double,Double)
        kmp = read kmpstr

    let stms = M.keys dmp

    let zss00 = dmp M.! 0

    putStrLn $ "\nSample Size: " ++ show (length zss00)
    putStrLn $ "Sequence Length: " ++ (show . length $ head zss00)
    putStrLn $ "Number of Neurons: " ++ (show . length . head $ head zss00)

    forM_ stms $ \stm -> do

        let zss0 = dmp M.! stm

        let n = length . head $ head zss0

        case someNatVal (fromIntegral n) of

          SomeNat (Proxy :: Proxy n) -> do

              let zss :: [[Response n]]
                  zss = strengthenObservations zss0

              let (cnt,cori,sori) = kmp M.! stm

              (mu,sd,cmu,csd,mu',sd') <- fitData zss

              putStrLn $!! concat
                  [ "\nContrast: ", show cnt
                  , ", Centre Ori: ", show cori
                  , ", Surround Ori: ", show sori
                  ,"\nSwitching CV Log-Likelihood: ", show mu, " ± ", show sd
                  ,"\nCoM CV Log-Likelihood: ", show cmu, " ± ", show csd
                  ,"\nStatic CV Log-Likelihood: ", show mu', " ± ", show sd']

printer x = do
    print =<< x
    x

fitData
    :: forall n . KnownNat n
    => [[Response n]]
    -> IO (Double,Double,Double,Double,Double,Double)
fitData zss0 = do
    zss1 <- realize $ shuffleList zss0
    let zss = map ((!!3) . S.toList) <$> zss1
        tvzss = kFold 10 zss

    putStrLn "Sample Trains:"
    print $ length zss

    llss <- forConcurrently tvzss $ \(tzss,vzss) -> do

        let nits = 50
            eps = 0.05
            nstps = 100

        ltnt0 :: Natural # LatentProcess Tensor Tensor Poisson (Categorical 1) Poisson (Categorical 1)
            <- realize $ uniformInitialize (-0.01,0.01)
        cltnt0 :: Natural # Affine Tensor Poisson (LatentProcess Tensor Tensor Poisson (Categorical 0) Poisson (Categorical 0)) (Categorical 7)
            <- realize $ uniformInitialize (-0.01,0.01)
            --cltnt0 =
                --let (ehrm,trns) = split ltnt0
                --    (pstr,nz) = split $ transposeHarmonium ehrm
                --    mapper z = Point $ coordinates z S.++ S.singleton (-1)
                --    ncz = mapReplicatedPoint mapper nz
                -- in join (transposeHarmonium $ join pstr ncz) trns
        let gp = defaultAdamPursuit

        let psn :: Natural # Poisson
            psn = toNatural . averageSufficientStatistic $ concat tzss

        let psns :: [Natural # Poisson]
            psns = toNatural . averageSufficientStatistic <$> L.transpose tzss


        let ltnts = take nits $ iterate (latentProcessExpectationMaximizationAscent eps nstps gp tzss) ltnt0
        let ltnts' = take nits $ iterate (latentProcessExpectationMaximization tzss) ltnt0
            --vll,vll' :: Double
            --vll = maximum $ average . (`logObservableDensities` vzss) <$> ltnts
            vll = (*8) . average . logDensities psn $ concat vzss
            cvll = sum . map average $ zipWith logDensities psns $ L.transpose vzss
            vll' = maximum $ average . (`logObservableDensities` vzss) <$> ltnts'
        --print $ toMean psn
        --print $ toMean <$> psns
        --let vll' = 0

            cltnts = take nits $ iterate
                (timeDependentLatentProcessExpectationMaximizationAscent eps nstps gp tzss) cltnt0
            --cvll = maximum $ average . (\cltnt -> timeDependentConjugatedSmoothingLogDensity cltnt <$> vzss) <$> cltnts

        return (vll,cvll,vll')

    let (lls,clls,lls') = unzip3 llss
        (mu,vr) = estimateMeanVariance lls
        (cmu,cvr) = estimateMeanVariance clls
        (mu',vr') = estimateMeanVariance lls'

    return (mu,sqrt vr,cmu,sqrt cvr,mu',sqrt vr')


--- Main ---


main :: IO ()
main = processData 1

    --let em (prr',trns',emsn') = stateSpaceExpectationMaximization prr' trns' emsn' zss

    --    hmms = take 500 $ iterate em (prr0,trns0,emsn0)

    --let lls (prr',trns',emsn') =
    --        average $ conjugatedFilteringLogDensity trns' emsn' prr' <$> zss

    --mapM_ (print . lls) hmms

    --putStrLn "\nModels:"
    --putStrLn "\nInitial:"
    --printHMM $ head hmms
    --putStrLn "\nLearned:"
    --printHMM $ last hmms

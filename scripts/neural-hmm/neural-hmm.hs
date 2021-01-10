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

-- Learning

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)

type Neurons n = Replicated n Poisson
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

    forConcurrently_ stms $ \stm -> do

        let zss0 = dmp M.! stm

        let n = length . head $ head zss0

        case someNatVal (fromIntegral n) of

          SomeNat (Proxy :: Proxy n) -> do

              let zss :: [[Response n]]
                  zss = strengthenObservations zss0

              let (cnt,cori,sori) = kmp M.! stm

              (mu,sd,mu',sd') <- realize $ fitData zss

              putStrLn $!! concat
                  [ "\nContrast: ", show cnt
                  , ", Centre Ori: ", show cori
                  , ", Surround Ori: ", show sori
                  ,"\nSwitching CV Log-Likelihood: ", show mu, " ± ", show sd
                  ,"\nStatic CV Log-Likelihood: ", show mu', " ± ", show sd']

fitData
    :: forall n r . KnownNat n
    => [[Response n]]
    -> Random r (Double,Double,Double,Double)
fitData zss0 = do
    zss <- shuffleList zss0
    let tvzss = kFold 5 zss

    llss <- forM tvzss $ \(tzss,vzss) -> do

        ltnt0 :: Natural # LatentProcess Tensor Tensor (Neurons n) (Categorical 1)
            <- uniformInitialize (-0.01,0.01)

        ltnt0' :: Natural # LatentProcess Tensor Tensor (Neurons n) (Categorical 0)
            <- uniformInitialize (-0.01,0.01)

        let ltnts = take 50 $ iterate (latentProcessExpectationMaximization tzss) ltnt0
            ltnts' = take 50 $ iterate (latentProcessExpectationMaximization tzss) ltnt0'
            vll,vll' :: Double
            vll = maximum $ average . (`logObservableDensities` vzss) <$> ltnts
            vll' = maximum $ average . (`logObservableDensities` vzss) <$> ltnts'

        return (vll,vll')

    let (lls,lls') = unzip llss
        (mu,vr) = estimateMeanVariance lls
        (mu',vr') = estimateMeanVariance lls'

    return (mu,sqrt vr,mu',sqrt vr')


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

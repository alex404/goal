#! /usr/bin/env stack
-- stack runghc --

{-# LANGUAGE ScopedTypeVariables,TypeOperators,TypeFamilies,FlexibleContexts,DataKinds #-}

--- Imports ---


-- Goal --

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S
import qualified Data.Map as M


-- Learning

alg :: (Double,GradientPursuit,Int)
alg = (0.05,defaultAdamPursuit,100)

type Response n = S.Vector n Int

strengthenNeuralData :: KnownNat n => [[[Int]]] -> [[[Response n]]]
strengthenNeuralData nsss =
     map (fromJust . S.fromList <$> ks) ss

--printHMM
--    :: ( KnownNat x, KnownNat z )
--    => ( Natural # Categorical x
--       , Natural # Categorical x <* Categorical x
--       , Natural # Categorical x <* Categorical z )
--    -> IO ()
--printHMM (prr',trns',emsn') = do
--    putStrLn "Prior: "
--    print . S.toList $ categoricalWeights prr'
--    putStrLn "Transitions: "
--    mapM_ print $ S.toList . categoricalWeights <$> trns' >$>* xspc
--    putStrLn "Emissions: "
--    mapM_ print $ S.toList . categoricalWeights <$> emsn' >$>* xspc

fitData :: Int -> IO ()
fitData ssn = do

    let flnm = "gratings-session-" ++ show ssn ++ "-map"

    putStrLn $ "\nProcessing File: " ++ flnm

    trns0 :: Natural # Affine Tensor (Categorical 1) (Categorical 1)
        <- realize $ uniformInitialize (-1,1)
    emsn0 :: Natural # Affine Tensor (Categorical 2) (Categorical 1)
        <- realize $ uniformInitialize (-1,1)
    prr0 :: Natural # Categorical 2 <- realize $ uniformInitialize (-1,1)

    dmpstr <- readFile flnm

    let dmp :: M.Map Int [[[Int]]]
        dmp = read dmpstr



    print "foo"


--- Main ---


main :: IO ()
main = undefined

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

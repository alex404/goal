#! stack runghc

{-# LANGUAGE ScopedTypeVariables,TypeApplications,DeriveGeneric,GADTs,FlexibleContexts,TypeOperators,DataKinds,KindSignatures,TypeFamilies,NoStarIsType,UndecidableInstances,TypeSynonymInstances,FlexibleInstances,MultiParamTypeClasses
   #-}

{-# OPTIONS_GHC -fplugin=GHC.TypeLits.KnownNat.Solver -fplugin=GHC.TypeLits.Normalise -fconstraint-solver-iterations=10 #-}

--- Imports ---

import Goal.Core
import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Goal.Core.Vector.Storable as S
import qualified Data.List as L


--- Globals ---


type N = 440
type M = 20
type K = 10

ldpth,csvpth :: FilePath
ldpth = "data"
csvpth = ldpth ++ "/kernels_selected.csv"

-- Training --

parseCSV :: String -> [S.Vector N Double]
parseCSV csvstr = do
    csv <- lines csvstr
    let xs = read $ '[' : csv ++ "]"
    return . fromJust $ S.fromList xs

-- | This is a selfcontained implementation of PCA-EM from Bishop, against
-- which I can benchmark the EF-based version.

ns :: Source # Normal
ns = fromTuple (0,1)

randomMOG :: Random (Natural # Mixture (FullNormal M) K)
randomMOG = do
    muss <- S.replicateM $ initialize ns
    let nrms = S.map (flip join (snd . split $ standardNormal)) muss
    return $ joinNaturalMixture nrms 0

-- Optimization

nepchs :: Int
nepchs = 100

eps :: Double
eps = 3e-4

nstps :: Int
nstps = 500

mogEM
    :: ( Square Natural f (MVNMean 20), Transition Source Mean (MultivariateNormal f 20)
       , Bilinear Source f (MVNMean 20) (MVNMean 20)
       , ExponentialFamily (LocationShape (MVNMean M) (f (MVNMean M) (MVNMean M))) )
    => [S.Vector M Double]
    -> Natural # Mixture (MultivariateNormal f M) K
    -> Natural #  Mixture (MultivariateNormal f M) K
mogEM zs igh0 =
    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps

hmogEM
    :: [S.Vector N Double]
    -> Natural # IsotropicHMoG N M K
    -> Natural # IsotropicHMoG N M K
hmogEM zs igh0 =
    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps

dmogEM
    :: [S.Vector N Double]
    -> Natural # DiagonalHMoG N M K
    -> Natural # DiagonalHMoG N M K
dmogEM zs igh0 =
    expectationMaximizationAscent eps defaultAdamPursuit zs igh0 !! nstps

--conjugatedPotentialTest hrm = do
--    let (lkl,nw) = split hrm
--        (rho0,rprms) = conjugationParameters lkl
--     in 0*potentialTest (nw + rprms) + rho0
--
--potentialTest p =
--    let (nmu,nsgma) = split p
--        (insgma,dtmnt,sgn) = S.inverseLogDeterminant . negate $ 2 * nsgma
--     in -0.25 * S.dotProduct nmu (S.matrixVectorMultiply insgma nmu) -0.5 * dtmnt



--- Main ---


main :: IO ()
main = do

    csvstr <- readFile csvpth

    let smps0 = parseCSV csvstr
        n = length smps0
    print $ "Number of samples: " ++ show n
    smps <- realize $ shuffleList smps0
    let tsmps = smps
        vsmps = tsmps
        --(tsmps,vsmps) = splitAt (round $ (0.8 :: Double) * fromIntegral n) smps

    let nvx :: Natural # IsotropicNormal N
        nvx = toNatural $ averageSufficientStatistic tsmps

    let nvx' :: Natural # DiagonalNormal N
        nvx' = toNatural $ averageSufficientStatistic tsmps

    (lds0 :: Natural # Tensor (MVNMean N) (MVNMean M))
        <- realize $ uniformInitialize (-1,1)

    let npca0 :: Natural # PrincipleComponentAnalysis N M
        npca0 = join nvx lds0
        nfa0 :: Natural # FactorAnalysis N M
        nfa0 = join nvx' lds0

    let emnpcas = take nepchs $ iterate (linearModelExpectationMaximization tsmps) npca0
        emnfas = take nepchs $ iterate (linearModelExpectationMaximization tsmps) nfa0

    let nlls = zip
            (logLikelihood vsmps . linearModelObservableDistribution <$> emnpcas)
            (logLikelihood vsmps . linearModelObservableDistribution <$> emnfas)

    putStrLn "Natural PCA/FA LL Ascent:"
    mapM_ print nlls

    let npca1 = last emnpcas
        nfa1 = last emnfas
        prjctn1,prjctn1' :: Natural # Affine Tensor (MVNMean M) (FullNormal M) (MVNMean N)
        prjctn1 = fst . split . transposeHarmonium . joinConjugatedHarmonium npca1 $ standardNormal
        prjctn1' = fst . split . transposeHarmonium . joinConjugatedHarmonium nfa1 $ standardNormal

        tprjcts1 = coordinates . fst . split . toSource <$> prjctn1 >$>* tsmps
        vprjcts1 = coordinates . fst . split . toSource <$> prjctn1 >$>* vsmps
        tprjcts1' = coordinates . fst . split . toSource <$> prjctn1' >$>* tsmps
        vprjcts1' = coordinates . fst . split . toSource <$> prjctn1' >$>* vsmps

    mog0 <- realize $ randomMOG

    let mogs = take nepchs $ iterate (mogEM tprjcts1) mog0
        mogs' = take nepchs $ iterate (mogEM tprjcts1') mog0

    let moglls = logLikelihood vprjcts1 <$> mogs
        moglls' = logLikelihood vprjcts1' <$> mogs'
        hmoglls = logLikelihood tsmps . joinConjugatedHarmonium npca1 <$> mogs
        hmoglls' = logLikelihood tsmps . joinConjugatedHarmonium nfa1 <$> mogs'
    putStrLn "Mog vs HMog LL Ascent:"
    mapM_ print $ L.zip4 moglls hmoglls moglls' hmoglls'
    --mapM_ print $ L.zip moglls hmoglls

    let hmog = joinConjugatedHarmonium nfa1 $ head mogs'
    --print $ logLikelihood tsmps hmog
    --print . any (isInfinite) $ listCoordinates hmog
    --let prms :: (Double,Natural # FullNormal M)
    --    prms = conjugationParameters . fst $ split hmog

    --    (thts,tht3) = split . fst $ split hmog
    --    (tht1,tht2) = split thts
    --    (itht20,lndt,_) = inverseLogDeterminant . negate $ 2 .> tht2
    --    itht2 = -2 .> itht20
    --    tht21 = itht2 >.> tht1
    --    rho0 = -0.25 * (tht1 <.> tht21) -0.5 * lndt
    --    rho1 = -0.5 .> (transpose tht3 >.> tht21)
    --    rho2 = -0.25 .> changeOfBasis tht3 itht2

    --print lndt


    let mog1 = last mogs
        mog1' = last mogs'
      --hmog1 = joinConjugatedHarmonium npca1 mog1
      --hmog1' = joinConjugatedHarmonium nfa1 mog1'


    let emhmogs = take nepchs . iterate (hmogEM tsmps) $ joinConjugatedHarmonium npca1 mog1
    let emhmogs' = take nepchs . iterate (dmogEM tsmps) $ joinConjugatedHarmonium nfa1 mog1'
        hmoglls1 = logLikelihood tsmps <$> emhmogs
        hmoglls1' = logLikelihood tsmps <$> emhmogs'

    putStrLn "True HMog LL Ascent:"
    mapM_ print $ zip hmoglls1 hmoglls1'

    putStrLn "Information Gains:"
    print $ (last hmoglls1 - last hmoglls,last hmoglls1' - last hmoglls')


    --let hmog2 = last emhmogs
    --    (npca2,mog2) = splitConjugatedHarmonium hmog2
    --let hmog2' = last emhmogs'
    --    (nfa2,mog2') = splitConjugatedHarmonium hmog2'

    --    prjctn2,prjctn2' :: Natural # Affine Tensor (MVNMean M) (MultivariateNormal M) (MVNMean N)
    --    prjctn2 = fst . split . transposeHarmonium . joinConjugatedHarmonium npca2 . transition
    --        $ joinMultivariateNormal 0 S.matrixIdentity
    --    prjctn2' = fst . split . transposeHarmonium . joinConjugatedHarmonium nfa2 . transition
    --        $ joinMultivariateNormal 0 S.matrixIdentity

    --    tprjcts2 = fst . splitMultivariateNormal . toSource <$> prjctn2 >$>* tsmps
    --    tprjcts2' = fst . splitMultivariateNormal . toSource <$> prjctn2' >$>* tsmps


    --let prjctss = [tprjcts1,tprjcts1',tprjcts2,tprjcts2']
    --    clstrmp1 = fst $ split mog1
    --    clstrmp2 = fst $ split mog2
    --    clstrmp1' = fst $ split mog1'
    --    clstrmp2' = fst $ split mog2'
    --    clstrmps = [clstrmp1,clstrmp1',clstrmp2,clstrmp2']
    --    nms = ["standard-pca","standard-fa","joint-pca","joint-fa"]

    --forM_ (zip3 prjctss clstrmps nms) $ \(prjcts,clstrmp,nm) -> do

    --    let catmp :: M.Map Int [S.Vector M Double]
    --        catmp = M.fromListWith (++) $ zip tcats $ (:[]) <$> prjcts
    --        kys = M.keys catmp

    --    sequence_ $ do
    --        ky <- kys
    --        let pnts = S.toList <$> catmp M.! ky
    --            mvn = toSource $ clstrmp >.>* ky
    --            xys = bivariateNormalConfidenceEllipse 100 1 mvn
    --        return $ do
    --            goalExport ldpth ("projection-cat-" ++ show ky) $ pnts
    --            goalExport ldpth ("clustering-cat-" ++ show ky) $ xys
    --    runGnuplotWithVariables ldpth "projection" [("nm",nm)]

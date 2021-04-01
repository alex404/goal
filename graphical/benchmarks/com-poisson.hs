{-# LANGUAGE DataKinds,ScopedTypeVariables,FlexibleContexts,TypeOperators,TypeFamilies #-}

import Goal.Geometry
import Goal.Probability
import Goal.Graphical

import qualified Criterion.Main as C


--- Globals ---


type N = 101
type K = 10

type PoissonMixture = Mixture (Replicated N Poisson) K
type CoMPoissonMixture = AffineMixture (Replicated N Poisson) (Replicated N CoMPoisson) K

mixtureStep
    :: (LegendreExponentialFamily f, Legendre f, ExpectationMaximization f)
    => [Observation f] -> Natural # f -> Natural # f
mixtureStep zs mxmdl =
    expectationMaximizationAscent 2e-3 defaultAdamPursuit zs mxmdl !! 20

comDeviations
    :: Double
    -> Double
    -> Mean # CoMPoisson
    -> Natural # CoMPoisson
    -> (Double,Double,Double)
comDeviations err lgprt0 mcm0 ncm =
    let lgprt = comPoissonLogPartitionSum err ncm
        mcm = comPoissonMeans err ncm
        [muerr,nuerr] = listCoordinates $ mcm0 - mcm
     in (lgprt0 - lgprt, muerr,nuerr)

comSDs
    :: [Double]
    -> [Mean # CoMPoisson]
    -> [Natural # CoMPoisson]
    -> Double
    -> (Double,Double,Double)
comSDs lgprts0 mcm0s ncms err =
    let (lgprterrs,muerrs,nuerrs) = unzip3 $ zipWith3 (comDeviations err) lgprts0 mcm0s ncms
        lgprtsd = snd $ estimateMeanVariance lgprterrs
        musd = snd $ estimateMeanVariance muerrs
        nusd = snd $ estimateMeanVariance nuerrs
     in (lgprtsd,musd,nusd)


--- Main ---


main :: IO ()
main = do

    mx0 :: Natural # PoissonMixture
        <- realize $ uniformInitialize (-2,2)

    zxs <- realize $ sample 100 mx0
    let zs = fst <$> zxs

    let (nyx,ny) = split $ transposeHarmonium mx0
        cmx0 :: Natural # CoMPoissonMixture
        cmx0 = transposeHarmonium . join nyx
            $ mapReplicatedPoint (`join` (-1.5)) ny

    let cms :: [Source # CoMPoisson]
        cms = do
            mu <- [0.5,2,20]
            nu <- [0.3,1,10]
            return $ fromTuple (mu,nu)
        ncms = toNatural <$> cms

    let err0 = 1e-20
    let lgprt0s = comPoissonLogPartitionSum err0 <$> ncms
        mcm0s = comPoissonMeans err0 <$> ncms

    let errs = [1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16]
        comFun = comSDs lgprt0s mcm0s ncms

    sequence_ $ do
        err <- errs
        return $ do
            putStrLn $ concat ["Error: ", show err]
            putStrLn $ concat ["Deviations: ", show $ comFun err, "\n" ]


    --- Criterion ---

    putStrLn "\n"


    C.defaultMain
       [ C.bench "com-error-1e-2" $ C.nf comFun 1e-2
       , C.bench "com-error-1e-4" $ C.nf comFun 1e-4
       , C.bench "com-error-1e-6" $ C.nf comFun 1e-6
       , C.bench "com-error-1e-8" $ C.nf comFun 1e-8
       , C.bench "com-error-1e-10" $ C.nf comFun 1e-10
       , C.bench "com-error-1e-12" $ C.nf comFun 1e-12
       , C.bench "com-error-1e-14" $ C.nf comFun 1e-14
       , C.bench "expectation-step" $ C.nf (expectationStep zs) mx0
       , C.bench "com-expectation-step" $ C.nf (expectationStep zs) cmx0
       , C.bench "transition" $ C.nf toMean mx0
       , C.bench "com-transition" $ C.nf toMean cmx0
       , C.bench "mixture-step" $ C.nf (mixtureStep zs) mx0
       , C.bench "com-mixture-step" $ C.nf (mixtureStep zs) cmx0 ]
